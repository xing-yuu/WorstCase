#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <wcso/mesh/MeshIO.h>
#include <wcso/optim/BilevelSolver.h>
#include <wcso/remesh/MeshImprove.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace fs = std::filesystem;
using nlohmann::json;
namespace {
json read_json(const fs::path &path) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open " + path.string());
    json j;
    input >> j;
    return j;
}
void write_json(const fs::path &path, const json &j) {
    const fs::path temporary = path.string() + ".tmp";
    {
        std::ofstream out(temporary);
        if (!out)
            throw std::runtime_error("Cannot write " + temporary.string());
        out << std::setw(2) << j << '\n';
        out.flush();
        if (!out)
            throw std::runtime_error("Failed writing " + temporary.string());
    }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), path.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        throw std::runtime_error("Cannot atomically replace " + path.string());
#else
    fs::rename(temporary, path);
#endif
}
template <class Matrix> json matrix_json(const Matrix &m) {
    json result = json::array();
    for (int i = 0; i < m.rows(); ++i) {
        json row = json::array();
        for (int j = 0; j < m.cols(); ++j)
            row.push_back(m(i, j));
        result.push_back(row);
    }
    return result;
}
template <class Scalar>
Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_from(const json &j, int columns) {
    if (!j.is_array())
        throw std::runtime_error("Malformed checkpoint matrix");
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m(static_cast<int>(j.size()), columns);
    for (int i = 0; i < m.rows(); ++i) {
        if (!j[i].is_array() || j[i].size() != static_cast<size_t>(columns))
            throw std::runtime_error("Malformed checkpoint matrix row");
        for (int k = 0; k < columns; ++k)
            m(i, k) = j[i][k].get<Scalar>();
    }
    return m;
}
json checkpoint(const wcso::TetMesh &mesh, const wcso::Surface &target, const wcso::RunState &state,
                const wcso::Normalize &normalization, const wcso::Config &config) {
    std::ostringstream random_state;
    random_state << state.rng;
    json samples = json::array();
    for (const auto &s : state.icp.samples)
        samples.push_back({s.face, s.barycentric[0], s.barycentric[1], s.barycentric[2]});
    return {{"schema", 1},
            {"config", config.to_json()},
            {"vertices", matrix_json(mesh.vertices)},
            {"rest_vertices", matrix_json(mesh.rest_vertices)},
            {"tets", matrix_json(mesh.tets)},
            {"target_vertices", matrix_json(target.vertices)},
            {"target_faces", matrix_json(target.faces)},
            {"offset", {normalization.offset[0], normalization.offset[1], normalization.offset[2]}},
            {"scale", normalization.scale},
            {"completed_outer", state.completed_outer},
            {"eta", state.eta},
            {"initial_boundary_stress", state.initial_boundary_stress},
            {"initial_all_stress", state.initial_all_stress},
            {"sigma_star", state.sigma_star},
            {"phase1_steps", state.phase1_steps},
            {"failed_steps", state.failed_steps},
            {"reached_phase2", state.reached_phase2},
            {"pressure", matrix_json(state.pressure)},
            {"rng", random_state.str()},
            {"samples", samples},
            {"samples_per_face", state.icp.samples_per_face}};
}
void restore(const json &j, wcso::TetMesh &mesh, wcso::Surface &target, wcso::RunState &state,
             wcso::Normalize &normalization) {
    if (j.at("schema").get<int>() != 1)
        throw std::runtime_error("Unsupported checkpoint schema");
    mesh.vertices = matrix_from<double>(j.at("vertices"), 3);
    mesh.tets = matrix_from<int>(j.at("tets"), 4);
    mesh.rebuild_topology();
    mesh.rest_vertices = matrix_from<double>(j.at("rest_vertices"), 3);
    if (mesh.rest_vertices.rows() != mesh.vertices.rows() || !mesh.rest_vertices.allFinite())
        throw std::runtime_error("Invalid checkpoint contact rest geometry");
    target.vertices = matrix_from<double>(j.at("target_vertices"), 3);
    target.faces = matrix_from<int>(j.at("target_faces"), 3);
    wcso::validate_surface(target);
    for (int i = 0; i < 3; ++i)
        normalization.offset[i] = j.at("offset").at(i).get<double>();
    normalization.scale = j.at("scale").get<double>();
    if (!(normalization.scale > 0) || !normalization.offset.allFinite())
        throw std::runtime_error("Invalid checkpoint normalization");
    state.completed_outer = j.at("completed_outer").get<int>();
    state.eta = j.at("eta").get<double>();
    state.sigma_star = j.at("sigma_star").get<double>();
    state.phase1_steps = j.at("phase1_steps").get<int>();
    state.failed_steps = j.at("failed_steps").get<int>();
    state.reached_phase2 = j.at("reached_phase2").get<bool>();
    state.pressure = matrix_from<double>(j.at("pressure"), 1);
    state.initial_boundary_stress = j.at("initial_boundary_stress").get<double>();
    state.initial_all_stress = j.at("initial_all_stress").get<double>();
    if (state.completed_outer < 0 || state.phase1_steps < 0 || state.failed_steps < 0 ||
        !(state.eta > 0 && state.eta < 1) || !(state.sigma_star > 0) ||
        !(state.initial_boundary_stress > 0) || !(state.initial_all_stress > 0) ||
        !state.pressure.allFinite() ||
        (state.pressure.size() != 0 &&
         state.pressure.size() != static_cast<int>(mesh.surface_to_volume.size())))
        throw std::runtime_error("Invalid checkpoint scheduling state");
    std::istringstream rng(j.at("rng").get<std::string>());
    if (!(rng >> state.rng))
        throw std::runtime_error("Invalid checkpoint RNG");
    auto temporary_rng = state.rng;
    state.icp = wcso::make_samples(mesh.surface(), 1, temporary_rng);
    state.icp.samples.clear();
    state.icp.samples_per_face = j.at("samples_per_face").get<std::vector<int>>();
    if (state.icp.samples_per_face.size() != static_cast<size_t>(mesh.faces.rows()))
        throw std::runtime_error("Checkpoint sampling topology mismatch");
    std::vector<int> actual(static_cast<size_t>(mesh.faces.rows()), 0);
    for (const auto &sample : j.at("samples")) {
        const int face = sample.at(0).get<int>();
        Eigen::Vector3d bary(sample.at(1).get<double>(), sample.at(2).get<double>(),
                             sample.at(3).get<double>());
        if (face < 0 || face >= mesh.faces.rows() || !bary.allFinite() || bary.minCoeff() < 0 ||
            std::abs(bary.sum() - 1) > 1e-12)
            throw std::runtime_error("Invalid checkpoint sample");
        ++actual[face];
        state.icp.samples.push_back({face, bary});
    }
    if (actual != state.icp.samples_per_face || *std::min_element(actual.begin(), actual.end()) <= 0)
        throw std::runtime_error("Checkpoint sample count mismatch");
}
} // namespace
int main(int argc, char **argv) {
    fs::path input, output, config_path, resume_path;
    json overrides = json::object();
    bool analyze_only = false;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            auto value = [&]() {
                if (++i >= argc)
                    throw std::invalid_argument("Missing value after " + arg);
                return std::string(argv[i]);
            };
            if (arg == "--help" || arg == "-h") {
                std::cout << "wcso --input MODEL.{obj,ply,stl,mesh,msh} --output DIR [--config FILE]\n"
                             "     [--threads N] [--seed N] [--solver cg|gcr] [--save-every N] [--debug]\n"
                             "     [--outer-iterations N] [--inner-iterations N] [--analyze-only]\n"
                             "wcso --resume CHECKPOINT.json --output NEW_DIR [--outer-iterations N]\n";
                return 0;
            } else if (arg == "--input")
                input = value();
            else if (arg == "--output")
                output = value();
            else if (arg == "--config")
                config_path = value();
            else if (arg == "--resume")
                resume_path = value();
            else if (arg == "--threads")
                overrides["threads"] = json::parse(value());
            else if (arg == "--seed")
                overrides["seed"] = json::parse(value());
            else if (arg == "--solver")
                overrides["linear_solver"] = value();
            else if (arg == "--save-every")
                overrides["save_every"] = json::parse(value());
            else if (arg == "--outer-iterations")
                overrides["outer_iterations"] = json::parse(value());
            else if (arg == "--inner-iterations")
                overrides["inner_iterations"] = json::parse(value());
            else if (arg == "--debug")
                overrides["debug"] = true;
            else if (arg == "--analyze-only")
                analyze_only = true;
            else
                throw std::invalid_argument("Unknown argument: " + arg);
        }
        if (output.empty() || (input.empty() && resume_path.empty()) ||
            (!input.empty() && !resume_path.empty()))
            throw std::invalid_argument("Provide --output and exactly one of --input or --resume");
        const bool resume = !resume_path.empty();
        json saved, settings = json::object();
        if (resume) {
            saved = read_json(resume_path);
            settings = saved.at("config");
        }
        if (!config_path.empty())
            settings.update(read_json(config_path));
        settings.update(overrides);
        const wcso::Config config = wcso::Config::from_json(settings);
        if (fs::exists(output) && (!fs::is_directory(output) || !fs::is_empty(output)))
            throw std::runtime_error("Output must be a new or empty directory");
        if (resume) {
            auto old = saved.at("config"), effective = config.to_json();
            for (const std::string key :
                 {"outer_iterations", "threads", "save_every", "checkpoint_every", "debug"}) {
                old.erase(key);
                effective.erase(key);
            }
            if (old != effective)
                throw std::invalid_argument("Resume must retain the saved numerical parameters");
        }
#ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(config.threads);
#endif
        wcso::TetMesh mesh;
        wcso::Surface target;
        wcso::Normalize normalization;
        wcso::RunState state;
        if (resume)
            restore(saved, mesh, target, state, normalization);
        else if (input.extension() == ".mesh" || input.extension() == ".msh") {
            mesh = wcso::read_tet_mesh(input);
            normalization = wcso::Normalize::from_vertices(mesh.vertices);
            mesh.set_vertices(normalization.apply(mesh.vertices));
            mesh.rebuild_topology();
            target = mesh.surface();
        } else {
            target = wcso::read_surface(input);
            normalization = wcso::Normalize::from_vertices(target.vertices);
            target.vertices = normalization.apply(target.vertices);
            wcso::validate_surface(target);
            if (wcso::has_intersections(target))
                throw std::runtime_error("Input surface self-intersects");
            mesh = wcso::tetrahedralize(target, config);
        }
        if (wcso::has_intersections(mesh.surface()))
            throw std::runtime_error("Volume mesh boundary self-intersects");
        fs::create_directories(output);
        write_json(output / "config_used.json", config.to_json());
        write_json(output / "normalization.json",
                   {{"offset", {normalization.offset[0], normalization.offset[1], normalization.offset[2]}},
                    {"scale", normalization.scale},
                    {"stress_units", "normalized analysis units"}});
        if (!resume && input.extension() != ".mesh" && input.extension() != ".msh") {
            auto generated = mesh;
            generated.vertices = normalization.inverse(generated.vertices);
            wcso::write_tet_mesh(output / "generated_input.mesh", generated);
        }
        if (analyze_only) {
            const wcso::FEM fem(mesh, config.lame_lambda, config.lame_mu);
            std::mt19937_64 random(config.seed);
            // Match a fresh optimization run's deterministic RNG consumption.
            (void)wcso::make_samples(mesh.surface(), config.sampling_density, random);
            const auto load = wcso::worst_case_load(fem, config, random);
            const auto stress = wcso::stress_barrier(mesh, fem, config, load.pressure, false);
            write_json(output / "analysis.json", {{"max_boundary_stress", stress.max_boundary},
                                                  {"max_all_stress", stress.max_all},
                                                  {"sigma_bar", stress.sigma_bar},
                                                  {"stress_energy", stress.value},
                                                  {"eigenvalues", matrix_json(load.eigenvalues)},
                                                  {"compliances", matrix_json(load.compliances)},
                                                  {"selected", load.selected},
                                                  {"eigen_backward_errors", matrix_json(load.residuals)},
                                                  {"lambda_max", load.lambda_max},
                                                  {"penalty", load.penalty},
                                                  {"boundary_tets", mesh.boundary_tets.size()},
                                                  {"all_tets", mesh.tets.rows()}});
            std::cout << "Analysis complete: boundary stress=" << stress.max_boundary
                      << ", all stress=" << stress.max_all << '\n';
            return 0;
        }
        std::ofstream history(output / "history.csv");
        history << std::setprecision(17)
                << "outer,inner,phase,accepted,sigma_bar,sigma_star,eta,max_boundary_stress,max_all_stress,F,"
                   "E_ss,E_si,E_ef,E_hd,E_lds,seconds,min_quality,event\n";
        wcso::BilevelSolver solver(std::move(mesh), target, config, std::move(state), resume);
        json final_statistics;
        solver.run([&](const wcso::Progress &p, const wcso::TetMesh &current,
                       const wcso::RunState &current_state) {
            history << p.outer << ',' << p.inner << ',' << p.phase << ',' << p.accepted << ',' << p.sigma_bar
                    << ',' << p.sigma_star << ',' << p.eta << ',' << p.max_boundary << ',' << p.max_all << ','
                    << p.objective << ',' << p.e_ss << ',' << p.e_si << ',' << p.e_ef << ',' << p.e_hd << ','
                    << p.e_lds << ',' << p.seconds << ',' << p.min_quality << ",\"" << p.event << "\"\n";
            history.flush();
            if (!history)
                throw std::runtime_error("Failed writing history.csv");
            if (config.debug || p.event == "outer_complete" || p.event == "final_analysis")
                std::cout << "outer=" << p.outer << " event=" << p.event
                          << " boundary_stress=" << p.max_boundary << " all_stress=" << p.max_all
                          << " eta=" << current_state.eta << " seconds=" << p.seconds << std::endl;
            if (p.event == "outer_complete") {
                if (config.checkpoint_every > 0 && p.outer % config.checkpoint_every == 0)
                    write_json(output / "checkpoint.json",
                               checkpoint(current, target, current_state, normalization, config));
                if (config.save_every > 0 && p.outer % config.save_every == 0) {
                    fs::create_directories(output / "snapshots");
                    auto s = current.surface();
                    s.vertices = normalization.inverse(s.vertices);
                    wcso::write_surface(
                        output / "snapshots" / ("iteration_" + std::to_string(p.outer) + ".obj"), s);
                }
            }
            if (p.event == "final_analysis")
                final_statistics = {{"completed_outer", p.outer},
                                    {"max_boundary_stress", p.max_boundary},
                                    {"max_all_stress", p.max_all},
                                    {"sigma_bar", p.sigma_bar},
                                    {"E_hd", p.e_hd},
                                    {"initial_boundary_stress", current_state.initial_boundary_stress},
                                    {"initial_all_stress", current_state.initial_all_stress},
                                    {"boundary_stress_reduction_percent",
                                     100 * (1 - p.max_boundary / current_state.initial_boundary_stress)},
                                    {"all_stress_reduction_percent",
                                     100 * (1 - p.max_all / current_state.initial_all_stress)},
                                    {"E_lds", p.e_lds},
                                    {"seconds", p.seconds},
                                    {"min_quality", p.min_quality},
                                    {"termination", current_state.eta > config.eta_stop
                                                        ? "eta_limit"
                                                        : "outer_iteration_limit"}};
        });
        auto result = solver.mesh().surface();
        result.vertices = normalization.inverse(result.vertices);
        wcso::write_surface(output / "final.obj", result);
        auto volume = solver.mesh();
        volume.vertices = normalization.inverse(volume.vertices);
        wcso::write_tet_mesh(output / "final.mesh", volume);
        write_json(output / "summary.json", final_statistics);
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "wcso: " << e.what() << '\n';
        return 1;
    }
}
