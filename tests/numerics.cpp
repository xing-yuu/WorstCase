#include "../src/remesh/BoundaryPreservation.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <set>
#include <stdexcept>
#include <wcso/energy/Energies.h>
#include <wcso/mesh/MeshIO.h>
#include <wcso/remesh/MeshImprove.h>

namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}
wcso::TetMesh grid(int subdivisions) {
    const int side = subdivisions + 1;
    auto index = [side](int x, int y, int z) { return x + side * (y + side * z); };
    wcso::TetMesh mesh;
    mesh.vertices.resize(side * side * side, 3);
    for (int z = 0; z < side; ++z)
        for (int y = 0; y < side; ++y)
            for (int x = 0; x < side; ++x)
                mesh.vertices.row(index(x, y, z)) << double(x) / subdivisions, double(y) / subdivisions,
                    double(z) / subdivisions;
    mesh.tets.resize(6 * subdivisions * subdivisions * subdivisions, 4);
    const int permutations[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    int t = 0;
    for (int z = 0; z < subdivisions; ++z)
        for (int y = 0; y < subdivisions; ++y)
            for (int x = 0; x < subdivisions; ++x)
                for (const auto &p : permutations) {
                    int coordinates[3] = {x, y, z};
                    mesh.tets(t, 0) = index(x, y, z);
                    for (int k = 0; k < 3; ++k) {
                        ++coordinates[p[k]];
                        mesh.tets(t, k + 1) = index(coordinates[0], coordinates[1], coordinates[2]);
                    }
                    ++t;
                }
    mesh.rebuild_topology();
    return mesh;
}
Eigen::VectorXd random_vector(int size, std::mt19937_64 &rng) {
    Eigen::VectorXd v(size);
    for (int i = 0; i < size; ++i)
        v[i] = 2 * std::generate_canonical<double, 53>(rng) - 1;
    return v;
}
void check_gradient(const std::string &name, const wcso::TetMesh &mesh,
                    const std::function<wcso::EnergyValue(const wcso::TetMesh &, bool)> &evaluate,
                    double tolerance = 2e-5) {
    const auto analytical = evaluate(mesh, true);
    require(std::isfinite(analytical.value) && analytical.gradient.allFinite(),
            name + ": nonfinite evaluation");
    std::mt19937_64 rng(538);
    double worst = 0;
    for (int trial = 0; trial < 6; ++trial) {
        auto direction = random_vector(static_cast<int>(mesh.vertices.size()), rng);
        for (int i = 0; i < mesh.vertices.rows(); ++i)
            if (mesh.volume_to_surface[i] < 0)
                direction.segment<3>(3 * i).setZero();
        direction.normalize();
        const double exact = analytical.gradient.dot(direction);
        double best = std::numeric_limits<double>::infinity();
        for (double step : {2e-5, 5e-6, 1e-6}) {
            auto plus = mesh, minus = mesh;
            for (int i = 0; i < mesh.vertices.rows(); ++i) {
                plus.vertices.row(i) += step * direction.segment<3>(3 * i).transpose();
                minus.vertices.row(i) -= step * direction.segment<3>(3 * i).transpose();
            }
            const double finite = (evaluate(plus, false).value - evaluate(minus, false).value) / (2 * step);
            best = std::min(best,
                            std::abs(finite - exact) / std::max({1e-7, std::abs(exact), std::abs(finite)}));
        }
        worst = std::max(worst, best);
    }
    std::cout << name << " directional_gradient_relative_error=" << worst << '\n';
    require(worst < tolerance, name + ": finite-difference gradient mismatch");
}
void test_fem() {
    auto mesh = grid(2);
    const double lame_lambda = 1.0, lame_mu = 0.3;
    wcso::FEM fem(mesh, lame_lambda, lame_mu);
    const Eigen::MatrixXd K(fem.K);
    require((K - K.transpose()).norm() < 1e-12, "K symmetry");
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(K);
    require((eigen.eigenvalues().array().abs() < 1e-9).count() == 6, "Exactly six rigid modes");
    require((K * fem.rigid).norm() < 1e-12, "Rigid-mode annihilation");
    Eigen::VectorXd displacement(3 * mesh.vertices.rows());
    const double strain = 1e-4;
    const double lateral_strain = -lame_lambda / (2 * (lame_lambda + lame_mu)) * strain;
    const double axial_stress = (lame_lambda * (1 + 2 * lateral_strain / strain) + 2 * lame_mu) * strain;
    for (int i = 0; i < mesh.vertices.rows(); ++i)
        displacement.segment<3>(3 * i) =
            Eigen::Vector3d(strain * mesh.vertices(i, 0), lateral_strain * mesh.vertices(i, 1),
                            lateral_strain * mesh.vertices(i, 2));
    const auto stresses = fem.stresses(mesh, displacement);
    for (int i = 0; i < stresses.rows(); ++i) {
        require(std::abs(stresses(i, 0) - axial_stress) < 1e-14, "Uniaxial bar stress");
        require(stresses.row(i).tail<5>().norm() < 1e-14, "Constant-strain patch");
    }
    require(std::abs(mesh.surface_areas().sum() - 6) < 1e-12, "Lumped surface mass conservation");
    const auto kernel = wcso::element_kernel(wcso::tet_vertices(mesh, 0), true);
    for (int k = 0; k < 12; ++k) {
        auto x = wcso::tet_vertices(mesh, 0), y = x;
        x(k / 3, k % 3) += 1e-6;
        y(k / 3, k % 3) -= 1e-6;
        const wcso::Matrix612 finite =
            (wcso::element_kernel(x, false).B - wcso::element_kernel(y, false).B) / 2e-6;
        require((finite - kernel.dB[k]).norm() < 1e-8, "Analytic B derivative");
    }
    wcso::Config c;
    c.linear_tolerance = 1e-11;
    std::mt19937_64 rng(1);
    const auto exact = fem.project(random_vector(static_cast<int>(K.rows()), rng));
    const Eigen::VectorXd rhs = fem.K * exact;
    Eigen::VectorXd cg;
    for (const std::string method : {"cg", "gcr"}) {
        c.linear_solver = method;
        wcso::ProjectedSolver solver(fem, c);
        const auto solved = solver.solve(rhs);
        require((solved.x - exact).norm() / exact.norm() < 1e-8, method + " manufactured solution");
        if (method == "cg")
            cg = solved.x;
        else
            require((cg - solved.x).norm() < 1e-8, "CG/GCR agreement");
        std::cout << method << " residual=" << solved.residual << " iterations=" << solved.iterations << '\n';
    }
}
void test_energies() {
    auto mesh = grid(3);
    require(mesh.boundary_tets.size() < static_cast<size_t>(mesh.tets.rows()),
            "Test must contain nonboundary tetrahedra");
    wcso::Config c;
    c.linear_tolerance = 1e-11;
    std::mt19937_64 rng(33);
    auto p = random_vector(static_cast<int>(mesh.surface_to_volume.size()), rng);
    p /= std::sqrt(p.dot(mesh.surface_areas().cwiseProduct(p)));
    auto stress = [&](const wcso::TetMesh &m, bool g) {
        wcso::FEM fem(m, 1, .3);
        return wcso::stress_barrier(m, fem, c, p, g);
    };
    const auto before = stress(mesh, false);
    double sum = before.sigma_bar;
    for (int t : mesh.boundary_tets) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(
            wcso::stress_tensor(before.stresses.row(t).transpose()));
        for (double v : eig.eigenvalues())
            sum -= c.gamma * (std::log(before.sigma_bar - v) + std::log(before.sigma_bar + v));
    }
    require(std::abs(sum - before.value) < 1e-9, "Boundary-only energy sum");
    check_gradient("E_ss_minimized", mesh, stress, 2e-4);
    const double bound = before.max_boundary + 1.0;
    check_gradient(
        "E_ss_fixed", mesh,
        [&](const wcso::TetMesh &m, bool g) {
            wcso::FEM fem(m, 1, .3);
            return wcso::stress_barrier(m, fem, c, p, g, bound);
        },
        2e-4);
    check_gradient("E_ef", mesh,
                   [](const wcso::TetMesh &m, bool g) { return wcso::element_flip(m, .02, g); });
    auto context = wcso::make_samples(mesh.surface(), 40, rng);
    auto target = mesh.surface();
    target.vertices.col(0).array() -= .025;
    wcso::update_correspondences(context, mesh.surface(), target);
    check_gradient("E_hd", mesh,
                   [&](const wcso::TetMesh &m, bool g) { return wcso::shape_fidelity(m, context, g); });
    check_gradient("E_lds", mesh,
                   [&](const wcso::TetMesh &m, bool g) { return wcso::local_smoothing(m, context, .1, g); });
    auto rotation = Eigen::AngleAxisd(.44, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
    auto rotated = mesh;
    rotated.vertices = mesh.vertices * rotation.transpose();
    require(std::abs(stress(rotated, false).max_boundary - before.max_boundary) < 1e-8,
            "Principal stress rotation invariance");
    wcso::TetMesh contact;
    contact.vertices.resize(8, 3);
    contact.vertices << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, .5, .1, .1, -.01, .9, .1, -.01, .1, .9, -.01, .1, .1,
        -.5;
    contact.tets.resize(2, 4);
    contact.tets << 0, 1, 2, 3, 4, 5, 6, 7;
    contact.rebuild_topology();
    require(!wcso::has_intersections(contact.surface()), "Separated contact fixture");
    require(wcso::self_intersection(contact, .1, false).value > 0, "Active collision barrier fixture");
    check_gradient(
        "E_si", contact, [](const wcso::TetMesh &m, bool g) { return wcso::self_intersection(m, .1, g); },
        2e-4);
    auto end = contact.surface().vertices;
    for (int i = 4; i < 8; ++i)
        end(i, 2) += .6;
    require(!wcso::collision_free(contact.surface(), end), "CCD crossing detection");
}
void test_inner() {
    auto mesh = grid(3);
    mesh.vertices.col(1) *= .35;
    mesh.vertices.col(2) *= .2;
    wcso::FEM fem(mesh, 1, .3);
    wcso::Config c;
    c.linear_tolerance = 1e-10;
    std::mt19937_64 rng(5);
    const auto load = wcso::worst_case_load(fem, c, rng);
    require(load.eigenvalues.size() == 3, "Three inner modes");
    require(std::abs(load.pressure.dot(fem.mass.cwiseProduct(load.pressure)) - 1) < 1e-10,
            "Pressure normalization");
    require(load.compliances[load.selected] >= load.compliances.maxCoeff() - 1e-10,
            "Maximum-compliance selection");
    const Eigen::MatrixXd G = fem.rigid.transpose() * fem.N;
    const Eigen::MatrixXd A =
        Eigen::MatrixXd(fem.N.transpose() * fem.K * fem.N) + load.penalty * G.transpose() * G;
    const Eigen::MatrixXd B = fem.mass.asDiagonal();
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> dense(A, B);
    require((load.eigenvalues - dense.eigenvalues().head(3)).norm() / load.eigenvalues.norm() < 2e-3,
            "Dense generalized eigenvalue reference");
    std::cout << "inner eigenvalues=" << load.eigenvalues.transpose()
              << " errors=" << load.residuals.transpose() << '\n';
}
wcso::TetMesh deformed_remesh_fixture() {
    auto mesh = grid(2);
    mesh.vertices(13, 0) = .92;
    const auto rotation = Eigen::AngleAxisd(.37, Eigen::Vector3d(1, 2, 3).normalized()).toRotationMatrix();
    mesh.vertices = (mesh.vertices * rotation.transpose()).eval();
    mesh.vertices.rowwise() += Eigen::Vector3d(.123456789, -.0723456789, .01123456789).transpose();
    return mesh;
}
void test_boundary_restoration() {
    const auto original = deformed_remesh_fixture();
    auto candidate = original;
    candidate.vertices = original.vertices.colwise().reverse().eval();
    candidate.tets = (original.vertices.rows() - 1 - original.tets.array()).matrix();
    candidate.rebuild_topology();
    const auto expected = candidate.vertices;
    for (int vertex : candidate.surface_to_volume)
        for (int d = 0; d < 3; ++d)
            candidate.vertices(vertex, d) =
                std::nextafter(candidate.vertices(vertex, d), std::numeric_limits<double>::infinity());
    require((candidate.vertices - expected).norm() > 0, "Nonbinary boundary roundoff fixture");
    wcso::detail::restore_required_boundary(original, candidate);
    require((candidate.vertices.array() == expected.array()).all(),
            "Restore exact boundary after roundoff and vertex renumbering");

    auto rejects = [](const wcso::TetMesh &source, wcso::TetMesh &test, const std::string &reason) {
        const auto saved = test.vertices;
        bool rejected = false;
        try {
            wcso::detail::restore_required_boundary(source, test);
        } catch (const std::runtime_error &error) {
            rejected = std::string(error.what()).find(reason) != std::string::npos;
        }
        require(rejected, "Reject invalid MMG boundary: " + reason);
        require((test.vertices.array() == saved.array()).all(), "Boundary rejection must be atomic");
    };
    candidate = original;
    candidate.vertices(candidate.surface_to_volume[0], 0) += 1e-8;
    rejects(original, candidate, "beyond floating-point roundoff");
    candidate = original;
    candidate.vertices.row(candidate.surface_to_volume[1]) =
        candidate.vertices.row(candidate.surface_to_volume[0]);
    rejects(original, candidate, "not bijective");
    candidate = original;
    candidate.faces.row(0) = candidate.faces.row(1);
    rejects(original, candidate, "triangle topology changed");
    auto ambiguous = original;
    ambiguous.vertices.row(ambiguous.surface_to_volume[1]) =
        ambiguous.vertices.row(ambiguous.surface_to_volume[0]);
    const int close_vertex = ambiguous.surface_to_volume[1];
    ambiguous.vertices(close_vertex, 0) =
        std::nextafter(ambiguous.vertices(close_vertex, 0), std::numeric_limits<double>::infinity());
    candidate = ambiguous;
    rejects(ambiguous, candidate, "ambiguous");
}
void test_mesh() {
    test_boundary_restoration();
    auto mesh = deformed_remesh_fixture();
    const auto normalize = wcso::Normalize::from_vertices(mesh.vertices);
    require((normalize.inverse(normalize.apply(mesh.vertices)) - mesh.vertices).norm() < 1e-14,
            "Normalization round trip");
    require(mesh.positive_volumes(), "Distorted MMG fixture");
    const double quality = wcso::worst_quality(mesh);
    const auto boundary = mesh.surface();
    wcso::Config c;
    c.tet_edge_length = .5;
    std::string diagnostic;
    const bool accepted = wcso::improve_mesh(mesh, c, diagnostic);
    std::cout << "MMG " << diagnostic << " quality=" << quality << " -> " << wcso::worst_quality(mesh)
              << '\n';
    require(accepted, "MMG should improve the distorted interior fixture");
    require(wcso::worst_quality(mesh) > quality, "MMG quality improvement");
    auto boundary_signature = [](const wcso::Surface &surface) {
        using Point = std::array<double, 3>;
        std::multiset<std::array<Point, 3>> faces;
        for (int f = 0; f < surface.faces.rows(); ++f) {
            std::array<Point, 3> face;
            for (int k = 0; k < 3; ++k)
                for (int d = 0; d < 3; ++d)
                    face[k][d] = surface.vertices(surface.faces(f, k), d);
            std::sort(face.begin(), face.end());
            faces.insert(face);
        }
        return faces;
    };
    require(boundary_signature(mesh.surface()) == boundary_signature(boundary),
            "MMG must preserve exact boundary coordinates and triangles");
    require((mesh.rest_vertices.array() == mesh.vertices.array()).all(),
            "Accepted MMG reference coordinates must include exact boundary restoration");
    const auto saved = mesh;
    c.tet_edge_length = -1;
    require(!wcso::improve_mesh(mesh, c, diagnostic), "Failed MMG call should roll back");
    require((mesh.vertices - saved.vertices).norm() == 0 && mesh.tets == saved.tets, "Atomic MMG rollback");
    auto invalid = boundary;
    invalid.faces.row(0).swap(invalid.faces.row(1));
    invalid.faces.row(0) = invalid.faces.row(1);
    bool rejected = false;
    try {
        wcso::validate_surface(invalid);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    require(rejected, "Reject duplicate surface faces");
}
} // namespace
int main() {
    try {
        test_fem();
        test_energies();
        test_inner();
        test_mesh();
        std::cout << "All numerical tests passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
