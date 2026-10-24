#include <chrono>
#include <cmath>
#include <stdexcept>
#include <wcso/optim/BilevelSolver.h>
#include <wcso/remesh/MeshImprove.h>

namespace wcso {
namespace {
struct Objective {
    Progress progress;
    Eigen::VectorXd gradient;
};
Objective evaluate(const TetMesh &mesh, const Config &config, const Eigen::VectorXd &pressure,
                   const ICPContext &context, int phase, double sigma_star, bool gradient) {
    Objective out;
    auto &p = out.progress;
    const auto flip = element_flip(mesh, config.d1, gradient);
    if (!std::isfinite(flip.value)) {
        p.objective = flip.value;
        return out;
    }
    const FEM fem(mesh, config.lame_lambda, config.lame_mu);
    const auto stress = stress_barrier(mesh, fem, config, pressure, gradient, phase == 2 ? sigma_star : -1);
    p.sigma_bar = stress.sigma_bar;
    p.max_boundary = stress.max_boundary;
    p.max_all = stress.max_all;
    if (!std::isfinite(stress.value)) {
        p.objective = stress.value;
        return out;
    }
    const auto contact = self_intersection(mesh, config.d0, gradient);
    p.e_ss = stress.value;
    p.e_si = contact.value;
    p.e_ef = flip.value;
    p.objective = config.w_ss * p.e_ss + config.w_si * p.e_si + config.w_ef * p.e_ef;
    if (gradient)
        out.gradient =
            config.w_ss * stress.gradient + config.w_si * contact.gradient + config.w_ef * flip.gradient;
    if (phase == 2) {
        const auto shape = shape_fidelity(mesh, context, gradient);
        const auto smooth = local_smoothing(mesh, context, config.d2, gradient);
        p.e_hd = shape.value;
        p.e_lds = smooth.value;
        p.objective += config.w_hd * p.e_hd + config.w_lds * p.e_lds;
        if (gradient)
            out.gradient += config.w_hd * shape.gradient + config.w_lds * smooth.gradient;
    }
    return out;
}
} // namespace
BilevelSolver::BilevelSolver(TetMesh mesh, Surface target, Config config, RunState state, bool resume)
    : mesh_(std::move(mesh)), target_(std::move(target)), config_(std::move(config)),
      state_(std::move(state)) {
    config_.validate();
    if (!resume) {
        state_.eta = config_.eta_initial;
        state_.rng.seed(config_.seed);
        state_.icp = make_samples(mesh_.surface(), config_.sampling_density, state_.rng);
    }
    update_correspondences(state_.icp, mesh_.surface(), target_);
}
void BilevelSolver::run(
    const std::function<void(const Progress &, const TetMesh &, const RunState &)> &callback) {
    const auto started = std::chrono::steady_clock::now();
    auto elapsed = [&]() {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    };
    for (int outer = state_.completed_outer;
         outer < config_.outer_iterations && state_.eta <= config_.eta_stop; ++outer) {
        const FEM initial(mesh_, config_.lame_lambda, config_.lame_mu);
        const auto worst = worst_case_load(initial, config_, state_.rng, &state_.pressure);
        state_.pressure = worst.pressure;
        const auto initial_stress = stress_barrier(mesh_, initial, config_, state_.pressure, false);
        const double reference_bar = initial_stress.sigma_bar;
        state_.sigma_star = state_.eta * reference_bar;
        state_.reached_phase2 = false;
        if (outer == 0) {
            state_.initial_boundary_stress = initial_stress.max_boundary;
            state_.initial_all_stress = initial_stress.max_all;
            Progress first;
            first.event = "initial_analysis";
            first.sigma_bar = initial_stress.sigma_bar;
            first.sigma_star = state_.sigma_star;
            first.eta = state_.eta;
            first.max_boundary = initial_stress.max_boundary;
            first.max_all = initial_stress.max_all;
            first.e_ss = initial_stress.value;
            first.e_si = self_intersection(mesh_, config_.d0, false).value;
            first.e_ef = element_flip(mesh_, config_.d1, false).value;
            first.objective =
                config_.w_ss * first.e_ss + config_.w_si * first.e_si + config_.w_ef * first.e_ef;
            first.seconds = elapsed();
            first.min_quality = worst_quality(mesh_);
            callback(first, mesh_, state_);
        }
        Progress latest;
        for (int inner = 0; inner < config_.inner_iterations && state_.eta <= config_.eta_stop; ++inner) {
            const FEM current(mesh_, config_.lame_lambda, config_.lame_mu);
            const auto current_stress = stress_barrier(mesh_, current, config_, state_.pressure, false);
            const int phase = current_stress.sigma_bar <= state_.sigma_star ? 2 : 1;
            if (phase == 2) {
                state_.reached_phase2 = true;
                state_.phase1_steps = state_.failed_steps = 0;
            }
            update_correspondences(state_.icp, mesh_.surface(), target_);
            auto value =
                evaluate(mesh_, config_, state_.pressure, state_.icp, phase, state_.sigma_star, true);
            if (!std::isfinite(value.progress.objective) || !value.gradient.allFinite())
                throw std::runtime_error("Non-finite objective or gradient at the current iterate");
            double max_gradient = 0;
            for (int id : mesh_.surface_to_volume)
                max_gradient = std::max(max_gradient, value.gradient.segment<3>(3 * id).norm());
            bool accepted = false;
            const double norm_squared = value.gradient.squaredNorm();
            double step = max_gradient > 0
                              ? std::min(1.0, config_.step_fraction * mesh_.min_edge_length() / max_gradient)
                              : 0;
            const auto previous_surface = mesh_.surface();
            for (int trial = 0;
                 trial < config_.max_backtracks && step * max_gradient > config_.step_tolerance; ++trial) {
                TetMesh candidate = mesh_;
                Eigen::MatrixXd coordinates = mesh_.vertices;
                for (int id : mesh_.surface_to_volume)
                    coordinates.row(id) -= step * value.gradient.segment<3>(3 * id).transpose();
                candidate.set_vertices(coordinates);
                if (candidate.positive_volumes() &&
                    collision_free(previous_surface, candidate.surface().vertices)) {
                    try {
                        auto proposed = evaluate(candidate, config_, state_.pressure, state_.icp, phase,
                                                 state_.sigma_star, false);
                        if (std::isfinite(proposed.progress.objective) &&
                            proposed.progress.objective <=
                                value.progress.objective - config_.armijo * step * norm_squared) {
                            mesh_ = std::move(candidate);
                            value.progress = proposed.progress;
                            accepted = true;
                            break;
                        }
                    } catch (const std::runtime_error &) {
                        // A trial with failed analysis is inadmissible; keep the current state.
                    }
                }
                step *= .5;
            }
            latest = value.progress;
            latest.outer = outer + 1;
            latest.inner = inner + 1;
            latest.phase = phase;
            latest.accepted = accepted;
            latest.sigma_star = state_.sigma_star;
            latest.eta = state_.eta;
            latest.seconds = elapsed();
            latest.min_quality = worst_quality(mesh_);
            latest.event = accepted ? "step_accepted" : "step_failed";
            if (phase == 1) {
                ++state_.phase1_steps;
                state_.failed_steps = accepted ? 0 : state_.failed_steps + 1;
                if (!state_.reached_phase2 && (state_.failed_steps >= config_.max_failed_steps ||
                                               state_.phase1_steps >= config_.max_phase1_steps)) {
                    state_.eta = (state_.eta + 1) / 2;
                    state_.sigma_star = state_.eta * reference_bar;
                    state_.phase1_steps = state_.failed_steps = 0;
                    latest.event += ";eta_relaxed";
                }
            }
            callback(latest, mesh_, state_);
        }
        if (config_.remesh && worst_quality(mesh_) < config_.quality_threshold) {
            std::string diagnostic;
            if (improve_mesh(mesh_, config_, diagnostic)) {
                state_.pressure.resize(0);
                state_.icp = make_samples(mesh_.surface(), config_.sampling_density, state_.rng);
                update_correspondences(state_.icp, mesh_.surface(), target_);
                const FEM improved(mesh_, config_.lame_lambda, config_.lame_mu);
                state_.pressure = worst_case_load(improved, config_, state_.rng).pressure;
                latest = evaluate(mesh_, config_, state_.pressure, state_.icp, 1, 0, false).progress;
                latest.outer = outer + 1;
            }
            Progress event;
            event.outer = outer + 1;
            event.seconds = elapsed();
            event.min_quality = worst_quality(mesh_);
            event.event = diagnostic;
            callback(event, mesh_, state_);
        }
        state_.completed_outer = outer + 1;
        latest.event = "outer_complete";
        latest.eta = state_.eta;
        latest.sigma_star = state_.sigma_star;
        latest.seconds = elapsed();
        latest.min_quality = worst_quality(mesh_);
        callback(latest, mesh_, state_);
    }
    // Reanalyse the final shape; this diagnostic must not alter resume RNG state.
    const FEM final_fem(mesh_, config_.lame_lambda, config_.lame_mu);
    auto diagnostic_rng = state_.rng;
    const auto worst = worst_case_load(final_fem, config_, diagnostic_rng, &state_.pressure);
    auto context = state_.icp;
    update_correspondences(context, mesh_.surface(), target_);
    auto final = evaluate(mesh_, config_, worst.pressure, context, 1, 0, false).progress;
    final.e_hd = shape_fidelity(mesh_, context, false).value;
    final.e_lds = local_smoothing(mesh_, context, config_.d2, false).value;
    final.outer = state_.completed_outer;
    final.event = "final_analysis";
    final.seconds = elapsed();
    final.min_quality = worst_quality(mesh_);
    final.eta = state_.eta;
    final.sigma_star = state_.sigma_star;
    callback(final, mesh_, state_);
}
} // namespace wcso
