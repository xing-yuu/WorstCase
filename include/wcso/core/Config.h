#pragma once
#include <cstdint>
#include <nlohmann/json.hpp>
#include <string>

namespace wcso {
/// Numerical parameters in longest-bounding-box-side normalized coordinates.
struct Config {
    double lame_lambda = 1.0, lame_mu = 0.3;
    double gamma = 0.5, nullspace_multiplier = 1e7;
    double eigen_tolerance = 1e-10, linear_tolerance = 1e-8;
    int eigen_max_iterations = 10000, linear_max_iterations = 3000, gcr_restart = 40;
    std::string linear_solver = "gcr";
    double w_ss = 1, w_si = 1, w_ef = 1, w_hd = 1, w_lds = 1;
    double d0 = 1e-3, d1 = 1e-3, d2 = 1e-3;
    double epsilon_d = 1e-4, quality_threshold = 0.1, tet_edge_length = 0.025;
    double eta_initial = 0.9, eta_stop = 0.99;
    int outer_iterations = 1000, inner_iterations = 5;
    int max_phase1_steps = 15, max_failed_steps = 2, max_backtracks = 50;
    double armijo = 1e-4, step_fraction = 0.25, step_tolerance = 1e-12;
    double sampling_density = 10000;
    int threads = 1, save_every = 0, checkpoint_every = 1;
    std::uint64_t seed = 1;
    bool remesh = true, debug = false;
    /// Reject unknown keys, wrong types, and out-of-range parameters.
    static Config from_json(const nlohmann::json &json);
    /// Serialize every effective numerical parameter for reproducibility.
    nlohmann::json to_json() const;
    /// Check admissible material, solver, and geometric parameter ranges.
    void validate() const;
};
} // namespace wcso
