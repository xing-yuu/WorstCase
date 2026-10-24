#pragma once
#include <wcso/energy/Energies.h>

namespace wcso {
/// Serializable scheduling state; counters include explicitly failed attempts.
struct RunState {
    int completed_outer = 0, phase1_steps = 0, failed_steps = 0;
    double eta = 0.9, sigma_star = 0;
    double initial_boundary_stress = 0, initial_all_stress = 0;
    bool reached_phase2 = false;
    Eigen::VectorXd pressure;
    std::mt19937_64 rng;
    ICPContext icp;
};
/// One accepted step or explicit failed attempt, for logging by the application.
struct Progress {
    int outer = 0, inner = 0, phase = 1;
    bool accepted = false;
    double sigma_bar = 0, sigma_star = 0, eta = 0, max_boundary = 0, max_all = 0;
    double objective = 0, e_ss = 0, e_si = 0, e_ef = 0, e_hd = 0, e_lds = 0;
    double seconds = 0, min_quality = 0;
    std::string event;
};
/// Phase I/II optimizer: owns its mesh and state, and performs no file I/O.
class BilevelSolver {
  public:
    /// Construct a fresh run or resume an explicitly supplied state.
    BilevelSolver(TetMesh mesh, Surface target, Config config, RunState state, bool resume = false);
    /// Run Algorithm 1; callbacks may write logs, snapshots, and checkpoints.
    void run(const std::function<void(const Progress &, const TetMesh &, const RunState &)> &callback);
    /// Final normalized volume mesh.
    const TetMesh &mesh() const {
        return mesh_;
    }
    /// Final scheduling and sampling state.
    const RunState &state() const {
        return state_;
    }

  private:
    TetMesh mesh_;
    Surface target_;
    Config config_;
    RunState state_;
};
} // namespace wcso
