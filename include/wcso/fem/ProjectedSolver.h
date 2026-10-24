#pragma once
#include <memory>
#include <wcso/core/Config.h>
#include <wcso/fem/Elasticity.h>

namespace wcso {
/// Iterative solve with independently recomputed projected relative residual.
struct SolveResult {
    Eigen::VectorXd x;
    int iterations = 0;
    double residual = 0;
};
/// AMG-preconditioned CG or restarted GCR, both on the rigid-mode complement.
class ProjectedSolver {
  public:
    /// Build an AMG preconditioner for this fixed FEM geometry.
    ProjectedSolver(const FEM &fem, const Config &config);
    ~ProjectedSolver();
    /// Solve K*x=P*b with x in range(P); throw on nonconvergence.
    SolveResult solve(const Eigen::VectorXd &rhs, const Eigen::VectorXd *initial = nullptr) const;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace wcso
