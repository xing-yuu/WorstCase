#pragma once
#include <random>
#include <wcso/fem/ProjectedSolver.h>

namespace wcso {
/// Three converged candidate modes and the selected maximum-compliance load.
struct WorstCaseResult {
    Eigen::VectorXd pressure, eigenvalues, compliances, residuals;
    int selected = 0;
    double lambda_max = 0, penalty = 0;
};
/// Solve the three-smallest-eigenpair surrogate, then compare true compliance.
WorstCaseResult worst_case_load(const FEM &fem, const Config &config, std::mt19937_64 &rng,
                                const Eigen::VectorXd *previous = nullptr);
} // namespace wcso
