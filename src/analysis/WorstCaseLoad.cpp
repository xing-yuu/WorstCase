#include <Eigen/Eigenvalues>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <algorithm>
#include <stdexcept>
#include <wcso/analysis/WorstCaseLoad.h>

namespace wcso {
namespace {
class PenaltyOperator {
  public:
    using Scalar = double;
    const Sparse &base;
    const Eigen::MatrixXd &G;
    double penalty;
    PenaltyOperator(const Sparse &b, const Eigen::MatrixXd &g, double m) : base(b), G(g), penalty(m) {}
    Eigen::Index rows() const {
        return base.rows();
    }
    Eigen::Index cols() const {
        return base.cols();
    }
    void perform_op(const double *in, double *out) const {
        const Eigen::Map<const Eigen::VectorXd> x(in, cols());
        Eigen::Map<Eigen::VectorXd> y(out, rows());
        y = base * x + penalty * G.transpose() * (G * x);
    }
};
Eigen::VectorXd random_vector(int n, std::mt19937_64 &rng) {
    Eigen::VectorXd v(n);
    for (int i = 0; i < n; ++i)
        v[i] = 2 * std::generate_canonical<double, 53>(rng) - 1;
    return v;
}
} // namespace
WorstCaseResult worst_case_load(const FEM &f, const Config &c, std::mt19937_64 &rng,
                                const Eigen::VectorXd *previous) {
    WorstCaseResult result;
    const int n = static_cast<int>(f.K.rows()), ns = static_cast<int>(f.N.cols());
    if (ns < 4)
        throw std::runtime_error("At least four pressure degrees of freedom are required");
    Spectra::SparseSymMatProd<double> kop(f.K);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> max_solver(kop, 1, std::min(n, 24));
    auto initial = random_vector(n, rng);
    max_solver.init(initial.data());
    if (max_solver.compute(Spectra::SortRule::LargestAlge, c.eigen_max_iterations, c.eigen_tolerance) != 1)
        throw std::runtime_error("Largest stiffness eigenvalue did not converge");
    result.lambda_max = max_solver.eigenvalues()[0];
    result.penalty = result.lambda_max * c.nullspace_multiplier;
    const Sparse base = f.N.transpose() * f.K * f.N;
    const Eigen::MatrixXd G = f.rigid.transpose() * f.N;
    Sparse mass(ns, ns);
    for (int i = 0; i < ns; ++i)
        mass.insert(i, i) = f.mass[i];
    PenaltyOperator op(base, G, result.penalty);
    Eigen::MatrixXd vectors;
    if (ns <= 48) {
        const Eigen::MatrixXd dense = Eigen::MatrixXd(base) + result.penalty * G.transpose() * G;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eig(dense, Eigen::MatrixXd(mass));
        if (eig.info() != Eigen::Success)
            throw std::runtime_error("Small generalized eigenproblem failed");
        result.eigenvalues = eig.eigenvalues().head(3);
        vectors = eig.eigenvectors().leftCols(3);
    } else {
        Spectra::SparseCholesky<double> mass_op(mass);
        Spectra::SymGEigsSolver<PenaltyOperator, Spectra::SparseCholesky<double>,
                                Spectra::GEigsMode::Cholesky>
            eig(op, mass_op, 3, std::min(ns, 48));
        initial = previous && previous->size() == ns ? *previous : random_vector(ns, rng);
        eig.init(initial.data());
        const auto converged = eig.compute(Spectra::SortRule::SmallestAlge, c.eigen_max_iterations,
                                           c.eigen_tolerance, Spectra::SortRule::SmallestAlge);
        if (converged != 3 || eig.info() != Spectra::CompInfo::Successful)
            throw std::runtime_error("The three smallest generalized eigenpairs did not converge");
        result.eigenvalues = eig.eigenvalues();
        vectors = eig.eigenvectors();
    }
    ProjectedSolver solver(f, c);
    result.compliances.resize(3);
    result.residuals.resize(3);
    for (int j = 0; j < 3; ++j) {
        vectors.col(j) /= std::sqrt(vectors.col(j).dot(f.mass.cwiseProduct(vectors.col(j))));
        Eigen::Index max_index;
        vectors.col(j).cwiseAbs().maxCoeff(&max_index);
        if (vectors(max_index, j) < 0)
            vectors.col(j) *= -1;
        const Eigen::VectorXd load = f.project(f.N * vectors.col(j));
        const auto solved = solver.solve(load);
        result.compliances[j] = load.dot(solved.x);
        Eigen::VectorXd applied(ns);
        op.perform_op(vectors.col(j).data(), applied.data());
        const Eigen::VectorXd bp = f.mass.cwiseProduct(vectors.col(j));
        // Backward error uses the full operator scale, including the large penalty.
        const double operator_bound = base.norm() + result.penalty * G.squaredNorm();
        result.residuals[j] = (applied - result.eigenvalues[j] * bp).norm() /
                              std::max(1e-30, operator_bound * vectors.col(j).norm() +
                                                  std::abs(result.eigenvalues[j]) * bp.norm());
        if (!std::isfinite(result.residuals[j]) ||
            result.residuals[j] > std::max(1e-8, 100 * c.eigen_tolerance))
            throw std::runtime_error("Generalized eigenpair failed backward-error validation");
    }
    Eigen::Index selected;
    result.compliances.maxCoeff(&selected);
    result.selected = static_cast<int>(selected);
    result.pressure = vectors.col(selected);
    return result;
}
} // namespace wcso
