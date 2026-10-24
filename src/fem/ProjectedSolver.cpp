#include <algorithm>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <wcso/fem/ProjectedSolver.h>

namespace wcso {
struct ProjectedSolver::Impl {
    using AMG = amgcl::amg<amgcl::backend::builtin<double>, amgcl::coarsening::smoothed_aggregation,
                           amgcl::relaxation::damped_jacobi>;
    const FEM &fem;
    Config config;
    std::unique_ptr<AMG> preconditioner;
    Impl(const FEM &f, const Config &c) : fem(f), config(c) {
        Eigen::SparseMatrix<double, Eigen::RowMajor> matrix = f.K;
        const double scale = f.K.diagonal().cwiseAbs().maxCoeff();
        // Regularization belongs only to the preconditioner, never the equation.
        for (int i = 0; i < matrix.rows(); ++i)
            matrix.coeffRef(i, i) += 1e-8 * scale;
        matrix.makeCompressed();
        AMG::params params;
        params.coarsening.aggr.block_size = 3;
        params.coarsening.aggr.eps_strong = 0;
        params.npre = params.npost = 2;
        params.relax.damping = .5;
        params.coarse_enough = 300;
        preconditioner = std::make_unique<AMG>(matrix, params);
    }
    Eigen::VectorXd precondition(const Eigen::VectorXd &r) const {
        std::vector<double> rhs(r.data(), r.data() + r.size()), x(static_cast<size_t>(r.size()), 0);
        preconditioner->apply(rhs, x);
        return fem.project(Eigen::Map<Eigen::VectorXd>(x.data(), r.size()));
    }
};
ProjectedSolver::ProjectedSolver(const FEM &fem, const Config &config)
    : impl_(std::make_unique<Impl>(fem, config)) {}
ProjectedSolver::~ProjectedSolver() = default;
SolveResult ProjectedSolver::solve(const Eigen::VectorXd &rhs, const Eigen::VectorXd *initial) const {
    const auto &f = impl_->fem;
    const auto &c = impl_->config;
    const Eigen::VectorXd b = f.project(rhs);
    const double bn = b.norm();
    SolveResult out;
    out.x = Eigen::VectorXd::Zero(b.size());
    if (initial && initial->size() == b.size())
        out.x = f.project(*initial);
    if (bn < 1e-30) {
        out.x.setZero();
        return out;
    }
    Eigen::VectorXd r = b - f.project(f.K * out.x);
    const double target = c.linear_tolerance * bn;
    if (c.linear_solver == "cg") {
        Eigen::VectorXd z = impl_->precondition(r), p = z;
        double rho = r.dot(z);
        while (r.norm() > target && out.iterations < c.linear_max_iterations) {
            const Eigen::VectorXd ap = f.project(f.K * p);
            const double denom = p.dot(ap);
            if (!(rho > 0 && denom > 0) || !std::isfinite(denom))
                throw std::runtime_error("Projected CG lost positive definiteness");
            const double alpha = rho / denom;
            out.x += alpha * p;
            r -= alpha * ap;
            ++out.iterations;
            if (r.norm() <= target)
                break;
            z = impl_->precondition(r);
            const double next = r.dot(z);
            p = z + (next / rho) * p;
            rho = next;
        }
    } else {
        std::vector<Eigen::VectorXd> directions, images;
        while (r.norm() > target && out.iterations < c.linear_max_iterations) {
            if (static_cast<int>(directions.size()) == c.gcr_restart) {
                directions.clear();
                images.clear();
                r = b - f.project(f.K * out.x);
            }
            Eigen::VectorXd p = impl_->precondition(r), ap = f.project(f.K * p);
            // Reorthogonalize for the large condition numbers caused by thin parts.
            for (int pass = 0; pass < 2; ++pass)
                for (size_t j = 0; j < images.size(); ++j) {
                    const double beta = images[j].dot(ap);
                    ap -= beta * images[j];
                    p -= beta * directions[j];
                }
            const double norm = ap.norm();
            if (!(norm > 1e-30) || !std::isfinite(norm))
                break;
            ap /= norm;
            p /= norm;
            const double alpha = ap.dot(r);
            out.x += alpha * p;
            r -= alpha * ap;
            ++out.iterations;
            images.push_back(std::move(ap));
            directions.push_back(std::move(p));
        }
    }
    out.x = f.project(out.x);
    out.residual = (f.project(f.K * out.x) - b).norm() / bn;
    if (!out.x.allFinite() || out.residual > 5 * c.linear_tolerance) {
        std::ostringstream diagnostic;
        diagnostic << std::scientific << std::setprecision(17)
                   << "Linear solve did not converge: method=" << c.linear_solver
                   << " iterations=" << out.iterations << "/" << c.linear_max_iterations
                   << " tolerance=" << c.linear_tolerance << " recursive_relative_residual=" << r.norm() / bn
                   << " true_projected_relative_residual=" << out.residual << " rhs_norm=" << bn;
        throw std::runtime_error(diagnostic.str());
    }
    return out;
}
} // namespace wcso
