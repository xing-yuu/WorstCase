#pragma once
#include <Eigen/Cholesky>
#include <Eigen/SparseCore>
#include <array>
#include <wcso/mesh/TetMesh.h>

namespace wcso {
using Sparse = Eigen::SparseMatrix<double>;
using Matrix6 = Eigen::Matrix<double, 6, 6>;
using Matrix612 = Eigen::Matrix<double, 6, 12>;
/// One shared linear-tetrahedron kernel for assembly, stresses, and derivatives.
struct Element {
    double volume = 0;
    Matrix612 B = Matrix612::Zero();
    std::array<Matrix612, 12> dB;
    Eigen::Matrix<double, 12, 1> dvolume = Eigen::Matrix<double, 12, 1>::Zero();
    /// Keep unused derivative storage initialized when copying value-only kernels.
    Element() {
        for (auto &derivative : dB)
            derivative.setZero();
    }
};
/// Isotropic 3D elasticity in engineering-shear Voigt order xx,yy,zz,xy,yz,zx.
Matrix6 elasticity_matrix(double lame_lambda, double lame_mu);
/// Evaluate shape-function gradients and optionally their exact derivatives.
Element element_kernel(const Eigen::Matrix<double, 4, 3> &vertices, bool derivatives);
/// Reconstruct the symmetric stress tensor from the six Voigt components.
Eigen::Matrix3d stress_tensor(const Eigen::Matrix<double, 6, 1> &stress);
/// FEM state for one immutable geometry revision.
struct FEM {
    Sparse K, N;
    Eigen::MatrixXd rigid, rigid_orthonormal;
    Eigen::MatrixXd rigid_gram_inverse;
    Eigen::VectorXd mass;
    Matrix6 D;
    std::vector<Element> elements;
    /// Assemble stiffness, unit-normal load map, masses, and rigid modes.
    FEM(const TetMesh &mesh, double lame_lambda, double lame_mu);
    /// Project a vector to the orthogonal complement of all rigid modes.
    Eigen::VectorXd project(const Eigen::VectorXd &vector) const;
    /// Return all six stress components per element for a displacement.
    Eigen::MatrixXd stresses(const TetMesh &mesh, const Eigen::VectorXd &displacement) const;
};
} // namespace wcso
