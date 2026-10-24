#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/QR>
#include <stdexcept>
#include <wcso/fem/Elasticity.h>

namespace wcso {
Matrix6 elasticity_matrix(double lame_lambda, double lame_mu) {
    if (!(lame_mu > 0 && 3 * lame_lambda + 2 * lame_mu > 0))
        throw std::invalid_argument("Invalid Lame parameters");
    Matrix6 D = Matrix6::Zero();
    D.topLeftCorner<3, 3>().setConstant(lame_lambda);
    for (int i = 0; i < 3; ++i)
        D(i, i) += 2 * lame_mu;
    for (int i = 3; i < 6; ++i)
        D(i, i) = lame_mu;
    return D;
}
namespace {
Matrix612 make_B(const Eigen::Matrix4d &inverse) {
    Matrix612 B = Matrix612::Zero();
    for (int i = 0; i < 4; ++i) {
        const double x = inverse(1, i), y = inverse(2, i), z = inverse(3, i);
        B(0, 3 * i) = x;
        B(1, 3 * i + 1) = y;
        B(2, 3 * i + 2) = z;
        B(3, 3 * i) = y;
        B(3, 3 * i + 1) = x;
        B(4, 3 * i + 1) = z;
        B(4, 3 * i + 2) = y;
        B(5, 3 * i) = z;
        B(5, 3 * i + 2) = x;
    }
    return B;
}
} // namespace
Element element_kernel(const Eigen::Matrix<double, 4, 3> &x, bool derivatives) {
    Element e;
    e.volume = tet_volume(x);
    if (!(e.volume > 1e-20))
        throw std::runtime_error("Inverted or degenerate tetrahedron");
    Eigen::Matrix4d M;
    M.col(0).setOnes();
    M.rightCols<3>() = x;
    const Eigen::Matrix4d inverse = M.inverse();
    e.B = make_B(inverse);
    if (derivatives) {
        e.dvolume = tet_volume_gradient(x);
        for (int i = 0; i < 4; ++i)
            for (int d = 0; d < 3; ++d)
                e.dB[3 * i + d] = make_B(-inverse.col(i) * inverse.row(d + 1));
    }
    return e;
}
Eigen::Matrix3d stress_tensor(const Eigen::Matrix<double, 6, 1> &s) {
    Eigen::Matrix3d tensor;
    tensor << s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2];
    return tensor;
}
FEM::FEM(const TetMesh &mesh, double lame_lambda, double lame_mu)
    : D(elasticity_matrix(lame_lambda, lame_mu)) {
    const int nv = static_cast<int>(mesh.vertices.rows()), nt = static_cast<int>(mesh.tets.rows());
    const int ns = static_cast<int>(mesh.surface_to_volume.size()), nr = 6 * mesh.component_count;
    if (!mesh.positive_volumes())
        throw std::runtime_error("Cannot assemble inverted volume mesh");
    elements.resize(static_cast<size_t>(nt));
    std::vector<std::array<double, 144>> stiffness(static_cast<size_t>(nt));
    // Each worker writes a disjoint element; sparse accumulation is deterministic.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int t = 0; t < nt; ++t) {
        elements[t] = element_kernel(tet_vertices(mesh, t), false);
        const auto &e = elements[t];
        const Eigen::Matrix<double, 12, 12> local = e.volume * e.B.transpose() * D * e.B;
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j)
                stiffness[t][12 * i + j] = local(i, j);
    }
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(static_cast<size_t>(nt) * 144);
    for (int t = 0; t < nt; ++t)
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j)
                entries.emplace_back(3 * mesh.tets(t, i / 3) + i % 3, 3 * mesh.tets(t, j / 3) + j % 3,
                                     stiffness[t][12 * i + j]);
    K.resize(3 * nv, 3 * nv);
    K.setFromTriplets(entries.begin(), entries.end());
    K.makeCompressed();
    const auto s = mesh.surface();
    Eigen::MatrixXd normals = Eigen::MatrixXd::Zero(ns, 3);
    for (int f = 0; f < s.faces.rows(); ++f) {
        const Eigen::Vector3d a = s.vertices.row(s.faces(f, 1)) - s.vertices.row(s.faces(f, 0));
        const Eigen::Vector3d b = s.vertices.row(s.faces(f, 2)) - s.vertices.row(s.faces(f, 0));
        const Eigen::Vector3d normal = a.cross(b);
        for (int k = 0; k < 3; ++k)
            normals.row(s.faces(f, k)) += normal.transpose();
    }
    entries.clear();
    for (int i = 0; i < ns; ++i) {
        const double norm = normals.row(i).norm();
        if (!(norm > 1e-16))
            throw std::runtime_error("Undefined surface vertex normal");
        normals.row(i) /= norm;
        for (int d = 0; d < 3; ++d)
            entries.emplace_back(3 * mesh.surface_to_volume[i] + d, i, normals(i, d));
    }
    N.resize(3 * nv, ns);
    N.setFromTriplets(entries.begin(), entries.end());
    mass = mesh.surface_areas();
    if (mass.minCoeff() <= 0)
        throw std::runtime_error("Nonpositive surface mass");
    rigid = Eigen::MatrixXd::Zero(3 * nv, nr);
    for (int i = 0; i < nv; ++i) {
        const int c = 6 * mesh.component[i];
        rigid.block<3, 3>(3 * i, c).setIdentity();
        const Eigen::Vector3d x = mesh.vertices.row(i);
        for (int d = 0; d < 3; ++d)
            rigid.block<3, 1>(3 * i, c + 3 + d) = Eigen::Vector3d::Unit(d).cross(x);
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(rigid);
    if (qr.rank() != nr)
        throw std::runtime_error("Degenerate rigid-mode basis");
    rigid_orthonormal = qr.householderQ() * Eigen::MatrixXd::Identity(3 * nv, nr);
    const Eigen::MatrixXd gram = rigid.transpose() * rigid;
    rigid_gram_inverse = gram.ldlt().solve(Eigen::MatrixXd::Identity(nr, nr));
}
Eigen::VectorXd FEM::project(const Eigen::VectorXd &v) const {
    return v - rigid_orthonormal * (rigid_orthonormal.transpose() * v);
}
Eigen::MatrixXd FEM::stresses(const TetMesh &mesh, const Eigen::VectorXd &u) const {
    Eigen::MatrixXd result(mesh.tets.rows(), 6);
    for (int t = 0; t < mesh.tets.rows(); ++t) {
        Eigen::Matrix<double, 12, 1> ue;
        for (int i = 0; i < 4; ++i)
            ue.segment<3>(3 * i) = u.segment<3>(3 * mesh.tets(t, i));
        result.row(t) = (D * elements[t].B * ue).transpose();
    }
    return result;
}
} // namespace wcso
