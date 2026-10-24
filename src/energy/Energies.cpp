#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <igl/AABB.h>
#include <limits>
#include <map>
#include <stdexcept>
#include <wcso/energy/Energies.h>

namespace wcso {
namespace {
constexpr double infinity = std::numeric_limits<double>::infinity();
void restrict_surface(const TetMesh &mesh, Eigen::VectorXd &gradient) {
    for (int i = 0; i < mesh.vertices.rows(); ++i)
        if (mesh.volume_to_surface[i] < 0)
            gradient.segment<3>(3 * i).setZero();
}
double stress_energy(const std::vector<Eigen::Vector3d> &spectra, double s, double gamma,
                     double *first = nullptr, double *second = nullptr) {
    double energy = s, d = 1, dd = 0;
    for (const auto &eigenvalues : spectra)
        for (int j = 0; j < 3; ++j) {
            const double minus = s - eigenvalues[j], plus = s + eigenvalues[j];
            if (!(minus > 0 && plus > 0))
                return infinity;
            energy -= gamma * (std::log(minus) + std::log(plus));
            d -= gamma * (1 / minus + 1 / plus);
            dd += gamma * (1 / (minus * minus) + 1 / (plus * plus));
        }
    if (first)
        *first = d;
    if (second)
        *second = dd;
    return energy;
}
double optimize_sigma(const std::vector<Eigen::Vector3d> &spectra, double maximum, double gamma) {
    const double margin = std::max(1e-12, 1e-12 * maximum);
    double low = maximum + margin;
    double high = maximum + 6 * gamma * static_cast<double>(spectra.size()) + 1;
    double s = (low + high) / 2;
    for (int iteration = 0; iteration < 200; ++iteration) {
        double derivative, curvature;
        stress_energy(spectra, s, gamma, &derivative, &curvature);
        if (std::abs(derivative) < 1e-12)
            return s;
        if (derivative < 0)
            low = s;
        else
            high = s;
        if (high - low < 1e-12 * std::max(1.0, std::abs(s)))
            return (low + high) / 2;
        const double candidate = s - derivative / curvature;
        s = std::isfinite(candidate) && candidate > low && candidate < high ? candidate : (low + high) / 2;
    }
    throw std::runtime_error("Stress-bound Newton solve did not converge");
}
std::array<Eigen::Vector3d, 3> area_gradient(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                                             const Eigen::Vector3d &c) {
    const Eigen::Vector3d normal = (b - a).cross(c - a).normalized();
    const Eigen::Vector3d gb = (c - a).cross(normal) / 2, gc = normal.cross(b - a) / 2;
    return {-gb - gc, gb, gc};
}
void add_surface_gradient(const TetMesh &mesh, int vertex, const Eigen::Vector3d &value, Eigen::VectorXd &g) {
    g.segment<3>(3 * mesh.surface_to_volume[vertex]) += value;
}
} // namespace
double barrier(double x, double d, double *derivative) {
    if (derivative)
        *derivative = 0;
    if (!(x > 0))
        return infinity;
    if (x >= d)
        return 0;
    const double delta = x - d, logarithm = std::log(x / d);
    if (derivative)
        *derivative = -2 * delta * logarithm - delta * delta / x;
    return -delta * delta * logarithm;
}
EnergyValue element_flip(const TetMesh &mesh, double d1, bool gradient) {
    EnergyValue out;
    if (gradient)
        out.gradient = Eigen::VectorXd::Zero(3 * mesh.vertices.rows());
    for (int t = 0; t < mesh.tets.rows(); ++t) {
        const auto x = tet_vertices(mesh, t);
        double derivative = 0;
        out.value += barrier(tet_volume(x), d1, &derivative);
        if (!std::isfinite(out.value))
            return out;
        if (gradient && derivative != 0) {
            const Eigen::Matrix<double, 12, 1> local = derivative * tet_volume_gradient(x);
            for (int i = 0; i < 4; ++i)
                out.gradient.segment<3>(3 * mesh.tets(t, i)) += local.segment<3>(3 * i);
        }
    }
    if (gradient)
        restrict_surface(mesh, out.gradient);
    return out;
}
StressValue stress_barrier(const TetMesh &mesh, const FEM &fem, const Config &c,
                           const Eigen::VectorXd &pressure, bool gradient, double fixed_sigma_bar) {
    StressValue out;
    if (gradient)
        out.gradient = Eigen::VectorXd::Zero(3 * mesh.vertices.rows());
    ProjectedSolver solver(fem, c);
    const Eigen::VectorXd force = fem.N * pressure;
    out.displacement = solver.solve(force).x;
    out.stresses = fem.stresses(mesh, out.displacement);
    std::vector<Eigen::Vector3d> spectra;
    std::vector<bool> boundary(static_cast<size_t>(mesh.tets.rows()), false);
    for (int t : mesh.boundary_tets)
        boundary[t] = true;
    for (int t = 0; t < mesh.tets.rows(); ++t) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(stress_tensor(out.stresses.row(t).transpose()));
        if (eig.info() != Eigen::Success)
            throw std::runtime_error("Stress tensor eigendecomposition failed");
        const double maximum = eig.eigenvalues().cwiseAbs().maxCoeff();
        out.max_all = std::max(out.max_all, maximum);
        if (boundary[t]) {
            spectra.push_back(eig.eigenvalues());
            out.max_boundary = std::max(out.max_boundary, maximum);
        }
    }
    if (spectra.empty())
        throw std::runtime_error("No boundary tetrahedra");
    out.sigma_bar =
        fixed_sigma_bar > 0 ? fixed_sigma_bar : optimize_sigma(spectra, out.max_boundary, c.gamma);
    out.value = stress_energy(spectra, out.sigma_bar, c.gamma);
    if (!gradient || !std::isfinite(out.value))
        return out;
    Eigen::VectorXd gu = Eigen::VectorXd::Zero(fem.K.rows());
    std::vector<Eigen::Matrix<double, 6, 1>> h(static_cast<size_t>(mesh.tets.rows()),
                                               Eigen::Matrix<double, 6, 1>::Zero());
    for (int t : mesh.boundary_tets) {
        const Eigen::Matrix3d sigma = stress_tensor(out.stresses.row(t).transpose());
        const Eigen::Matrix3d minus = out.sigma_bar * Eigen::Matrix3d::Identity() - sigma;
        const Eigen::Matrix3d plus = out.sigma_bar * Eigen::Matrix3d::Identity() + sigma;
        const Eigen::Matrix3d H = c.gamma * (minus.llt().solve(Eigen::Matrix3d::Identity()) -
                                             plus.llt().solve(Eigen::Matrix3d::Identity()));
        h[t] << H(0, 0), H(1, 1), H(2, 2), 2 * H(0, 1), 2 * H(1, 2), 2 * H(2, 0);
        const Eigen::Matrix<double, 12, 1> local = fem.elements[t].B.transpose() * fem.D * h[t];
        for (int i = 0; i < 4; ++i)
            gu.segment<3>(3 * mesh.tets(t, i)) += local.segment<3>(3 * i);
    }
    const Eigen::VectorXd adjoint = solver.solve(gu).x;
    for (int t = 0; t < mesh.tets.rows(); ++t) {
        const Element e = element_kernel(tet_vertices(mesh, t), true);
        Eigen::Matrix<double, 12, 1> u, a;
        for (int i = 0; i < 4; ++i) {
            u.segment<3>(3 * i) = out.displacement.segment<3>(3 * mesh.tets(t, i));
            a.segment<3>(3 * i) = adjoint.segment<3>(3 * mesh.tets(t, i));
        }
        const Eigen::Matrix<double, 6, 1> bu = e.B * u, ba = e.B * a;
        for (int k = 0; k < 12; ++k) {
            const int id = mesh.tets(t, k / 3);
            if (mesh.volume_to_surface[id] < 0)
                continue;
            const auto &db = e.dB[k];
            const double direct = h[t].dot(fem.D * (db * u));
            const double stiffness = e.dvolume[k] * ba.dot(fem.D * bu) +
                                     e.volume * ((db * a).dot(fem.D * bu) + ba.dot(fem.D * (db * u)));
            out.gradient[3 * id + k % 3] += direct - stiffness;
        }
    }
    const auto surface = mesh.surface();
    Eigen::MatrixXd normals = Eigen::MatrixXd::Zero(surface.vertices.rows(), 3);
    for (int t = 0; t < surface.faces.rows(); ++t) {
        const auto &F = surface.faces;
        const Eigen::Vector3d e1 = surface.vertices.row(F(t, 1)) - surface.vertices.row(F(t, 0));
        const Eigen::Vector3d e2 = surface.vertices.row(F(t, 2)) - surface.vertices.row(F(t, 0));
        for (int k = 0; k < 3; ++k)
            normals.row(F(t, k)) += e1.cross(e2).transpose();
    }
    Eigen::MatrixXd weights(normals.rows(), 3);
    for (int i = 0; i < normals.rows(); ++i) {
        const double norm = normals.row(i).norm();
        const Eigen::Vector3d n = normals.row(i) / norm;
        const Eigen::Vector3d a = adjoint.segment<3>(3 * mesh.surface_to_volume[i]);
        weights.row(i) = (pressure[i] * (a - n * n.dot(a)) / norm).transpose();
    }
    for (int t = 0; t < surface.faces.rows(); ++t) {
        const int a = surface.faces(t, 0), b = surface.faces(t, 1), d = surface.faces(t, 2);
        const Eigen::Vector3d w = weights.row(a) + weights.row(b) + weights.row(d);
        const Eigen::Vector3d e1 = surface.vertices.row(b) - surface.vertices.row(a);
        const Eigen::Vector3d e2 = surface.vertices.row(d) - surface.vertices.row(a);
        const Eigen::Vector3d gb = e2.cross(w), gd = w.cross(e1);
        add_surface_gradient(mesh, a, -gb - gd, out.gradient);
        add_surface_gradient(mesh, b, gb, out.gradient);
        add_surface_gradient(mesh, d, gd, out.gradient);
    }
    // Differentiate the moving null-space constraint through the augmented solve.
    const Eigen::VectorXd load_multiplier = fem.rigid_gram_inverse * (fem.rigid.transpose() * force);
    const Eigen::VectorXd adjoint_multiplier = fem.rigid_gram_inverse * (fem.rigid.transpose() * gu);
    for (int i : mesh.surface_to_volume) {
        const int index = 6 * mesh.component[i] + 3;
        out.gradient.segment<3>(3 * i) +=
            load_multiplier.segment<3>(index).cross(adjoint.segment<3>(3 * i)) +
            adjoint_multiplier.segment<3>(index).cross(out.displacement.segment<3>(3 * i));
    }
    restrict_surface(mesh, out.gradient);
    return out;
}
ICPContext make_samples(const Surface &s, double density, std::mt19937_64 &rng) {
    ICPContext context;
    context.samples_per_face.resize(static_cast<size_t>(s.faces.rows()));
    struct EdgeData {
        int a, b, c, d = -1;
    };
    std::map<std::pair<int, int>, EdgeData> edges;
    for (int f = 0; f < s.faces.rows(); ++f) {
        const Eigen::Vector3d a = s.vertices.row(s.faces(f, 0));
        const Eigen::Vector3d b = s.vertices.row(s.faces(f, 1));
        const Eigen::Vector3d c = s.vertices.row(s.faces(f, 2));
        const double area = (b - a).cross(c - a).norm() / 2;
        const double count = std::max(1.0, std::ceil(density * area));
        if (count > 1e6 || context.samples.size() + static_cast<size_t>(count) > 10000000)
            throw std::invalid_argument("Sampling density would exceed the sample budget");
        const int n = static_cast<int>(count);
        context.samples_per_face[f] = n;
        for (int j = 0; j < n; ++j) {
            const double r = std::sqrt(std::generate_canonical<double, 53>(rng));
            const double q = std::generate_canonical<double, 53>(rng);
            context.samples.push_back({f, Eigen::Vector3d(1 - r, r * (1 - q), r * q)});
        }
        for (int k = 0; k < 3; ++k) {
            const int i = s.faces(f, k), j = s.faces(f, (k + 1) % 3), opposite = s.faces(f, (k + 2) % 3);
            const std::pair<int, int> key = std::minmax(i, j);
            if (!edges.count(key))
                edges[key] = {i, j, opposite, -1};
            else {
                auto &e = edges.at(key);
                if (e.a != j || e.b != i || e.d >= 0)
                    throw std::runtime_error("Invalid oriented hinge");
                e.d = opposite;
            }
        }
    }
    for (const auto &item : edges) {
        const auto &e = item.second;
        if (e.d < 0)
            throw std::runtime_error("Open surface hinge");
        context.hinges.push_back({{e.a, e.b, e.c, e.d}});
    }
    return context;
}
void update_correspondences(ICPContext &context, const Surface &s, const Surface &target) {
    igl::AABB<Eigen::MatrixXd, 3> tree;
    tree.init(target.vertices, target.faces);
    Eigen::MatrixXd points(static_cast<int>(context.samples.size() + context.hinges.size()), 3);
    int row = 0;
    for (const auto &sample : context.samples) {
        points.row(row).setZero();
        for (int k = 0; k < 3; ++k)
            points.row(row) += sample.barycentric[k] * s.vertices.row(s.faces(sample.face, k));
        ++row;
    }
    for (const auto &h : context.hinges)
        points.row(row++) = (s.vertices.row(h.vertices[0]) + s.vertices.row(h.vertices[1])) / 2;
    Eigen::VectorXd distances;
    Eigen::VectorXi faces;
    Eigen::MatrixXd closest;
    tree.squared_distance(target.vertices, target.faces, points, distances, faces, closest);
    context.sample_targets = closest.topRows(static_cast<int>(context.samples.size()));
    context.hinge_targets = closest.bottomRows(static_cast<int>(context.hinges.size()));
}
EnergyValue shape_fidelity(const TetMesh &mesh, const ICPContext &context, bool gradient) {
    EnergyValue out;
    const auto s = mesh.surface();
    if (gradient)
        out.gradient = Eigen::VectorXd::Zero(3 * mesh.vertices.rows());
    const int nf = static_cast<int>(s.faces.rows());
    std::vector<double> averages(static_cast<size_t>(nf), 0);
    Eigen::MatrixXd face_gradients = Eigen::MatrixXd::Zero(nf, 9);
    for (int j = 0; j < static_cast<int>(context.samples.size()); ++j) {
        const auto &sample = context.samples[j];
        Eigen::Vector3d x = Eigen::Vector3d::Zero();
        for (int k = 0; k < 3; ++k)
            x += sample.barycentric[k] * s.vertices.row(s.faces(sample.face, k)).transpose();
        const Eigen::Vector3d delta = x - context.sample_targets.row(j).transpose();
        const double inverse_count = 1.0 / context.samples_per_face[sample.face];
        averages[sample.face] += delta.squaredNorm() * inverse_count;
        if (gradient)
            for (int k = 0; k < 3; ++k)
                face_gradients.block<1, 3>(sample.face, 3 * k) +=
                    (2 * inverse_count * sample.barycentric[k] * delta).transpose();
    }
    double area = 0, numerator = 0;
    Eigen::VectorXd darea, dnumerator;
    if (gradient) {
        darea = out.gradient;
        dnumerator = out.gradient;
    }
    for (int f = 0; f < nf; ++f) {
        const Eigen::Vector3d a = s.vertices.row(s.faces(f, 0)), b = s.vertices.row(s.faces(f, 1)),
                              c = s.vertices.row(s.faces(f, 2));
        const double local_area = (b - a).cross(c - a).norm() / 2;
        area += local_area;
        numerator += local_area * averages[f];
        if (gradient) {
            const auto ga = area_gradient(a, b, c);
            for (int k = 0; k < 3; ++k) {
                add_surface_gradient(mesh, s.faces(f, k), ga[k], darea);
                add_surface_gradient(mesh, s.faces(f, k),
                                     averages[f] * ga[k] +
                                         local_area * face_gradients.block<1, 3>(f, 3 * k).transpose(),
                                     dnumerator);
            }
        }
    }
    out.value = numerator / area;
    if (gradient)
        out.gradient = dnumerator / area - numerator / (area * area) * darea;
    return out;
}
EnergyValue local_smoothing(const TetMesh &mesh, const ICPContext &context, double d2, bool gradient) {
    EnergyValue out;
    if (gradient)
        out.gradient = Eigen::VectorXd::Zero(3 * mesh.vertices.rows());
    const auto s = mesh.surface();
    for (int j = 0; j < static_cast<int>(context.hinges.size()); ++j) {
        const auto &ids = context.hinges[j].vertices;
        std::array<Eigen::Vector3d, 4> x;
        for (int k = 0; k < 4; ++k)
            x[k] = s.vertices.row(ids[k]);
        const Eigen::Vector3d edge = x[1] - x[0], side1 = x[2] - x[0], side2 = x[3] - x[1];
        const Eigen::Vector3d cross1 = edge.cross(side1), cross2 = (-edge).cross(side2);
        const double l = edge.norm(), c1 = cross1.norm(), c2 = cross2.norm();
        if (std::min({l, c1, c2}) < 1e-16) {
            out.value = infinity;
            return out;
        }
        const Eigen::Vector3d en = edge / l, n1 = cross1 / c1, n2 = cross2 / c2;
        const double sin_angle = en.dot(n1.cross(n2)), cos_angle = n1.dot(n2);
        const double angle = std::atan2(sin_angle, cos_angle), area_sum = (c1 + c2) / 2;
        const double weight = 3 * l * l / area_sum;
        const Eigen::Vector3d offset = (x[0] + x[1]) / 2 - context.hinge_targets.row(j).transpose();
        const double distance = offset.norm(), ratio = distance / d2;
        const double activation = distance >= d2 ? 1 : ratio * ratio * (3 - 2 * ratio);
        const double activation_derivative = distance > 0 && distance < d2 ? 6 * ratio * (1 - ratio) / d2 : 0;
        const double shell = angle * angle * weight;
        out.value += activation * shell;
        if (!gradient || (activation == 0 && activation_derivative == 0))
            continue;
        for (int k = 0; k < 4; ++k)
            for (int axis = 0; axis < 3; ++axis) {
                std::array<Eigen::Vector3d, 4> dx;
                for (auto &v : dx)
                    v.setZero();
                dx[k][axis] = 1;
                const Eigen::Vector3d de = dx[1] - dx[0];
                const Eigen::Vector3d dc1 = de.cross(side1) + edge.cross(dx[2] - dx[0]);
                const Eigen::Vector3d dc2 = (-de).cross(side2) + (-edge).cross(dx[3] - dx[1]);
                const double dl = en.dot(de), da = (n1.dot(dc1) + n2.dot(dc2)) / 2;
                const Eigen::Vector3d dn1 = (dc1 - n1 * n1.dot(dc1)) / c1,
                                      dn2 = (dc2 - n2 * n2.dot(dc2)) / c2;
                const Eigen::Vector3d den = (de - en * dl) / l;
                const double dsin = den.dot(n1.cross(n2)) + en.dot(dn1.cross(n2) + n1.cross(dn2));
                const double dcos = dn1.dot(n2) + n1.dot(dn2);
                const double dangle =
                    (cos_angle * dsin - sin_angle * dcos) / (sin_angle * sin_angle + cos_angle * cos_angle);
                const double dweight = 6 * l * dl / area_sum - weight * da / area_sum;
                const double ddistance = distance > 0 ? offset.dot((dx[0] + dx[1]) / 2) / distance : 0;
                const double value = activation * (2 * angle * dangle * weight + angle * angle * dweight) +
                                     activation_derivative * ddistance * shell;
                out.gradient[3 * mesh.surface_to_volume[ids[k]] + axis] += value;
            }
    }
    return out;
}
double surface_deviation(const Surface &a, const Surface &b) {
    auto directed = [](const Surface &source, const Surface &target) {
        igl::AABB<Eigen::MatrixXd, 3> tree;
        tree.init(target.vertices, target.faces);
        Eigen::MatrixXd queries(source.vertices.rows() + 4 * source.faces.rows(), 3);
        queries.topRows(source.vertices.rows()) = source.vertices;
        int row = static_cast<int>(source.vertices.rows());
        for (int f = 0; f < source.faces.rows(); ++f) {
            queries.row(row++) =
                (source.vertices.row(source.faces(f, 0)) + source.vertices.row(source.faces(f, 1)) +
                 source.vertices.row(source.faces(f, 2))) /
                3;
            for (int k = 0; k < 3; ++k)
                queries.row(row++) = (source.vertices.row(source.faces(f, k)) +
                                      source.vertices.row(source.faces(f, (k + 1) % 3))) /
                                     2;
        }
        Eigen::VectorXd distances;
        Eigen::VectorXi indices;
        Eigen::MatrixXd closest;
        tree.squared_distance(target.vertices, target.faces, queries, distances, indices, closest);
        return std::sqrt(distances.maxCoeff());
    };
    return std::max(directed(a, b), directed(b, a));
}
} // namespace wcso
