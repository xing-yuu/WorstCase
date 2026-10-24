#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <wcso/mesh/TetMesh.h>

namespace wcso {
Eigen::MatrixXd Normalize::apply(const Eigen::MatrixXd &v) const {
    return (v.rowwise() - offset.transpose()) / scale;
}
Eigen::MatrixXd Normalize::inverse(const Eigen::MatrixXd &v) const {
    return (v * scale).rowwise() + offset.transpose();
}
Normalize Normalize::from_vertices(const Eigen::MatrixXd &v) {
    if (v.rows() < 4 || v.cols() != 3 || !v.allFinite())
        throw std::invalid_argument("Input requires at least four finite 3D vertices");
    Normalize n;
    n.offset = v.colwise().minCoeff();
    n.scale = (v.colwise().maxCoeff() - v.colwise().minCoeff()).maxCoeff();
    if (!(n.scale > 0))
        throw std::invalid_argument("Input has zero extent");
    return n;
}
double tet_volume(const Eigen::Matrix<double, 4, 3> &x) {
    const Eigen::Vector3d a = x.row(1) - x.row(0), b = x.row(2) - x.row(0), c = x.row(3) - x.row(0);
    return a.dot(b.cross(c)) / 6;
}
Eigen::Matrix<double, 12, 1> tet_volume_gradient(const Eigen::Matrix<double, 4, 3> &x) {
    const Eigen::Vector3d a = x.row(1) - x.row(0), b = x.row(2) - x.row(0), c = x.row(3) - x.row(0);
    Eigen::Matrix<double, 12, 1> g;
    g.segment<3>(3) = b.cross(c) / 6;
    g.segment<3>(6) = c.cross(a) / 6;
    g.segment<3>(9) = a.cross(b) / 6;
    g.head<3>() = -g.segment<3>(3) - g.segment<3>(6) - g.segment<3>(9);
    return g;
}
Eigen::Matrix<double, 4, 3> tet_vertices(const TetMesh &mesh, int index) {
    Eigen::Matrix<double, 4, 3> x;
    for (int i = 0; i < 4; ++i)
        x.row(i) = mesh.vertices.row(mesh.tets(index, i));
    return x;
}
double tet_quality(const Eigen::Matrix<double, 4, 3> &x) {
    double sum = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            sum += (x.row(i) - x.row(j)).squaredNorm();
    return sum > 0 ? 144 * std::sqrt(3.0) * tet_volume(x) / std::pow(sum, 1.5) : 0;
}
double worst_quality(const TetMesh &mesh) {
    double value = std::numeric_limits<double>::infinity();
    for (int i = 0; i < mesh.tets.rows(); ++i)
        value = std::min(value, tet_quality(tet_vertices(mesh, i)));
    return value;
}
void validate_surface(const Surface &s) {
    if (s.vertices.cols() != 3 || s.faces.cols() != 3 || s.faces.rows() < 4 || !s.vertices.allFinite())
        throw std::invalid_argument("Expected a finite triangulated closed surface");
    std::map<std::pair<int, int>, std::pair<int, int>> incidence;
    std::vector<std::map<int, std::vector<int>>> links(static_cast<size_t>(s.vertices.rows()));
    std::set<std::array<int, 3>> unique;
    for (int f = 0; f < s.faces.rows(); ++f) {
        std::array<int, 3> ids;
        for (int k = 0; k < 3; ++k) {
            ids[k] = s.faces(f, k);
            if (ids[k] < 0 || ids[k] >= s.vertices.rows())
                throw std::invalid_argument("Surface index out of range");
        }
        auto key = ids;
        std::sort(key.begin(), key.end());
        if (key[0] == key[1] || key[1] == key[2] || !unique.insert(key).second)
            throw std::invalid_argument("Repeated vertex or duplicate surface face");
        const Eigen::Vector3d a = s.vertices.row(ids[1]) - s.vertices.row(ids[0]);
        const Eigen::Vector3d b = s.vertices.row(ids[2]) - s.vertices.row(ids[0]);
        if (!(a.cross(b).norm() > 1e-16))
            throw std::invalid_argument("Degenerate surface triangle");
        for (int k = 0; k < 3; ++k) {
            const int i = ids[k], j = ids[(k + 1) % 3];
            auto &edge = incidence[std::minmax(i, j)];
            ++edge.first;
            edge.second += i < j ? 1 : -1;
            const int v = ids[(k + 2) % 3];
            links[v][i].push_back(j);
            links[v][j].push_back(i);
        }
    }
    for (const auto &edge : incidence)
        if (edge.second.first != 2 || edge.second.second != 0)
            throw std::invalid_argument("Surface must be watertight, manifold, and consistently oriented");
    for (const auto &link : links) {
        if (link.empty())
            throw std::invalid_argument("Surface contains unused vertices");
        std::set<int> visited;
        std::vector<int> pending{link.begin()->first};
        while (!pending.empty()) {
            const int v = pending.back();
            pending.pop_back();
            if (!visited.insert(v).second)
                continue;
            if (link.at(v).size() != 2)
                throw std::invalid_argument("Nonmanifold surface vertex");
            for (int neighbor : link.at(v))
                pending.push_back(neighbor);
        }
        if (visited.size() != link.size())
            throw std::invalid_argument("Disconnected face fans at a surface vertex");
    }
}
void TetMesh::rebuild_topology() {
    if (vertices.cols() != 3 || tets.cols() != 4 || tets.rows() == 0 || !vertices.allFinite())
        throw std::invalid_argument("Invalid volume mesh");
    struct FaceRecord {
        int count = 0, tet = -1;
        std::array<int, 3> oriented{};
    };
    std::map<std::array<int, 3>, FaceRecord> records;
    std::set<std::array<int, 4>> unique;
    std::vector<int> parent(static_cast<size_t>(vertices.rows()));
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](int v) {
        while (parent[v] != v) {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    };
    std::vector<bool> used(static_cast<size_t>(vertices.rows()), false);
    for (int t = 0; t < tets.rows(); ++t) {
        std::array<int, 4> ids;
        for (int k = 0; k < 4; ++k) {
            ids[k] = tets(t, k);
            if (ids[k] < 0 || ids[k] >= vertices.rows())
                throw std::invalid_argument("Tet index out of range");
            used[ids[k]] = true;
        }
        auto key = ids;
        std::sort(key.begin(), key.end());
        if (std::adjacent_find(key.begin(), key.end()) != key.end() || !unique.insert(key).second)
            throw std::invalid_argument("Repeated vertex or duplicate tetrahedron");
        const double vol = tet_volume(tet_vertices(*this, t));
        if (!std::isfinite(vol) || std::abs(vol) < 1e-20)
            throw std::invalid_argument("Degenerate tetrahedron");
        if (vol < 0)
            std::swap(tets(t, 0), tets(t, 1));
        for (int k = 1; k < 4; ++k)
            parent[find(tets(t, k))] = find(tets(t, 0));
        for (int opposite = 0; opposite < 4; ++opposite) {
            std::array<int, 3> f;
            int at = 0;
            for (int k = 0; k < 4; ++k)
                if (k != opposite)
                    f[at++] = tets(t, k);
            const Eigen::Vector3d a = vertices.row(f[1]) - vertices.row(f[0]);
            const Eigen::Vector3d b = vertices.row(f[2]) - vertices.row(f[0]);
            const Eigen::Vector3d inward = vertices.row(tets(t, opposite)) - vertices.row(f[0]);
            if (a.cross(b).dot(inward) > 0)
                std::swap(f[1], f[2]);
            auto sorted = f;
            std::sort(sorted.begin(), sorted.end());
            auto &record = records[sorted];
            if (record.count == 1)
                for (int k = 0; k < 3; ++k)
                    if (record.oriented[k] == f[0] && record.oriented[(k + 1) % 3] == f[1])
                        throw std::invalid_argument(
                            "Overlapping tetrahedra on the same side of a shared face");
            if (++record.count > 2)
                throw std::invalid_argument("Nonmanifold tetrahedral face");
            record.tet = t;
            record.oriented = f;
        }
    }
    if (std::find(used.begin(), used.end(), false) != used.end())
        throw std::invalid_argument("Volume mesh contains unused vertices");
    component.resize(static_cast<size_t>(vertices.rows()));
    std::map<int, int> component_ids;
    for (int i = 0; i < vertices.rows(); ++i) {
        const int root = find(i);
        if (!component_ids.count(root))
            component_ids[root] = static_cast<int>(component_ids.size());
        component[i] = component_ids[root];
    }
    component_count = static_cast<int>(component_ids.size());
    std::set<int> boundary_vertices, selected_tets;
    std::vector<std::array<int, 3>> boundary_faces;
    for (const auto &item : records)
        if (item.second.count == 1) {
            const auto &record = item.second;
            boundary_faces.push_back(record.oriented);
            selected_tets.insert(record.tet);
            for (int v : record.oriented)
                boundary_vertices.insert(v);
        }
    surface_to_volume.assign(boundary_vertices.begin(), boundary_vertices.end());
    boundary_tets.assign(selected_tets.begin(), selected_tets.end());
    volume_to_surface.assign(static_cast<size_t>(vertices.rows()), -1);
    for (int i = 0; i < static_cast<int>(surface_to_volume.size()); ++i)
        volume_to_surface[surface_to_volume[i]] = i;
    faces.resize(static_cast<int>(boundary_faces.size()), 3);
    std::set<std::pair<int, int>> edge_set;
    for (int f = 0; f < faces.rows(); ++f) {
        for (int k = 0; k < 3; ++k)
            faces(f, k) = volume_to_surface[boundary_faces[f][k]];
        for (int k = 0; k < 3; ++k)
            edge_set.insert(std::minmax(faces(f, k), faces(f, (k + 1) % 3)));
    }
    edges.resize(static_cast<int>(edge_set.size()), 2);
    int index = 0;
    for (const auto &edge : edge_set) {
        edges(index, 0) = edge.first;
        edges(index++, 1) = edge.second;
    }
    validate_surface(surface());
    rest_vertices = vertices;
    ++geometry_version;
    ++topology_version;
}
void TetMesh::set_vertices(const Eigen::MatrixXd &v) {
    if (v.rows() != vertices.rows() || v.cols() != 3 || !v.allFinite())
        throw std::invalid_argument("Invalid updated coordinates");
    vertices = v;
    ++geometry_version;
}
Surface TetMesh::surface() const {
    Surface s;
    s.faces = faces;
    s.vertices.resize(static_cast<int>(surface_to_volume.size()), 3);
    for (int i = 0; i < s.vertices.rows(); ++i)
        s.vertices.row(i) = vertices.row(surface_to_volume[i]);
    return s;
}
Eigen::VectorXd TetMesh::surface_areas() const {
    const auto s = surface();
    Eigen::VectorXd areas = Eigen::VectorXd::Zero(s.vertices.rows());
    for (int f = 0; f < faces.rows(); ++f) {
        const Eigen::Vector3d a = s.vertices.row(faces(f, 1)) - s.vertices.row(faces(f, 0));
        const Eigen::Vector3d b = s.vertices.row(faces(f, 2)) - s.vertices.row(faces(f, 0));
        const double share = a.cross(b).norm() / 6;
        for (int k = 0; k < 3; ++k)
            areas[faces(f, k)] += share;
    }
    return areas;
}
double TetMesh::min_edge_length() const {
    const auto s = surface();
    double length = std::numeric_limits<double>::infinity();
    for (int i = 0; i < edges.rows(); ++i)
        length = std::min(length, (s.vertices.row(edges(i, 0)) - s.vertices.row(edges(i, 1))).norm());
    return length;
}
bool TetMesh::positive_volumes() const {
    for (int t = 0; t < tets.rows(); ++t)
        if (!(tet_volume(tet_vertices(*this, t)) > 1e-20))
            return false;
    return true;
}
} // namespace wcso
