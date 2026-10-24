#pragma once
#include <Eigen/Core>
#include <cstdint>
#include <string>
#include <vector>

namespace wcso {
/// Affine coordinate transform; stress remains in normalized analysis units.
struct Normalize {
    Eigen::Vector3d offset = Eigen::Vector3d::Zero();
    double scale = 1;
    /// Map original vertices to coordinates with maximum bounding-box side one.
    Eigen::MatrixXd apply(const Eigen::MatrixXd &vertices) const;
    /// Recover the original coordinate system.
    Eigen::MatrixXd inverse(const Eigen::MatrixXd &vertices) const;
    /// Compute a transform, rejecting zero extent or non-finite input.
    static Normalize from_vertices(const Eigen::MatrixXd &vertices);
};
/// A consistently oriented surface mesh.
struct Surface {
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi faces;
};
/// Volume mesh with deterministic boundary maps and revision-based cache keys.
class TetMesh {
  public:
    Eigen::MatrixXd vertices;
    Eigen::MatrixXd rest_vertices;
    Eigen::MatrixXi tets, faces, edges;
    std::vector<int> surface_to_volume, volume_to_surface, boundary_tets, component;
    int component_count = 0;
    std::uint64_t geometry_version = 0, topology_version = 0;
    /// Validate connectivity, orient tets, and extract outward boundary faces.
    void rebuild_topology();
    /// Update coordinates without changing connectivity.
    void set_vertices(const Eigen::MatrixXd &updated);
    /// Return the current boundary mesh in deterministic vertex order.
    Surface surface() const;
    /// Area-lumped vertex masses: one third of incident triangle areas.
    Eigen::VectorXd surface_areas() const;
    /// Smallest surface edge length, in normalized length units.
    double min_edge_length() const;
    /// Whether every signed tetrahedron volume is strictly positive.
    bool positive_volumes() const;
};
/// Signed volume; ordering agrees with the FEM shape-function matrix.
double tet_volume(const Eigen::Matrix<double, 4, 3> &vertices);
/// All twelve derivatives of signed tetrahedron volume.
Eigen::Matrix<double, 12, 1> tet_volume_gradient(const Eigen::Matrix<double, 4, 3> &vertices);
/// Paper Sec. 6 quality, without silently renormalizing its constant.
double tet_quality(const Eigen::Matrix<double, 4, 3> &vertices);
/// Minimum paper quality over all elements; does not mutate the mesh.
double worst_quality(const TetMesh &mesh);
/// Gather one tetrahedron's vertices.
Eigen::Matrix<double, 4, 3> tet_vertices(const TetMesh &mesh, int index);
/// Check a triangle surface for finite, nondegenerate, manifold, closed geometry.
void validate_surface(const Surface &surface);
} // namespace wcso
