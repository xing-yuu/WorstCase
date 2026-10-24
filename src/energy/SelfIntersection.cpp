#include <igl/edges.h>
#include <ipc/collisions/collisions.hpp>
#include <ipc/ipc.hpp>
#include <ipc/potentials/barrier_potential.hpp>
#include <limits>
#include <wcso/energy/Energies.h>

namespace wcso {
EnergyValue self_intersection(const TetMesh &mesh, double d0, bool gradient) {
    const auto surface = mesh.surface();
    Eigen::MatrixXd rest(surface.vertices.rows(), 3);
    for (int i = 0; i < rest.rows(); ++i)
        rest.row(i) = mesh.rest_vertices.row(mesh.surface_to_volume[i]);
    ipc::CollisionMesh collision_mesh(rest, mesh.edges, surface.faces);
    ipc::Collisions collisions;
    collisions.build(collision_mesh, surface.vertices, d0);
    const ipc::BarrierPotential potential(d0);
    EnergyValue out;
    out.value = potential(collisions, collision_mesh, surface.vertices);
    if (gradient) {
        out.gradient = Eigen::VectorXd::Zero(3 * mesh.vertices.rows());
        const Eigen::VectorXd local = potential.gradient(collisions, collision_mesh, surface.vertices);
        for (int i = 0; i < surface.vertices.rows(); ++i)
            out.gradient.segment<3>(3 * mesh.surface_to_volume[i]) = local.segment<3>(3 * i);
    }
    return out;
}
bool collision_free(const Surface &start, const Eigen::MatrixXd &end_vertices) {
    Eigen::MatrixXi edges;
    igl::edges(start.faces, edges);
    ipc::CollisionMesh mesh(start.vertices, edges, start.faces);
    return ipc::is_step_collision_free(mesh, start.vertices, end_vertices);
}
bool has_intersections(const Surface &surface) {
    Eigen::MatrixXi edges;
    igl::edges(surface.faces, edges);
    ipc::CollisionMesh mesh(surface.vertices, edges, surface.faces);
    return ipc::has_intersections(mesh, surface.vertices);
}
} // namespace wcso
