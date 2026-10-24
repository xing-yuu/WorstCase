#include "BoundaryPreservation.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <libmmg3d.h>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
#include <wcso/energy/Energies.h>
#include <wcso/remesh/MeshImprove.h>

namespace wcso {
namespace {
struct MmgMesh {
    MMG5_pMesh mesh = nullptr;
    MMG5_pSol metric = nullptr;
    MmgMesh() {
        MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
    }
    ~MmgMesh() {
        MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mesh, MMG5_ARG_ppMet, &metric, MMG5_ARG_end);
    }
};
void check(int result) {
    if (result != 1)
        throw std::runtime_error("MMG API call failed");
}
} // namespace
void detail::restore_required_boundary(const TetMesh &original, TetMesh &candidate) {
    const auto boundary = original.surface();
    const int count = static_cast<int>(original.surface_to_volume.size());
    if (candidate.surface_to_volume.size() != original.surface_to_volume.size())
        throw std::runtime_error("MMG boundary vertex count changed");
    if (count == 0 || !boundary.vertices.allFinite() || !candidate.vertices.allFinite())
        throw std::runtime_error("MMG boundary coordinates are invalid");

    // MMG scales and unscales even required vertices. This bound accounts only
    // for floating-point roundoff; it is independent of the geometry tolerance.
    const double scale = std::max(
        boundary.vertices.cwiseAbs().maxCoeff(),
        (boundary.vertices.colwise().maxCoeff() - boundary.vertices.colwise().minCoeff()).maxCoeff());
    if (!(scale > 0) || !std::isfinite(scale))
        throw std::runtime_error("MMG boundary coordinate scale is invalid");
    const double roundoff = 32 * std::numeric_limits<double>::epsilon() * scale;
    using Point = std::array<double, 3>;
    std::map<Point, int> source;
    for (int i = 0; i < count; ++i) {
        const Point point{boundary.vertices(i, 0), boundary.vertices(i, 1), boundary.vertices(i, 2)};
        if (!source.emplace(point, i).second)
            throw std::runtime_error("MMG boundary vertex mapping is ambiguous");
    }

    std::vector<int> mapping(count, -1);
    std::vector<bool> used(count, false);
    const double infinity = std::numeric_limits<double>::infinity();
    for (int i = 0; i < count; ++i) {
        const Eigen::Vector3d point = candidate.vertices.row(candidate.surface_to_volume[i]).transpose();
        auto entry = source.lower_bound(Point{point[0] - roundoff, -infinity, -infinity});
        for (; entry != source.end() && entry->first[0] <= point[0] + roundoff; ++entry) {
            if (std::abs(entry->first[0] - point[0]) > roundoff ||
                std::abs(entry->first[1] - point[1]) > roundoff ||
                std::abs(entry->first[2] - point[2]) > roundoff)
                continue;
            if (mapping[i] >= 0)
                throw std::runtime_error("MMG boundary vertex mapping is ambiguous");
            mapping[i] = entry->second;
        }
        if (mapping[i] < 0)
            throw std::runtime_error("MMG boundary vertex moved beyond floating-point roundoff");
        if (used[mapping[i]])
            throw std::runtime_error("MMG boundary vertex mapping is not bijective");
        used[mapping[i]] = true;
    }

    using Face = std::array<int, 3>;
    auto canonical_faces = [](const Eigen::MatrixXi &faces, const std::vector<int> *indices) {
        std::multiset<Face> canonical;
        for (int f = 0; f < faces.rows(); ++f) {
            Face face;
            for (int k = 0; k < 3; ++k)
                face[k] = indices ? (*indices)[faces(f, k)] : faces(f, k);
            std::sort(face.begin(), face.end());
            canonical.insert(face);
        }
        return canonical;
    };
    if (canonical_faces(original.faces, nullptr) != canonical_faces(candidate.faces, &mapping))
        throw std::runtime_error("MMG boundary triangle topology changed");

    // Restore exact coordinates only after the complete correspondence passes.
    for (int i = 0; i < count; ++i)
        candidate.vertices.row(candidate.surface_to_volume[i]) = boundary.vertices.row(mapping[i]);
}
bool improve_mesh(TetMesh &mesh, const Config &config, std::string &diagnostic) {
    try {
        MmgMesh mmg;
        check(MMG3D_Set_meshSize(mmg.mesh, static_cast<MMG5_int>(mesh.vertices.rows()),
                                 static_cast<MMG5_int>(mesh.tets.rows()), 0,
                                 static_cast<MMG5_int>(mesh.faces.rows()), 0, 0));
        for (int i = 0; i < mesh.vertices.rows(); ++i) {
            check(MMG3D_Set_vertex(mmg.mesh, mesh.vertices(i, 0), mesh.vertices(i, 1), mesh.vertices(i, 2), 0,
                                   i + 1));
            if (mesh.volume_to_surface[i] >= 0)
                check(MMG3D_Set_requiredVertex(mmg.mesh, i + 1));
        }
        for (int t = 0; t < mesh.tets.rows(); ++t)
            check(MMG3D_Set_tetrahedron(mmg.mesh, mesh.tets(t, 0) + 1, mesh.tets(t, 1) + 1,
                                        mesh.tets(t, 2) + 1, mesh.tets(t, 3) + 1, 0, t + 1));
        for (int f = 0; f < mesh.faces.rows(); ++f) {
            check(MMG3D_Set_triangle(mmg.mesh, mesh.surface_to_volume[mesh.faces(f, 0)] + 1,
                                     mesh.surface_to_volume[mesh.faces(f, 1)] + 1,
                                     mesh.surface_to_volume[mesh.faces(f, 2)] + 1, 0, f + 1));
            check(MMG3D_Set_requiredTriangle(mmg.mesh, f + 1));
        }
        check(MMG3D_Set_iparameter(mmg.mesh, mmg.metric, MMG3D_IPARAM_verbose, -1));
        check(MMG3D_Set_dparameter(mmg.mesh, mmg.metric, MMG3D_DPARAM_hausd, config.epsilon_d));
        check(MMG3D_Set_dparameter(mmg.mesh, mmg.metric, MMG3D_DPARAM_hmax, config.tet_edge_length));
        if (MMG3D_mmg3dlib(mmg.mesh, mmg.metric) != MMG5_SUCCESS)
            throw std::runtime_error("MMG could not improve the mesh");
        MMG5_int nv, nt;
        check(MMG3D_Get_meshSize(mmg.mesh, &nv, &nt, nullptr, nullptr, nullptr, nullptr));
        TetMesh candidate;
        candidate.vertices.resize(nv, 3);
        candidate.tets.resize(nt, 4);
        for (int i = 0; i < nv; ++i)
            check(MMG3D_Get_vertex(mmg.mesh, &candidate.vertices(i, 0), &candidate.vertices(i, 1),
                                   &candidate.vertices(i, 2), nullptr, nullptr, nullptr));
        for (int t = 0; t < nt; ++t) {
            MMG5_int a, b, c, d;
            check(MMG3D_Get_tetrahedron(mmg.mesh, &a, &b, &c, &d, nullptr, nullptr));
            candidate.tets.row(t) << a - 1, b - 1, c - 1, d - 1;
        }
        candidate.rebuild_topology();
        detail::restore_required_boundary(mesh, candidate);
        if (!candidate.positive_volumes())
            throw std::runtime_error("MMG candidate has nonpositive tetrahedron volume");
        if (has_intersections(candidate.surface()))
            throw std::runtime_error("MMG candidate surface self-intersects");
        if (worst_quality(candidate) <= worst_quality(mesh))
            throw std::runtime_error("MMG candidate did not improve worst element quality");
        candidate.rest_vertices = candidate.vertices;
        candidate.geometry_version = mesh.geometry_version + 1;
        candidate.topology_version = mesh.topology_version + 1;
        mesh = std::move(candidate);
        diagnostic = "remesh_accepted";
        return true;
    } catch (const std::exception &e) {
        diagnostic = std::string("remesh_rolled_back: ") + e.what();
        return false;
    }
}
} // namespace wcso
