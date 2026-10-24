#include <stdexcept>
#include <wcso/remesh/MeshImprove.h>
#ifdef WCSO_WITH_TETWILD
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/Logger.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/TriangleInsertion.h>
#include <geogram/basic/common.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/process.h>
#include <oneapi/tbb/global_control.h>
#endif

namespace wcso {
TetMesh tetrahedralize(const Surface &surface, const Config &config) {
#ifdef WCSO_WITH_TETWILD
    GEO::initialize();
    GEO::Numeric::random_reset();
    struct ThreadLimit {
        GEO::index_t previous = GEO::Process::max_threads();
        explicit ThreadLimit(int count) {
            GEO::Process::set_max_threads(static_cast<GEO::index_t>(count));
        }
        ~ThreadLimit() {
            GEO::Process::set_max_threads(previous);
        }
    } geogram_threads(config.threads);
    oneapi::tbb::global_control tbb_threads(oneapi::tbb::global_control::max_allowed_parallelism,
                                            static_cast<size_t>(config.threads));
    floatTetWild::Logger::init(false, "");
    std::vector<floatTetWild::Vector3> vertices;
    std::vector<floatTetWild::Vector3i> faces;
    for (int i = 0; i < surface.vertices.rows(); ++i)
        vertices.push_back(surface.vertices.row(i));
    for (int i = 0; i < surface.faces.rows(); ++i)
        faces.push_back(surface.faces.row(i));
    GEO::Mesh geometry;
    std::vector<int> tags;
    floatTetWild::MeshIO::load_mesh(vertices, faces, geometry, tags);
    floatTetWild::AABBWrapper tree(geometry);
    floatTetWild::Mesh mesh;
    auto &parameters = mesh.params;
    parameters.eps_rel = config.epsilon_d / tree.get_sf_diag();
    parameters.ideal_edge_length_rel = config.tet_edge_length / tree.get_sf_diag();
    parameters.is_quiet = true;
    // Upstream otherwise shuffles input faces with std::random_device.
    parameters.not_sort_input = true;
    parameters.num_threads = static_cast<unsigned int>(config.threads);
    if (!parameters.init(tree.get_sf_diag()))
        throw std::runtime_error("fTetWild parameter initialization failed");
    tags.resize(faces.size(), 0);
    floatTetWild::simplify(vertices, faces, tags, tree, parameters, true);
    tree.init_b_mesh_and_tree(vertices, faces, mesh);
    std::vector<bool> inserted(faces.size(), false);
    floatTetWild::FloatTetDelaunay::tetrahedralize(vertices, faces, tree, mesh, inserted);
    floatTetWild::insert_triangles(vertices, faces, tags, mesh, inserted, tree, false);
    floatTetWild::optimization(vertices, faces, tags, inserted, mesh, tree, {{1, 1, 1, 1}});
    floatTetWild::correct_tracked_surface_orientation(mesh, tree);
    floatTetWild::filter_outside(mesh);
    TetMesh result;
    floatTetWild::MatrixXs V;
    // Winding-number filtering marks is_removed, not is_outside.
    floatTetWild::MeshIO::extract_volume_mesh(mesh, V, result.tets, false);
    result.vertices = V.cast<double>();
    result.rebuild_topology();
    return result;
#else
    (void)surface;
    (void)config;
    throw std::runtime_error(
        "Surface input requires WCSO_WITH_TETWILD=ON; this build accepts .mesh/.msh input");
#endif
}
} // namespace wcso
