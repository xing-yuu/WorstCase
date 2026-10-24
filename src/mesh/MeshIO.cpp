#include <fstream>
#include <igl/read_triangle_mesh.h>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <wcso/mesh/MeshIO.h>

namespace wcso {
Surface read_surface(const std::filesystem::path &path) {
    Surface s;
    if (!igl::read_triangle_mesh(path.string(), s.vertices, s.faces))
        throw std::runtime_error("Cannot read surface: " + path.string());
    if (path.extension() == ".stl") {
        // STL is triangle soup. Weld only exactly equal coordinates, without
        // introducing an arbitrary geometric tolerance or moving the surface.
        if (!s.vertices.allFinite())
            throw std::runtime_error("Nonfinite STL coordinates");
        std::map<std::array<double, 3>, int> ids;
        std::vector<Eigen::Vector3d> points;
        std::vector<int> mapping(static_cast<size_t>(s.vertices.rows()));
        for (int i = 0; i < s.vertices.rows(); ++i) {
            const std::array<double, 3> point{s.vertices(i, 0), s.vertices(i, 1), s.vertices(i, 2)};
            auto entry = ids.emplace(point, static_cast<int>(points.size()));
            if (entry.second)
                points.push_back(s.vertices.row(i));
            mapping[i] = entry.first->second;
        }
        for (int f = 0; f < s.faces.rows(); ++f)
            for (int k = 0; k < 3; ++k)
                s.faces(f, k) = mapping.at(s.faces(f, k));
        s.vertices.resize(static_cast<int>(points.size()), 3);
        for (int i = 0; i < s.vertices.rows(); ++i)
            s.vertices.row(i) = points[i];
    }
    return s;
}
namespace {
TetMesh read_medit(const std::filesystem::path &path) {
    std::ifstream file(path);
    if (!file)
        throw std::runtime_error("Cannot open " + path.string());
    std::ostringstream clean;
    std::string line;
    while (std::getline(file, line))
        clean << line.substr(0, line.find('#')) << '\n';
    std::istringstream in(clean.str());
    TetMesh mesh;
    std::string section;
    const std::map<std::string, int> ignored_widths = {{"Triangles", 4},
                                                       {"Edges", 3},
                                                       {"Corners", 1},
                                                       {"Ridges", 1},
                                                       {"RequiredVertices", 1},
                                                       {"RequiredEdges", 1},
                                                       {"RequiredTriangles", 1},
                                                       {"RequiredTetrahedra", 1},
                                                       {"Normals", 3},
                                                       {"NormalAtVertices", 2},
                                                       {"NormalAtTriangleVertices", 3},
                                                       {"Tangents", 3},
                                                       {"TangentAtVertices", 2},
                                                       {"TangentAtEdgeVertices", 3}};
    while (in >> section) {
        if (section == "End")
            break;
        if (section == "MeshVersionFormatted") {
            int version;
            in >> version;
            if (version < 1 || version > 3)
                throw std::runtime_error("Unsupported Medit version");
        } else if (section == "Dimension") {
            int dimension;
            in >> dimension;
            if (dimension != 3)
                throw std::runtime_error("Expected a 3D Medit mesh");
        } else {
            int count;
            if (!(in >> count) || count < 0 || count > 100000000)
                throw std::runtime_error("Invalid Medit section count");
            if (section == "Vertices") {
                mesh.vertices.resize(count, 3);
                for (int i = 0; i < count; ++i) {
                    int reference;
                    in >> mesh.vertices(i, 0) >> mesh.vertices(i, 1) >> mesh.vertices(i, 2) >> reference;
                }
            } else if (section == "Tetrahedra") {
                mesh.tets.resize(count, 4);
                for (int i = 0; i < count; ++i) {
                    int reference;
                    for (int k = 0; k < 4; ++k) {
                        in >> mesh.tets(i, k);
                        --mesh.tets(i, k);
                    }
                    in >> reference;
                }
            } else {
                const auto width = ignored_widths.find(section);
                if (width == ignored_widths.end())
                    throw std::runtime_error("Unsupported Medit section: " + section);
                double discard;
                for (int i = 0; i < count; ++i)
                    for (int k = 0; k < width->second; ++k)
                        in >> discard;
            }
        }
        if (!in)
            throw std::runtime_error("Truncated Medit input");
    }
    return mesh;
}
TetMesh read_gmsh(const std::filesystem::path &path) {
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open " + path.string());
    std::string token;
    double version = 0;
    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::array<long long, 4>> tetrahedra;
    std::map<long long, int> node_map;
    while (in >> token) {
        if (token == "$MeshFormat") {
            int binary, size;
            in >> version >> binary >> size;
            if (binary != 0 || !(version == 2.2 || version == 4.1))
                throw std::runtime_error("Supported Gmsh formats: ASCII 2.2 and 4.1 with linear tetrahedra");
        } else if (token == "$Nodes") {
            if (version == 2.2) {
                size_t count;
                in >> count;
                for (size_t i = 0; i < count; ++i) {
                    long long tag;
                    Eigen::Vector3d x;
                    in >> tag >> x[0] >> x[1] >> x[2];
                    if (!node_map.emplace(tag, static_cast<int>(vertices.size())).second)
                        throw std::runtime_error("Duplicate Gmsh node tag");
                    vertices.push_back(x);
                }
            } else if (version == 4.1) {
                size_t blocks, total;
                long long minimum, maximum;
                in >> blocks >> total >> minimum >> maximum;
                for (size_t b = 0; b < blocks; ++b) {
                    int dimension, entity, parametric;
                    size_t count;
                    in >> dimension >> entity >> parametric >> count;
                    std::vector<long long> tags(count);
                    for (auto &tag : tags)
                        in >> tag;
                    for (size_t i = 0; i < count; ++i) {
                        Eigen::Vector3d x;
                        in >> x[0] >> x[1] >> x[2];
                        if (parametric)
                            for (int j = 0; j < dimension; ++j) {
                                double discard;
                                in >> discard;
                            }
                        if (!node_map.emplace(tags[i], static_cast<int>(vertices.size())).second)
                            throw std::runtime_error("Duplicate Gmsh node tag");
                        vertices.push_back(x);
                    }
                }
                if (vertices.size() != total)
                    throw std::runtime_error("Gmsh node count mismatch");
            } else
                throw std::runtime_error("Missing Gmsh format header");
        } else if (token == "$Elements") {
            auto parse_line = [&](const std::string &line, int known_type) {
                std::istringstream row(line);
                long long tag;
                int type = known_type;
                row >> tag;
                if (known_type < 0) {
                    int tags;
                    row >> type >> tags;
                    for (int i = 0; i < tags; ++i) {
                        long long discard;
                        row >> discard;
                    }
                }
                if (type != 4 && type != 1 && type != 2 && type != 3 && type != 15)
                    throw std::runtime_error("Unsupported Gmsh element type; provide linear tetrahedra and "
                                             "linear boundary elements");
                if (type == 4) {
                    std::array<long long, 4> tet;
                    for (auto &vertex : tet)
                        if (!(row >> vertex))
                            throw std::runtime_error("Truncated Gmsh tetrahedron");
                    tetrahedra.push_back(tet);
                }
            };
            if (version == 2.2) {
                size_t count;
                in >> count;
                std::string line;
                std::getline(in, line);
                for (size_t i = 0; i < count; ++i) {
                    std::getline(in, line);
                    parse_line(line, -1);
                }
            } else if (version == 4.1) {
                size_t blocks, total;
                long long minimum, maximum;
                in >> blocks >> total >> minimum >> maximum;
                for (size_t b = 0; b < blocks; ++b) {
                    int dimension, entity, type;
                    size_t count;
                    in >> dimension >> entity >> type >> count;
                    std::string line;
                    std::getline(in, line);
                    for (size_t i = 0; i < count; ++i) {
                        std::getline(in, line);
                        parse_line(line, type);
                    }
                }
            } else
                throw std::runtime_error("Missing Gmsh format header");
        }
        if (!in && !in.eof())
            throw std::runtime_error("Malformed Gmsh input");
    }
    if (tetrahedra.empty())
        throw std::runtime_error("Gmsh input contains no linear tetrahedra");
    TetMesh mesh;
    mesh.vertices.resize(static_cast<int>(vertices.size()), 3);
    for (int i = 0; i < mesh.vertices.rows(); ++i)
        mesh.vertices.row(i) = vertices[i];
    mesh.tets.resize(static_cast<int>(tetrahedra.size()), 4);
    for (int i = 0; i < mesh.tets.rows(); ++i)
        for (int k = 0; k < 4; ++k) {
            const auto found = node_map.find(tetrahedra[i][k]);
            if (found == node_map.end())
                throw std::runtime_error("Unknown Gmsh node tag");
            mesh.tets(i, k) = found->second;
        }
    return mesh;
}
} // namespace
TetMesh read_tet_mesh(const std::filesystem::path &path) {
    TetMesh mesh;
    if (path.extension() == ".msh")
        mesh = read_gmsh(path);
    else
        mesh = read_medit(path);
    // Meshing tools may retain unused nodes after deleting exterior elements.
    std::map<int, int> used;
    for (int t = 0; t < mesh.tets.rows(); ++t)
        for (int k = 0; k < 4; ++k) {
            const int id = mesh.tets(t, k);
            if (id < 0 || id >= mesh.vertices.rows())
                throw std::runtime_error("Tet index out of range");
            used[id] = 0;
        }
    Eigen::MatrixXd compact(static_cast<int>(used.size()), 3);
    int index = 0;
    for (auto &entry : used) {
        entry.second = index;
        compact.row(index++) = mesh.vertices.row(entry.first);
    }
    for (int t = 0; t < mesh.tets.rows(); ++t)
        for (int k = 0; k < 4; ++k)
            mesh.tets(t, k) = used.at(mesh.tets(t, k));
    mesh.vertices = compact;
    mesh.rebuild_topology();
    return mesh;
}
void write_tet_mesh(const std::filesystem::path &path, const TetMesh &mesh) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot write " + path.string());
    out << std::setprecision(17) << "MeshVersionFormatted 2\nDimension 3\nVertices\n"
        << mesh.vertices.rows() << '\n';
    for (int i = 0; i < mesh.vertices.rows(); ++i)
        out << mesh.vertices.row(i) << " 0\n";
    out << "Triangles\n" << mesh.faces.rows() << '\n';
    for (int f = 0; f < mesh.faces.rows(); ++f) {
        for (int k = 0; k < 3; ++k)
            out << mesh.surface_to_volume[mesh.faces(f, k)] + 1 << ' ';
        out << "0\n";
    }
    out << "Tetrahedra\n" << mesh.tets.rows() << '\n';
    for (int t = 0; t < mesh.tets.rows(); ++t)
        out << (mesh.tets.row(t).array() + 1).matrix() << " 0\n";
    out << "End\n";
    if (!out)
        throw std::runtime_error("Failed writing " + path.string());
}
void write_surface(const std::filesystem::path &path, const Surface &surface) {
    std::ofstream out(path);
    if (!out)
        throw std::runtime_error("Cannot write " + path.string());
    out << std::setprecision(17);
    for (int i = 0; i < surface.vertices.rows(); ++i)
        out << "v " << surface.vertices.row(i) << '\n';
    for (int i = 0; i < surface.faces.rows(); ++i)
        out << "f " << (surface.faces.row(i).array() + 1).matrix() << '\n';
    if (!out)
        throw std::runtime_error("Failed writing " + path.string());
}
} // namespace wcso
