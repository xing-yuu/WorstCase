#pragma once
#include <filesystem>
#include <wcso/mesh/TetMesh.h>

namespace wcso {
/// Read OBJ/PLY/STL. Closedness and intersection checks occur before meshing.
Surface read_surface(const std::filesystem::path &path);
/// Read Medit .mesh or ASCII Gmsh 2.2/4.1 linear tetrahedra.
TetMesh read_tet_mesh(const std::filesystem::path &path);
/// Write a Medit volume mesh, retaining full double precision.
void write_tet_mesh(const std::filesystem::path &path, const TetMesh &mesh);
/// Write a triangle OBJ, retaining full double precision.
void write_surface(const std::filesystem::path &path, const Surface &surface);
} // namespace wcso
