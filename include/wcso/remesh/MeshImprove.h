#pragma once
#include <wcso/core/Config.h>
#include <wcso/mesh/TetMesh.h>

namespace wcso {
/// Try MMG in memory; update only after geometry and quality checks pass.
bool improve_mesh(TetMesh &mesh, const Config &config, std::string &diagnostic);
/// Tetrahedralize a validated normalized surface using the optional fTetWild API.
TetMesh tetrahedralize(const Surface &surface, const Config &config);
} // namespace wcso
