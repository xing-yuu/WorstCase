#pragma once
#include <wcso/mesh/TetMesh.h>

namespace wcso::detail {
/// Restore required boundary coordinates after MMG's floating-point round trip.
/// Reject movement beyond roundoff, ambiguous correspondence, or changed faces.
/// The candidate is unchanged if validation fails.
void restore_required_boundary(const TetMesh &original, TetMesh &candidate);
} // namespace wcso::detail
