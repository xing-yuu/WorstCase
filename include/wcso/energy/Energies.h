#pragma once
#include <functional>
#include <wcso/analysis/WorstCaseLoad.h>

namespace wcso {
/// Value and optional full-volume gradient from the same evaluation routine.
struct EnergyValue {
    double value = 0;
    Eigen::VectorXd gradient;
};
/// Stress evaluation distinguishes boundary constraints from global diagnostics.
struct StressValue : EnergyValue {
    double sigma_bar = 0, max_boundary = 0, max_all = 0;
    Eigen::VectorXd displacement;
    Eigen::MatrixXd stresses;
};
/// Eq. (7), using only tets with a boundary face for value and all derivatives.
StressValue stress_barrier(const TetMesh &mesh, const FEM &fem, const Config &config,
                           const Eigen::VectorXd &pressure, bool gradient, double fixed_sigma_bar = -1);
/// Analytic locally supported barrier P(x,d); infinity for x<=0.
double barrier(double x, double d, double *derivative = nullptr);
/// Eq. E_ef over signed element volumes, with exact derivatives.
EnergyValue element_flip(const TetMesh &mesh, double d1, bool gradient);
/// IPC collision barrier E_si and correctly scattered surface gradient.
EnergyValue self_intersection(const TetMesh &mesh, double d0, bool gradient);
/// Conservative CCD for a linear surface trajectory.
bool collision_free(const Surface &start, const Eigen::MatrixXd &end_vertices);
/// Static surface self-intersection check, including nonadjacent triangles.
bool has_intersections(const Surface &surface);
/// Fixed barycentric samples used throughout a topology epoch.
struct Sample {
    int face = 0;
    Eigen::Vector3d barycentric;
};
/// A surface-edge hinge, oriented as faces (a,b,c) and (b,a,d).
struct Hinge {
    std::array<int, 4> vertices;
};
/// Immutable ICP targets during each value/gradient/line-search comparison.
struct ICPContext {
    std::vector<Sample> samples;
    std::vector<int> samples_per_face;
    std::vector<Hinge> hinges;
    Eigen::MatrixXd sample_targets, hinge_targets;
};
/// Create uniform area samples; only changes after a topology change.
ICPContext make_samples(const Surface &surface, double density, std::mt19937_64 &rng);
/// Refresh closest points outside the line search, keeping sample identities.
void update_correspondences(ICPContext &context, const Surface &current, const Surface &target);
/// Area-normalized quadrature of squared distance with area derivatives.
EnergyValue shape_fidelity(const TetMesh &mesh, const ICPContext &context, bool gradient);
/// Zero-rest-angle hinge energy with midpoint activation and all product terms.
EnergyValue local_smoothing(const TetMesh &mesh, const ICPContext &context, double d2, bool gradient);
/// Symmetric sampled surface deviation used to validate remeshing.
double surface_deviation(const Surface &first, const Surface &second);
} // namespace wcso
