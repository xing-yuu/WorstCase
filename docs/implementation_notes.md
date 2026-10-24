# Numerical specification and paper/code decisions

This is the author-approved corrected implementation. The supplied review was
an audit of a snapshot, not an instruction to preserve known numerical defects.
The changes below intentionally alter the old trajectory. `configs/paper.json`
names the author-approved corrected settings; it is not a literal transcription
of the manuscript and does not certify reproduction of every published
experiment. Author approval to retain a difference does not make that
difference disappear from the manuscript.

## Section 2 decisions

| Item | Implemented decision |
|---|---|
| Material | Lamé parameters lambda=1 and mu=0.3, matching the manuscript's stated convention. |
| Barrier gamma | 0.5. |
| Null-space penalty | `M = 1e7 * lambda_max_raw(K)`, recomputed for the current geometry. This approved choice differs from the manuscript's initial-stiffness maximum-eigenvalue prescription. |
| Eigen tolerance | 1e-10. |
| Inner problem | Compute three smallest generalized eigenpairs, mass-normalize each pressure, solve its projected displacement, select the greatest compliance. The manuscript describes a single smallest-eigenpair solution. |
| Eigen implementation | Spectra for the sparse problem, with Eigen's dense solver on small systems; the manuscript names ARPACK. |
| E_ss domain | **Only tets having at least one boundary face.** A tet touching only a boundary vertex or edge is not included. Scalar bound minimization, feasibility tests, energy, and direct stress derivatives use the same set. |
| Maximum stress | Maximum absolute eigenvalue of the symmetric stress tensor, including compressive stress. Log boundary and whole-volume maxima separately. |
| Gradient masking | Remove the old 10%-of-maximum filter. Only interior vertex DOFs, which are not optimization variables, are zeroed. |
| Phase II bound | Hold sigma_bar=sigma_star fixed and include E_ss in Eq. (8). |
| Phase II updates | Joint weighted sum of E_ss, E_si, E_ef, E_hd, E_lds; all weights default to one. |
| Phase I | E_ss + E_si + E_ef with unit weights. |
| Local smoothing | Direct squared signed dihedral angle with edge-midpoint activation and zero rest angle. Include angle, geometric weight, and activation derivatives. |
| d0, d1, d2 | All 1e-3. d1 thresholds **tet volume**, while d0/d2 threshold distances. |
| Mesher tolerance | epsilon_d=1e-4 after longest-side normalization. |
| Iteration limits | 1000 outer iterations and 5 inner attempts; no hidden loop resets. |
| Eta relaxation | Start at .9. Relax to (eta+1)/2 after two failed Phase I steps or 15 Phase I attempts without reaching Phase II; reset those counters on Phase II. Stop once eta exceeds .99. |
| Linear solves | Projected AMG-preconditioned CG and restarted GCR. `linear_solver` chooses one; default GCR. |
| Mesh improvement | Retain MMG3D, isolated behind an in-memory adapter. Lock the boundary and reject unsafe/non-improving candidates. This differs from the mesh-improvement method cited in the paper. |
| Quality | Implement the printed formula `144*sqrt(3)*V/(sum_{i<j} ||xi-xj||^2)^(3/2)` directly; remeshing trigger .1. A rejected remesh does not enforce a final minimum of .1. Under this formula a regular tet has quality 2, not 1. |

## Exact conventions

The scalar stress energy is

`sigma_bar - gamma * sum_{t in boundary_tets} log det(sigma_bar^2 I - sigma_t^2)`.

The implementation uses the three stress eigenvalues to evaluate the log
determinant and a safeguarded scalar Newton iteration in Phase I. Phase II
returns an infinite objective outside the strict barrier domain. The auxiliary
variable and measured principal stress are different quantities.

Although E_ss is restricted to boundary tets, equilibrium and the derivative
of stiffness use **all** tetrahedra. The stress adjoint also includes derivatives
of the geometry-dependent rigid-mode projection and the unit-normal load map.
The pressure eigenvector is held fixed during the outer update, as in the
approximate bilevel method; its eigenvector derivative is not included.

The load mass is barycentric lumped surface area: each incident triangle
contributes one third of its area. Summed masses therefore equal total surface
area. The original mean-of-adjacent-areas convention does not satisfy this
identity and has been corrected.

E_hd uses area-weighted quadrature:

`sum_t (A_t/n_t) sum_i ||x_ti - closest_ti||^2 / sum_t A_t`.

Sample counts are `max(1, ceil(rho*A_t))` at topology-epoch creation. Uniform
barycentric samples are generated once for that topology with one seeded RNG.
Correspondences update before an inner attempt; sample positions and
correspondences remain fixed during its line search. Derivatives include both
the triangle-area weights and total-area normalization. A remesh creates a new
sampling epoch. This convention avoids discontinuities from resampling or
changing integer sample counts inside a line search.

E_lds uses `3*l_e^2/(A_left+A_right)` times squared dihedral angle and
the smooth activation at the edge midpoint's distance to its fixed target
correspondence. IPC uses a fixed rest geometry for edge-edge mollifiers over
each topology epoch, and its returned gradient is scattered into volume DOFs.

AMG is built on a slightly regularized stiffness matrix for the
**preconditioner only**. CG/GCR solve the unregularized stiffness equation on
the complement of all component rigid modes. Multiple disjoint components
receive six rigid modes each. A failed iterative solve rejects the trial or
returns an explicit analysis error; it never silently substitutes a solution.

The generalized operator uses the raw translation/rotation basis in its penalty,
consistent with the retained formula. Projected linear solves orthonormalize
that basis. Residual diagnostics are normwise backward errors for the full
penalized eigenproblem; because its penalty is large, these are not estimates
of relative error in its smallest eigenvalues. A dense oracle test separately
checks those eigenvalues.

## Correctness audit: disposition of all 21 findings

| Review item | Resolution |
|---|---|
| 1-2 duplicated/misweighted gradients | One composite evaluator builds both value and gradient from identical weights. |
| 3 missing collision gradient | Return IPC potential gradient; test active contact and CCD crossing. |
| 4 inconsistent shell formulations | Remove LibShell; one dihedral evaluator serves value and gradient. |
| 5 debug I/O changes energy | No numerical file output; callback-based logging in the CLI. |
| 6 pre-existing directories | CLI creates the output tree. |
| 7 target overwritten | Keep the target in memory and checkpoint it explicitly. |
| 8 missing normalization | Normalize every input; store inverse transform; restore mesh output coordinates. |
| 9 wrong sample counts/normalization | Ceil of density times area, area-weighted integration, total-area normalization. |
| 10 unseeded/random line search | One saved RNG and fixed samples/correspondences during a search. |
| 11 shared stress race | Disjoint per-element assembly storage; deterministic serial sparse accumulation. |
| 12 break inside parallel region | Remove that loop; no unsupported OpenMP control flow. |
| 13 missing return | Remove the filter and its undefined-return function. |
| 14 abort on barrier failure | Feasibility checks and exception/error propagation. |
| 15 leaked materials | Remove shell material objects; owned numerical objects use RAII. |
| 16 incomplete mesh copy | Standard value semantics copy all mesh state. |
| 17 unordered surface IDs/leaks | Sorted indices and owned vectors; no raw allocation. |
| 18 mutating quality getter | Pure quality function on the project's coordinates. |
| 19 exits/temp files in remesher | RAII MMG adapter, no temporary files, rollback on failure. |
| 20 loop mutation for remeshing | In-memory remeshing at outer boundaries. |
| 21 byte-hash cache | Explicit geometry/topology versions and fresh immutable FEM snapshots; no coordinate hash cache. |

The source is a consolidated modular rewrite, not a sequence of isolated,
behavior-preserving commits for each audit item. No per-fix benchmark table
or historical-equivalence claim is made. Validation records the frozen old
snapshot run, mathematical unit tests, and corrected runs separately.

The supplied review's 21 old-code findings have implementation dispositions
above. This is not evidence that every recommendation in that PDF was adopted:
the optional surface build remains default-OFF, and MMG remains the mesh
improver. Nor does passing the regression suite establish the absence of all
bugs or complete reproduction of the manuscript's experiments.

## Dependencies and deployment differences

The core removes CGAL, CppAD, LibShell, copied AMGCL solver headers, and unrelated
graphics/optimization libraries. Optional fTetWild still requires GMP and uses
upstream Geogram/oneTBB state. Its adapter is intended for serial initialization,
not concurrent independent mesher instances in the same process.
Input-face sorting is disabled because the pinned upstream implementation
uses an unseeded hardware-random shuffle for that option. Geogram and oneTBB
thread limits follow the configuration. The CLI saves the generated volume
mesh, and serial repeatability is tested in two fresh processes.
Geogram 1.9.6 has a second hardware-random shuffle in `basic/algorithm.h`;
`cmake/PatchGeogram.cmake` replaces that helper's seed with 42 while retaining
its BSD header. The patch checks the pinned source before applying and is
idempotent. Meshing uses this fixed dependency seed; `seed` controls the
optimizer's eigenvector initialization and quadrature sequence.

`WCSO_WITH_TETWILD` defaults to **OFF** to keep the standard build independent
of an externally installed GMP. Enabling it supports surface inputs through
the library API. This is a packaging choice different from the review's
suggested default-ON setting.

MMG candidates must preserve the exact boundary triangles and coordinates;
this gives zero surface displacement on accepted remeshes. Required boundary
vertices are matched bijectively within 32 machine epsilons times the original
coordinate magnitude or bounding-box extent, then restored to their exact
original coordinates before volume, self-intersection, and quality checks.
This handles MMG's coordinate scaling round trip without relaxing the geometry
tolerance. Ambiguous mappings, larger movement, and changed boundary triangles
are rejected. The auxiliary
`surface_deviation` routine samples vertices, edge midpoints, and face centers;
it is a diagnostic, **not** a certified continuous Hausdorff bound.

The original chair asset is excluded from the public source tree because its
redistribution terms were not supplied. Local validation used the user's copy.
No unpublished manuscript PDF, review document, machine-specific build path,
old executable, or generated optimization result is bundled as source.
