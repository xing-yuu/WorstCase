# Validation record

Verified locally on 2026-10-21. The Windows machine uses an Intel Core
i7-13700KF. Builds use Release optimization and the revisions pinned in the
project. Local dependency caches/source overrides were used; this is not a
claim of an independently provisioned fresh CI runner.

This source version is ready for review and testing. Three selected full
main-method runs have completed: chair and plane on the earlier frozen Lamé
build, and leaf on the corrected MMG-boundary build. This is not a complete
three-model validation of one final executable. The MMG boundary correction
and higher-precision solver reporting passed the Windows core/optional and
WSL regression suites below. The limited beam/chair measurements remain
observations of the earlier Lamé build; they are distinct from the full runs.
Known manuscript and review-recommendation differences are documented in
`implementation_notes.md`; completion does not establish a bug-free program.

| Environment | Configuration | Result |
|---|---|---|
| Windows, MSVC 19.29 / VS2019, CMake 3.27 | Core, OpenMP ON | 3/3 CTest tests passed |
| Windows, same compiler | fTetWild ON, OpenMP ON | 4/4 CTest tests passed |
| Ubuntu in WSL2, GCC 11.4, CMake 3.27 | Core, OpenMP OFF | 3/3 CTest tests passed |
| GitHub Actions | Ubuntu/macOS/Windows core, Ubuntu surface job | Workflow supplied; not executed here |

The WSL build reused the pinned sources downloaded on Windows. MMG's shell
script was converted from the Windows checkout's CRLF to LF for this local
cross-platform test; archives fetched normally on Linux already use LF.
Our own compilation units have no reported compiler warnings in the final
build logs. Upstream libraries still emit warnings under some compilers.
The checked-in clang-format 23.1.3 style and source-hygiene check pass.

## Numerical tests

All FEM and optimization validation below uses Lamé parameters lambda=1 and
mu=0.3, matching the manuscript convention.

- FEM symmetry, exactly six rigid modes, rigid-mode annihilation,
  constant-strain/uniaxial stress, and analytic shape-function derivatives.
- Barycentric mass conservation and manufactured projected linear solves.
  CG residual: 1.37e-15; GCR residual: 9.30e-16 on the small fixture.
- Six directions with three finite-difference step sizes for each energy.
  The largest best-step relative error per energy on Windows was:

| Energy | Observed directional error |
|---|---:|
| E_ss, optimized scalar bound | 1.72e-5 |
| E_ss, fixed scalar bound | 2.31e-8 |
| E_ef | 1.04e-9 |
| E_hd | 8.60e-11 |
| E_lds | 1.48e-8 |
| E_si, active near-contact fixture | 9.69e-8 |

- E_ss boundary-only sum, principal-stress rotation invariance, active
  collision gradient, and CCD rejection of a crossing.
- Three eigenpairs, mass normalization, maximum-compliance selection, and
  comparison with a dense generalized eigensolver on a thin beam.
- MMG improves a rotated, translated interior-distorted fixture, restores its
  exact boundary coordinates and triangles, and rolls back a failed call.
  Boundary regression cases cover one-ULP changes with vertex renumbering,
  real displacement, nonbijective/ambiguous matches, changed triangle
  topology, atomic rejection, and consistent accepted reference coordinates.

## CLI and regression tests

The procedural beam has 80 vertices, 216 tetrahedra, and 112 boundary
tetrahedra. Both phases are exercised. Tests compare CG/GCR and thread settings,
reject invalid configurations and malformed checkpoints, confirm no input
modification, and check identical final OBJ bytes for continuous and resumed
runs. ASCII Gmsh 2.2/4.1 node tags, translated/scaled units, annotated Medit
input, and output overwrite protection are covered.

`tests/reference/beam.json` stores the initial analysis, scheduler events,
per-step energy history, and final summary. Cross-platform tolerances are
2e-4 for initial stress analysis and 5e-3 for the nonlinear trajectory; schedule
events must match exactly. Timing is excluded from regression comparison.

The Windows optional build tests OBJ, PLY, and STL through actual
tetrahedralization and stress analysis. It also compares two fresh serial
processes: their generated volume meshes are byte-identical after the
documented Geogram shuffle patch.

The benchmark and plotting scripts were exercised on the included beam.
Its smoke run took 0.151 seconds wall time and
16.7 MiB peak working set in a sampled Windows measurement.
The final critical-load boundary stress was
496.456752, from
889.583409.
These are coarse-fixture results, not estimates for production models.

## Supplied chair validation

Input: the supplied `067_rocking_chair_mmg.mesh`, with 28155 tets
and 9160 boundary-face tets. Configuration: corrected
paper defaults, seed 20260926, one OpenMP thread, two outer iterations and five
inner attempts per outer iteration. All ten attempts completed.
The final result was reanalyzed with a newly computed critical load.

| Metric | Initial | Final |
|---|---:|---:|
| Maximum absolute principal stress, boundary tets | 116340.887112 | 84608.030575 |
| Maximum absolute principal stress, all tets | 116340.887112 | 84608.030575 |
| E_hd | 0 | 3.14019499e-09 |

Boundary stress reduction: **27.276%**.
Whole-volume stress reduction: **27.276%**.
Recorded solver time, including final analysis: 395.206 seconds.
The initial analysis uses the same original mesh, normalization, sampling/RNG
initialization, and configuration as the run; its scalar bound matches the
first recorded sigma_star/eta within 1e-7.

The frozen old snapshot also completed its configured Nu=2, N=5 run in
932.470 seconds. Its old stress measure is an axis-dependent diagonal
maximum and its defaults/loop semantics differ. It is not used to compute a
stress-improvement claim for the rewrite. Timings were collected while other
verification work was active and are not controlled performance comparisons.
The original ZIP SHA-256 is recorded in `validation_results.json`.

The chair mesh and generated chair results are excluded from this source
package because redistribution terms for the supplied asset were not provided.
The scalar validation record is included.

## Full main-method runs and version attribution

All selected runs used Lamé lambda=1, mu=0.3, the same numerical configuration,
and the supplied original inputs. They exited with code 0, completed the
sequential outer counts below, and stopped at `eta_limit` with eta=0.99375.
The configured 1000 outer iterations are a maximum, not a completion target.
Final stress means maximum absolute principal stress from `final_analysis`,
after recomputing the critical load on the final shape. Values are in
normalized analysis units; boundary and whole-volume scopes are distinct.

| Selected run | Outer iterations | Final boundary stress | Boundary reduction | Final all-tet stress | All-tet reduction | Minimum quality |
|---|---:|---:|---:|---:|---:|---:|
| Rocking chair, earlier frozen build | 121 | 20782.173145 | 82.136827% | 24095.172626 | 79.289162% | 0.4772693331 |
| Plane, earlier frozen build | 43 | 1096367.740959 | 73.573869% | 3337942.261802 | 34.451895% | 0.0839047297 |
| Leaf, corrected build, fresh full run | 31 | 140264.143544 | 97.513878% | 140264.143544 | 97.513878% | 0.1865300515 |

Chair and plane are selected from `outputs/paper-experiments-lame-rerun`.
Their executable SHA-256 is
`b38f9885364b100c412eb23c49f35907ee070c07df9cfda9bb819566735dbf2c`.
The corrected leaf is selected from
`outputs/paper-experiments-lame-leaf-remesh-fix`; its executable SHA-256 is
`a2f89e1b3e782d70ade8f24e3a18440a4d3a8d34ca1d9a77e1dc57cced093d93`.
The corrected leaf's frozen source ZIP SHA-256 is
`4fc59c8aaa637066f6fc9ea6e8c18fbc6ef196d6a4eeb853e59b9836684f2154`;
that immutable snapshot predates these documentation updates.

All selected final tetrahedra have positive volume. Plane has one tetrahedron
below the 0.1 remeshing trigger. Its only remesh attempt, at outer 43, was
rolled back with a generic boundary-preservation or quality-validation message.
The earlier binary does not identify the precise rejection cause; it is not
established that this plane candidate was rejected for the leaf roundoff bug.
Positive volumes do not imply that the quality threshold was met. Corrected
leaf accepted one remesh at outer 7, with no remesh rollback.

Measured process wall times were 19314.335876 s (chair), 7520.026748 s (plane),
and 2322.706766 s (corrected leaf), on an i7-13700KF with 32 GB RAM. The original
three tasks began together on 2026-10-21; the failed leaf was replaced later
by a fresh run. Diagnostic builds and tests also competed for resources.
These are observations of overlapping runs and two program versions, not
a controlled same-version three-task parallel-performance benchmark.
The combined Chinese report and detailed audit are stored locally outside the
source package in `outputs/paper-experiments-lame-final/REPORT_ZH.md`.

## Preserved leaf failure, repair, and recovery

The original leaf in `outputs/paper-experiments-lame-rerun` exited with code 1
during outer iteration 15 after inner step 2, at 2026-10-21 16:30:48 +08:00.
Its last complete checkpoint is outer 14. Eight MMG candidates were rejected
in outer iterations 7-14. This failed run has no final analysis or final stress
result and is excluded from the selected completed runs.

This failure belongs to the frozen build before the boundary correction,
preserved as `provenance/source.zip` in that experiment directory (SHA-256
`a132cd3b10c3c613e3b9fea8ea68f98b489ce108ec74009861af8db0247cb1b5`).
Its executable SHA-256 is
`b38f9885364b100c412eb23c49f35907ee070c07df9cfda9bb819566735dbf2c`.
The original evidence and binary are unchanged.

Diagnosis using the actual outer-14 checkpoint and frozen MMG DLL confirmed
that exact boundary-coordinate comparison rejects floating-point round-trip
changes. Both meshes had 2830 boundary vertices with a bijective matching and
identical boundary triangle topology. There were 1454 changed coordinate
entries; maximum boundary displacement was 1.2490e-16. Restoring the original
boundary coordinates exactly gave minimum quality 0.23944070154826996, from
0.042314374738294976, with all orientations positive. The original comparison
rejects before intersection validation.

The correction is implemented: a bijective correspondence within a
machine-roundoff bound and identical boundary triangle topology are required
before restoring the original boundary coordinates exactly. The bound is
independent of the configured geometry tolerance. Positive-volume,
self-intersection and strict quality-improvement checks remain required.
The corrected production C++ path accepted the actual outer-14 checkpoint
with exit code 0 and `remesh_accepted`; its IPC intersection check executed.
Exact boundary coordinates/triangles and reference-coordinate consistency
passed. Minimum quality improved from 0.04231437473829487 to
0.2394407015482699, with minimum volume 5.8350893766187475e-9 and all volumes
positive. Vertices changed from 3665 to 3596 and tetrahedra from 12349 to
11953. Local evidence is
`work/leaf-diagnostics/mmg_boundary_leaf_outer14_cpp.json`. The new boundary
regressions passed in Windows core (3/3), Windows optional (4/4) and WSL core
(3/3) suites. These checkpoint checks validate remeshing; the subsequent full
optimization recovery has separate evidence below.

A diagnostic reproduction located the terminating solve in the stress
adjoint: GCR used 3000 of 3000 iterations. The recursive relative residual
was 1.49514260022547979e-6 and the final true projected relative residual was
1.49514615118338084e-6, against tolerance 1e-8 (RHS norm
0.0883146662722424). This was genuine iteration-budget exhaustion, with no
early false convergence observed. The reproduction preserves the original
configuration and matches the last steps apart from timing; its local record
is `work/leaf-diagnostics/linear_failure.json`. Higher-precision error
reporting is implemented and the regression suites passed; the numerical
solver algorithm, tolerance and iteration budget are unchanged.

After the regression and actual-checkpoint checks, corrected leaf started
fresh from the original input at 2026-10-21 16:59:46 +08:00 and finished at
17:38:29 +08:00. Its numerical configuration is identical to the original
formal run. It completed 31 outer iterations, including terminal critical-load
analysis, and its selected results are reported above. The degraded old
checkpoint was not resumed. Original failure evidence remains in
`LEAF_FAILURE.json`, `LEAF_FAILURE.md`, and the diagnosis files in the failed
experiment directory, without replacing that original run's data.

The earlier E=1, nu=0.3 results are historical and do not validate this material
version. A prior Lamé attempt was interrupted before completion, with no
observed solver exit code; it is not counted as a completed experiment.

## Remaining empirical limits

Chair and plane have not been rerun with the corrected leaf executable.
Thus the three selected completed results do not fully validate the current
code version on all three models. No other paper baseline, ablation,
quadratic-FEM comparison, or physical experiment is claimed. There is no
per-bug isolated-commit benchmark table or
behavior-preserving equivalence claim. macOS and the optional Linux surface
build have CI definitions but were not run locally; hosted GitHub Actions was
not executed here. The implementation retains documented differences from
the manuscript: current-geometry 1e7-scaled penalty, three-candidate load
selection, Spectra instead of ARPACK, and MMG instead of the cited mesh
improver. These author-approved decisions must not be presented as literal
agreement with every manuscript statement or every review recommendation.
