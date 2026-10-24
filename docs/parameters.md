# Parameters and units

All fields can appear in a JSON object. Omitted fields use the defaults below.
`Config::from_json` rejects unknown keys, wrong types, nonfinite values, and
invalid ranges. Command-line values override JSON values.

The mesh is translated by its bounding-box minimum and divided by the longest
bounding-box side. Length and volume thresholds apply after this transform.
Reported stress remains in normalized analysis units. Normalization metadata
is written with every run.

| JSON key | Symbol | Default | Meaning / allowed values |
|---|---|---|---|
| `lame_lambda` | lambda | `1.0` | First Lamé parameter in normalized analysis units. |
| `lame_mu` | mu | `0.3` | Second Lamé parameter (shear modulus); positive. |
| `gamma` | gamma | `0.5` | Positive stress log-barrier coefficient. |
| `nullspace_multiplier` | M / lambda_max_raw(K) | `10000000.0` | Positive multiplier of the largest stiffness eigenvalue. |
| `eigen_tolerance` | epsilon_p | `1e-10` | Generalized eigensolver stopping tolerance, in (0,1). |
| `linear_tolerance` | - | `1e-08` | Relative projected linear residual tolerance, in (0,1). |
| `eigen_max_iterations` | - | `10000` | Maximum eigen iterations; positive integer. |
| `linear_max_iterations` | - | `3000` | Maximum CG/GCR iterations; positive integer. |
| `gcr_restart` | - | `40` | Positive GCR basis length before restart. |
| `linear_solver` | - | `"gcr"` | `gcr` or `cg`; both projected and AMG-preconditioned. |
| `w_ss` | - | `1.0` | Positive stress-energy weight. |
| `w_si` | - | `1.0` | Positive contact-energy weight. |
| `w_ef` | - | `1.0` | Positive element-flip barrier weight. |
| `w_hd` | w | `1.0` | Nonnegative shape-fidelity weight, active in Phase II. |
| `w_lds` | w | `1.0` | Nonnegative local-smoothing weight, active in Phase II. |
| `d0` | d_0 | `0.001` | Positive contact activation distance; normalized length. |
| `d1` | d_1 | `0.001` | Positive element-flip activation **volume**; normalized length cubed. |
| `d2` | d_2 | `0.001` | Positive smoothing activation distance; normalized length. |
| `epsilon_d` | epsilon_d | `0.0001` | Positive mesher envelope/hausd parameter; normalized length. |
| `quality_threshold` | Q_min | `0.1` | Trigger remeshing below this quality; [0,2.01). |
| `tet_edge_length` | - | `0.025` | Positive target mesher edge length; normalized length. |
| `eta_initial` | eta | `0.9` | Initial target-bound multiplier; 0<eta_initial<eta_stop<1. |
| `eta_stop` | - | `0.99` | Terminate when relaxed eta exceeds this value. |
| `outer_iterations` | N_u | `1000` | Positive maximum number of completed outer iterations. |
| `inner_iterations` | N | `5` | Positive maximum inner attempts per outer iteration. |
| `max_phase1_steps` | - | `15` | Relax eta after this many Phase I attempts without reaching Phase II. |
| `max_failed_steps` | - | `2` | Relax eta after this many consecutive failed Phase I attempts. |
| `max_backtracks` | - | `50` | Positive maximum halvings per line search. |
| `armijo` | - | `0.0001` | Sufficient-decrease constant in (0,1). |
| `step_fraction` | - | `0.25` | Maximum initial displacement relative to shortest edge; (0,0.5]. |
| `step_tolerance` | - | `1e-12` | Positive minimum proposed maximum displacement; normalized length. |
| `sampling_density` | rho | `10000.0` | Samples per normalized area; count=max(1,ceil(rho*A)); (0,1e8]. |
| `threads` | - | `1` | 1..256 OpenMP threads if available. Mesher internals may use their own runtime. |
| `save_every` | - | `0` | Surface snapshot interval in outer iterations; 0 disables. |
| `checkpoint_every` | - | `1` | Checkpoint interval in outer iterations; 0 disables. |
| `seed` | - | `1` | Unsigned 64-bit seed for pressure initialization and quadrature sampling. |
| `remesh` | - | `true` | Enable in-memory MMG quality improvement at outer boundaries. |
| `debug` | - | `false` | Log every attempt, in addition to outer/final progress. |

`configs/paper.json` gives the author-approved corrected settings. The penalty,
candidate eigenpair policy, and implementation dependencies deliberately differ
from the manuscript; see `implementation_notes.md`.
`configs/smoke.json` deliberately reduces iteration count and quadrature density.
Every run saves all effective fields, including defaults, in `config_used.json`.
No legacy environment variables configure numerical behavior.

The E_ss domain is fixed to boundary-face tetrahedra by design and is not a
configuration option. The inner eigenpair count is fixed at three.
