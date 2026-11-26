# Yin–Yang tests

New Mac Studio / new developer: start with the repository-root
[handoff guide](../../YINYANG_MAC_STUDIO_HANDOFF.md) for the current
capabilities and a cache-free bootstrap.

See `doc/yinyang_v1.md` for scope, geometry, restart and validation details.
Never compile or run inside this source tree. `run_validation.py` copies the
source and MHD case to a new external directory, checks MPI/restart equivalence,
and optionally runs convergence/perturbation tests and all six full-library builds.

`geometry_preflight.py` is an independent geometry-only probe; it does not validate
VACPP expansion, the Fortran exchange, or MHD evolution.

`mhd/` supplies uniform, entropy-pulse and blast initial conditions. `hd/` supplies
uniform-flow/perturbation initial conditions. `rho/` is stationary scalar
preservation. Diagnostics in `check_mhd.py` require NumPy and native snapshots.

For communication-only changes, `compare_exchange.py BASELINE_RUN CANDIDATE_RUN`
compares all stored fields exactly and checks donor reduction diagnostics and
runtime stderr. See `doc/yinyang_exchange_validation.md` for the deduplication
implementation and measured costs.

`interface_transport/` adds exact forward/reverse entropy and force-free magnetic
transport, three-resolution seam convergence, interior controls and seam-inclusive
divergence diagnostics. Results: `doc/yinyang_interface_validation.md`.

The `glm_amr/` regression exercises two-level interior refinement and
coarsening, migration, paired restart and fixed-grid and native spherical
controls. The `glm_amr_cross/` regression exercises refinement on the seam,
the one-level stencil bound, mixed-level transfer and physical radial corners.
Both use the same GLM AMR path.


The [CT AMR fixture](ct_amr/README.md) covers interior refine/evolve/coarsen,
fully matched fine levels, cross-level transport, a dynamic L4 seam and a
prescribed inner-wall edge test at a coarse/fine seam junction.
