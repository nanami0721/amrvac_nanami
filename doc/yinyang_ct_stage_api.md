# Experimental staged CT API

This is the RK-stage foundation for overset electric-field exchange. Yin–Yang
now has a separate, explicit `yy_ct_experimental` integration path described below.
The default remains disabled; background-field splitting remains unsupported.

## Configuration and lifetime

`ct_stage_split=.true.` in `&methodlist` selects a two-pass update. The default is
false, preserving the original block-at-a-time route. The experimental route
requires 3D, global unsplit SSPRK3 stepping, TVDLF or HLL fluxes and UCT-HLL.
The standalone staged route uses a fixed level-one mesh. The explicit
Yin–Yang CT route also accepts native AMR. Cost-weighted balancing and local
timesteps are outside this configuration.

Each active local block owns cached numerical face fluxes and time-integrated
edge electromotive forces. Arrays are allocated outside threaded kernels and
reused across stages. Workspace is not snapshot state and is regenerated after
restart. It adds storage proportional to local block count; the memory watchdog
must cover the complete MPI process tree.

## Sequence

1. Convert the RK input state, reconstruct faces and compute Riemann fluxes and
   UCT-HLL edge line integrals for every local active block. Do not update the
   destination state or invoke geometric/physical source updates in this pass.
2. Call `usr_ct_stage(qt,qdt,sCT,s,fluxes)` once per rank, outside OpenMP loops.
   Different ranks may own different numbers of blocks. Yin–Yang edge exchange
   follows this callback on the world communicator.
3. Apply the discrete curl to face fields, then the original finite-volume update,
   geometric source, face-to-center synchronization, small-value handling and
   physical sources in their original order. Store final edge integrals for the
   existing same-patch CT correction path.
4. Apply existing flux/edge corrections and fill ghost cells.

`phys_get_ct_edges` computes edge integrals with the existing face-update signature;
`phys_apply_ct_edges(ixI,ixO,fE,s)` applies their discrete curl. UCT-HLL registers
both. The original `phys_update_faces` entry remains as a wrapper calling both,
so original call sites and other CT choices keep their existing behavior.

`ct_stage_flux` holds `face(spatial,nwflux,ndim)` and `edge(spatial,sdim:3)` arrays.
Only active local blocks are populated. Entries outside the generated stencil
are initialized to zero; zero initialization does not make them valid donors.
The exchange callback may change valid edge line integrals. States and face
fluxes are read-only by contract. The integrals already include substage dt and
physical edge length: do not multiply by those factors again. The old per-block
`usr_set_electric_field` hook runs during edge generation, before this callback.

Register `usr_ct_stage` unconditionally in `usr_init`; `&methodlist` has not yet
been read at that point. Check configuration-dependent requirements in
`usr_set_parameters`, after the main parameter reader. The callback is invoked
only when staging is enabled.

The existing compute-only load-balance timer excludes the new callback duration.
No speedup is claimed: staging trades persistent workspace and a second primitive
conversion for an all-block communication boundary.

## Verification definitions

`tests/yinyang/ct_preflight/stage_suite.py EXECUTABLE EXTERNAL_DIRECTORY` runs under
its process-tree watchdog. It checks native versus staged snapshots, including
staggered data, by complete-file SHA-256 equality at 1/2/3 ranks. Twenty blocks
produce an uneven partition at three ranks. Cases include uniform static fields
and nonuniform discrete-curl fields with finite velocity. A half-strength EMF
applied through the stage hook must match the old per-block hook and must differ
from the unmodified solution. Contact CT with staging must be rejected.

The test does not establish Yin–Yang interface consistency, ghost-face divergence
constraints, overlap convergence or global magnetic-flux conservation.

## Validation evidence

The 2026-11-17 temporary test directories are no longer available. Do not use
those paths as retained evidence. Fresh debug/optimized suites and dimension
builds passed and are recorded in the persistent external directory
`/Users/nanami/workspace/2608_uawsom/yinyang-ct-stage-20261118`.

The reproducible entry point is
`tests/yinyang/ct_preflight/build_stage_suite.py EXTERNAL_DIRECTORY`, with optional
`--cache-source EXTERNAL_SOURCE` and `--resume`. Builds use a copied source tree,
serial compilation and per-job process-tree memory guards. Resume requires exact
source/test manifest agreement. `status.json` records completed tasks, while
`manifest.json`, `source.diff` and `commit.txt` identify the source snapshot.
Full accepted results are listed in `yinyang_fixed_grid_progress.md`.

On 2026-11-18 both suites passed 19 runs and 13 binary comparisons each; all six
full debug dimension libraries compiled. Three-rank staged CT restart was binary
identical to uninterrupted evolution. Existing Yin–Yang GLM matched its retained
baseline at 2/4 ranks. See the progress document and external `final-audit.json`
for divergence and memory measurements. Actual Yin–Yang CT remains disabled by default.

## Explicit Yin–Yang transfer components (2026-11-18)

`yy_prepare_ct_edges` and `yy_exchange_ct_edges` now provide a caller-owned plan
and an adapter for `ct_stage_flux%edge`. The latter is intended for this stage
callback and packs all donors before overwriting receiver edges. It does not
register itself. `yy_ct_close_ghosts` provides angular ghost constraints after
magnetic ghost faces have been populated. See
`tests/yinyang/ct_coupling/README.md` and the progress document for tested contracts
and the subsequent full-MHD integration described below. CT is disabled by default.

## Magnetic face population (2026-11-18)

`yy_prepare_ct_ghosts`/`yy_fill_ct_ghosts` now supply the face population step before
angular divergence closure. They protect active/boundary faces, use only physical
face donors and support radial walls via scalar Br interpolation. They are
explicit APIs and do not register production ghost hooks. Refer to
`tests/yinyang/ct_ghost_fill/README.md` for coverage, accuracy and remaining
integration requirements.

## Same-patch ordering (2026-11-18)

`yy_prepare_ct_patch` caches geometrically owned shared-edge copies and two
magnetic face halo schedules. `yy_ordered_ct_edges` synchronizes duplicate edges,
exchanges owner-only overset receivers, then synchronizes the updated duplicates.
`yy_ordered_ct_ghosts` populates physical-domain same-patch halos, fills and
constrains overset angular halos, then completes mixed block/boundary corners.
All ranks call these wrappers outside block/OpenMP loops. Active magnetic faces
and physical radial ghost cells are protected. See
`tests/yinyang/ct_patch_sync/README.md` for the interface and probe tests.

These are explicit integration components. At that component milestone, production hook registration,
radial-wall EMF closure, physical radial boundary sequencing, and cell-centered
B/energy consistency still precede full evolving Yin–Yang CT activation.

## Radial boundary and variable completion (2026-11-18)

The opt-in `radial_extrapolate=.true.` argument to `yy_prepare_ct_edges` closes
E_r sampling at physical radial endpoints using the nearest two physical radial
samples. It admits no point outside the shell and does not weaken angular donor
coverage. Existing calls retain their strict no-extrapolation behavior.

After ordinary and overset **cell** ghost exchange, explicit face/state completion
is now available in this order:

```fortran
call yy_complete_ct_faces(patch_plan,ghost_plan,s)
call yy_sync_ct_state(s,mhd_ct_sync_center)
```

The first call finishes the ordered angular face exchange, then extends radial
tangential faces and closes radial normal ghost fluxes at zero discrete divergence.
The second uses the native face-to-center operator: active total energy is fixed,
angular ghost thermal energy is retained, then consistent boundary cell data are
copied into radial ghosts and their magnetic energy is compensated. All radial
normal boundary anchors remain fixed during ghost fill. This is the continuous
radial extrapolation/zero-divergence ghost policy; no line-tied electric-field
condition or characteristic wind outflow is imposed. See
`tests/yinyang/ct_boundary_sync/README.md` for supported energy modes and tests.

The later experimental integration section describes stage/getbc registration
and static MHD evolution. See the adaptive integration section below for
paired restart and the radial-boundary scope.

When wiring the production path, retain the existing EoS/cache refresh lifecycle
after B/energy completion. The probe tests use classical total-energy GLM MHD;
LTE/radiation-specific auxiliary-cache validation is not covered by this suite.

## Actual experimental Yin–Yang CT integration

The explicit opt-in configuration is:

```fortran
&methodlist
 ct_stage_split=T
 yy_ct_experimental=T
 time_stepper='threestep'
/
&mhd_list
 typedivbfix='ct'
 type_ct='uct_hll'
 boundary_divbfix=6*.false.
/
```

The mesh must satisfy the 3D/global-step restrictions. Classical total and
internal energy are supported. The outer radial boundary uses `cont`. The
inner radial boundary accepts either `cont` or the paired fluid and edge
callbacks described in [the Yin–Yang usage guide](yinyang_v1.md#supported-scope).
The restriction on additional MHD physics applies across all mesh levels.

`advect1` invokes `yy_runtime_ct_edges` after the user stage callback and before
any face curl update. The adapter caches geometric-owner angular interface edges,
including radial endpoints, and applies the tested same-patch/overset sequence.
`getbc` supplies radial cell donors before angular primitive exchange, then calls
`yy_runtime_ct_ghosts` to fill faces and synchronize cell variables. It bypasses
the old radial physical-boundary callback only in this explicitly selected mode.
On a fixed single-level mesh, the edge plan checks whether every sampled edge
is read from its geometric owner. If so, the runtime exchanges overset edges
directly and distributes the updated owner values afterwards. The ordinary
sibling face fill in `getbc` supplies the interior ghosts before overset face
interpolation. The final extended face synchronization still fills mixed corners.
Adaptive meshes retain the preliminary edge and face synchronizations.
The existing post-getbc EoS refresh is retained. There are no active-field resets.

`tests/yinyang/ct_static` exercises the real MHD time integrator with zero and
uniform Cartesian magnetic fields. A native spherical zero-field control measures
the existing geometric-source imbalance of uniform pressure. Native divB,
magnetic drift, analytic-field error, spurious speed, positivity and callback
counts are recorded over time. This test does not qualify magnetic transport,
low-beta evolution, restart or alternative radial boundary models.


`tests/yinyang/ct_transport` adds actual weak compact-loop advection across the
seam in both directions, with interior controls, angular refinement, amplitude
halving and centered/face MPI consistency. Native divB remains near roundoff;
substantial diffusion occurs in both seam and interior cases. This establishes
stable transport, not high-accuracy, full-order, low-beta or restart qualification.
See the dated transport section in `yinyang_fixed_grid_progress.md`.


## Adaptive CT stage and ghosts

Set `refine_max_level>1` to allow native Yin–Yang refinement. The sparse
plan uses actual patch-local leaf IDs and world ranks. An all-world
invalidation follows the native regrid and ownership update, so both
patches rebuild even when only one tree changes. Restart rebuilds the plan
from both restored trees. Version-2 paired adaptive CT snapshots also
support a different valid even MPI world size.

1. Balance the proposed leaf layout using the geometry of the receiver's
   face and edge stencils. This includes angular ghost halos at the
   radial wall and at native block interfaces.
2. Synchronize same-level edge copies under the geometric owner. In the
   opposite patch, build the receiving level's virtual donor-lattice
   samples. Equal-level requests take native data. Coarse aligned edge
   requests sum available fine segments where they own the physical
   edge. Fine-to-coarse face transfer adds four magnetic subface fluxes.
3. Prolong coarse data with tensor products of one-dimensional nodal
   and integral weights. The corresponding nodal difference and
   integrated maps commute on a complete donor-aligned stencil. The
   wider pair is used on a uniform coarse donor patch. A mixed donor patch
   uses the same linear/constant pair on every face and edge, so crossing
   a level transition does not change the transfer within a child cell.
   A face conversion uses the native area.
   An edge conversion uses the native length and retains the substage
   time factor supplied by UCT.
4. Reconstruct at the rotated receiver and project onto its normal or
   tangent. All active levels pack their native edge data before any
   overset receiving edge is written. The physical inner tangential
   prescription is excluded from receiving edge updates.
5. Apply the discrete curl with the final stage edge integrals, followed
   by native flux and edge reflux. Complete angular face halos, close
   their discrete magnetic divergence, and then fill the radial face
   halos and synchronize cell-centered fields and energy.

The donor-aligned map is distinct from rotated projection between
patches. Their composition has a measurable curl commutator. The
experimental face-moment option keeps its existing default and scope.

## Consistent native EMF representation (2026-11-19)

`yy_prepare_ct_edges` now defaults to native midpoint quadrature. Both donor
and receiver values represent `dt * edge_length * E_tangent(edge_midpoint)`.
Donors are divided only by edge length, reconstructed at the receiving edge
midpoint in Cartesian components, projected on its local tangent, and multiplied
by the native edge length. Zero dt remains safe. This is a consistent discrete
integral approximation, not an exact line integral or a higher-order UCT scheme.
Native spherical CT and its interior EMFs are unchanged.

An optional `native_midpoint=.false.` in this explicit API retains the previous
two-node receiver quadrature for diagnostic comparisons. At the lower-level
`yy_ct_edge_build` API, the historical default remains two-node quadrature;
the mesh adapter passes its selection explicitly. Selection occurs at plan
construction, with one or two samples/weights per receiver, no new per-stage
allocation and no string dispatch. Same-patch ordering and face closure are
unchanged. See `tests/yinyang/ct_accuracy` for matched legacy/current operator
tests and `ct_transport/run_representation.py` for paired evolved tests.

## Experimental face-average solenoidal reconstruction (2026-11-19)

`yy_ct_face_moments=T` in `&methodlist` opts into the new face reconstruction.
It requires `yy_ct_experimental=T`; the default is **false** because the tested
nonlinear transport cases did not improve. The explicit component API is
`yy_prepare_ct_ghosts(plan,s,face_moments=.true.)`, independently selectable in
analytic tests. Omitting the optional argument preserves the legacy point-based
face transfer. Edge midpoint representation remains enabled as previously.

The new plan fits B=b+G*(x-x_target)/h with trace(G)=0 to 24 **face averages**
(11 coefficients). Surface moments use four-point Gauss quadrature in each face
coordinate and native spherical face areas. A reorthogonalized QR solves the
small geometric fit once per receiver at initialization. Trilinear locality
weights have a fixed positive 0.01 floor before taking square roots for weighted
least squares; the floor avoids rank loss at exact lattice nodes. Radial/angular
node ties are canonicalized within geometry tolerances, preventing roundoff
from changing the positive-weight stencil. Ill-conditioned fits or failed
moment reproduction abort explicitly rather than silently changing methods.

Receiver radial/tangential plans both draw all three donor component families.
`plan%faces%radial` labels the **receiver** family; its requests no longer imply
component1. Packing must use request row5. Radial-wall fits use physical faces
only, with one-sided stencils for radially centered tangential components. No
new per-stage fit or allocation is introduced; radial sampling/communication
increases from eight to 24 samples per receiver. This is not a performance win.

Each face is fitted independently. Although each polynomial is solenoidal, the
collection of reconstructed faces is not guaranteed discretely divergence-free.
The existing outward ghost flux closure remains necessary and can amplify
nonlinear interpolation errors. This option is not a limiter, a proof of global
flux conservation, or a validated improvement for shocks/low-beta eruptions.
