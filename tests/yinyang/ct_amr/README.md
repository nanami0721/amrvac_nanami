# Yin–Yang CT AMR regression cases

Build with an external copy of AMRVAC and point `AMRVAC_DIR` at that
copy. Use `setup.pl -d=3 -arch=debug` or `-arch=default` in a separate
external case directory. Copy `mod_usr.t` and the selected parameter
file there before building. Keep each MPI run and its stderr in a new
output directory.

- `amrvac.par` refines an interior level-one block after iteration 1,
  evolves on the finer tree, and requests coarsening at iteration 4.
  Run with two MPI ranks. Repeat with `mhd_internal_e=T` in `&mhd_list`
  for the alternative energy formulation.
- `dynamic_wall.par` refines a block touching the inner wall and theta
  seam, allows level 4, and requests coarsening from iteration 4. The
  stage test measures any change in the nonzero prescribed tangential
  wall EMFs and checks that the opposite patch's plan uses two fine
  integrals with the native coarse measure where a coarse edge shares
  their line. It reports `CT AMR wall pairs`. A second check applies the
  donor-aligned prolongation to eight virtual children of actual coarse
  boundary cells and compares each child's oriented face-flux sum with
  one eighth of its parent flux sum. The fixture also checks active
  divergence and thermal pressure. Use two MPI ranks, then continue the paired
  iteration-3 snapshot with four ranks using the literal filename
  placeholder `-if 'ct-amr_{patch}0003.dat'`.
- `run_transport.py EXECUTABLE EXTERNAL_DIRECTORY` takes the externally
  built `tests/yinyang/ct_transport/mod_usr.t` executable. It runs both
  directions for 60 steps, comparing an all-L2 layout with the equivalent
  uniform mesh and comparing 2 and 4 ranks. It also evolves L1 on Yin and
  L2 on Yang. The runner preserves the MPI logs and writes the matched
  active cell and face comparisons to `comparison.json`.

The `CT AMR` lines report iteration, number of staged calls, leaf counts
in each patch by level, maximum active and full-block divergence, and
minimum thermal pressure. The checks use the code's staggered divergence.
