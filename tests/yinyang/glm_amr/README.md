# Yin–Yang two-level GLM MHD regression

The four-rank case starts with two equal level-one patches. Yin refines one
interior root into eight children at step 1, retains them at step 2 and
coarsens them at step 3. Yang stays at level one throughout. The chosen
interior root lies outside the angular receiver and donor stencils. The
test checks that these stencils stay at level one while actual Yin donors
move between ranks during refinement and coarsening. It also checks that
Yang donors keep their owners.

The user module reseeds the test state after every advance so the expected
donor primitives can be evaluated analytically, including after migration and
restart. The scalars and global Cartesian velocity and magnetic field are
bilinear in the donor's native angles. The reference evaluates those functions
at the transformed receiver coordinate and verifies their rotation into the
receiver basis, with separate radial boundary clamping. The case checks all
angular overset ghosts, including radial and angular corners, to an absolute
tolerance of `2e-13`. The prescribed-state comparison isolates AMR and
overset exchange behavior.

The suite also runs an initially mixed grid, restarts from the mixed step-1
checkpoint, checks fixed-level GLM and native spherical AMR, and verifies
rejection of unsupported CT AMR, a level above the native maximum and forced
I/O refinement. The companion `tests/yinyang/glm_amr_cross` suite exercises
refinement on the angular interface with the same AMR implementation.

Make a source copy in an external task directory and copy this case into a
sibling directory. Set `AMRVAC_DIR` to the source copy and build the case using the normal setup
entry point. Use `-arch=debug` for the checked build or `-arch=default` for
the optimized build.

```sh
"$AMRVAC_DIR/setup.pl" -d=3 -arch=debug
make -j1
python3 "$AMRVAC_DIR/tests/yinyang/glm_amr/run.py" ./amrvac ./regression-output
```

The output directory must be new and outside the source copy. On macOS add
`--memory-guard "$AMRVAC_DIR/tests/yinyang/interface_transport/convergence/memory_guard.py"`
to apply the project's memory thresholds to each MPI run. Use `--mpiexec` to
supply launcher options if required. Each case retains stdout and stderr,
and the successful suite writes `results.json`.
