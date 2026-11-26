#!/usr/bin/env python3
"""Compare uniform and matched L2 CT meshes and run both directions on a mixed seam."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


def final_fields(directory, family, key_columns):
    ranks = directory.glob(f"{family}-rank-*.txt")
    rows = np.concatenate([np.loadtxt(path, ndmin=2) for path in ranks])
    order = np.lexsort(rows[:, :key_columns].T[::-1])
    return rows[order]


def compare(a, b, family, key_columns):
    first = final_fields(a, family, key_columns)
    second = final_fields(b, family, key_columns)
    np.testing.assert_array_equal(first[:, :key_columns], second[:, :key_columns])
    error = np.max(np.abs(first[:, key_columns:] - second[:, key_columns:]))
    if error > 5e-12:
        raise AssertionError(f"{family} differs by {error:.4g} in {a.name} / {b.name}")
    return float(error)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exe", type=Path, help="externally built tests/yinyang/ct_transport executable")
    parser.add_argument("out", type=Path, help="new external directory for the eight runs")
    parser.add_argument("--mpi", default="mpirun")
    args = parser.parse_args()
    fixture = Path(__file__).resolve().parent
    source = fixture.parents[2]
    guard = fixture.parent / "interface_transport" / "convergence" / "memory_guard.py"
    executable = args.exe.resolve(strict=True)
    out = args.out.resolve()
    if source in (out, *out.parents) or source in executable.parents:
        parser.error("build and run outside the source tree")
    out.mkdir()
    base = (fixture / "amrvac.par").read_text()
    base = base.replace("it_max=6", "it_max=60").replace("ditsave_dat=5", "ditsave_dat=60")
    base = base.replace("dtpar=0.001d0", "dtpar=0.01d0")
    base = base.replace("block_nx2=4", "block_nx2=8").replace("block_nx3=4", "block_nx3=16")
    base = base.replace("flux_scheme=20*'tvdlf'", "flux_scheme=20*'hll'")
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    comparison = {}
    for direction in (-1, 1):
        runs = []
        for label, ranks, cells, level, patches in [
            ("uniform", 2, (16, 32, 64), 1, 0),
            ("matched", 2, (8, 16, 32), 2, 3),
            ("matched", 4, (8, 16, 32), 2, 3),
            ("mixed", 2, (8, 16, 32), 2, 2),
        ]:
            case = out / f"{label}-np{ranks}-direction{direction}"
            case.mkdir()
            par = base.replace("refine_max_level=2", f"refine_max_level={level}")
            for axis, old, new in zip(range(1, 4), (12, 16, 48), cells):
                par = par.replace(f"domain_nx{axis}={old}", f"domain_nx{axis}={new}")
            par += f"\n&transport_test\n direction={direction}\n refined_patches={patches}\n/\n"
            (case / "amrvac.par").write_text(par)
            with (case / "run.log").open("w") as log:
                subprocess.run(
                    [sys.executable, str(guard), "--out", "memory", "--limit-mib", "1024",
                     "--timeout-seconds", "180", "--", args.mpi, "--oversubscribe",
                     "-np", str(ranks), str(executable)],
                    cwd=case, stdout=log, stderr=log, env=env, check=True,
                )
            if (case / "memory" / "command.err").stat().st_size:
                raise RuntimeError(f"Unexpected runtime stderr in {case}")
            history = np.loadtxt(case / "transport-history.txt")
            if history[-1, 0] != 60 or abs(history[-1, 1] - 0.6) > 1e-12:
                raise AssertionError(f"Incomplete advection run in {case}")
            runs.append(case)
        for control, candidate in [(runs[0], runs[1]), (runs[1], runs[2])]:
            comparison[candidate.name] = {
                "max_cell_difference": compare(control, candidate, "transport", 4),
                "max_face_difference": compare(control, candidate, "faces", 8),
            }
        (out / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")
    print(out / "comparison.json")


if __name__ == "__main__":
    main()
