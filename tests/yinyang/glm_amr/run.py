#!/usr/bin/env python3
"""Run the patch-local GLM AMR regression with an externally built executable."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def run(case, ranks, source, expected=None, restart=None):
    run_dir = args.workdir / case
    run_dir.mkdir()
    (run_dir / "amrvac.par").write_text(source)
    cmd = [*args.mpiexec, "-np", str(ranks), str(args.executable)]
    if restart is not None:
        cmd += ["-if", str(restart)]
    logs = run_dir
    if args.memory_guard is not None:
        logs = run_dir / "memory"
        cmd = [
            sys.executable, str(args.memory_guard), "--out", str(logs),
            "--limit-mib", "1024", "--timeout-seconds", "120", "--", *cmd,
        ]
    with (run_dir / "command.log").open("w") as out, (run_dir / "command.err").open("w") as err:
        proc = subprocess.run(cmd, cwd=run_dir, env=env, stdout=out, stderr=err, timeout=150)
    stdout = (logs / "command.log").read_text()
    stderr = logs / "command.err"
    if expected is not None:
        if proc.returncode == 0 or expected not in stdout:
            raise RuntimeError(f"{case} failed to reject its input. See {logs}")
        print(f"{case} rejected", flush=True)
        return []
    finished = stdout.count("Finished AMRVAC")
    if proc.returncode != 0 or finished != (1 if case == "native" else 2) or stderr.stat().st_size:
        raise RuntimeError(f"{case} did not finish cleanly. See {logs}")
    report = []
    for line in stdout.splitlines():
        if line.startswith("YY AMR "):
            step, ghosts, yin, yang, migrated, error = line[7:].split()
            report.append(
                dict(it=int(step), ghosts=int(ghosts), yin_fine=int(yin), yang_fine=int(yang),
                     moved_yin_donors=int(migrated), max_ghost_error=float(error))
            )
        elif line.startswith("NATIVE AMR "):
            step, fine = line[11:].split()
            report.append(dict(it=int(step), fine=int(fine)))
    print(f"{case} passed", flush=True)
    return report


def check_steps(report, steps, fine):
    assert [row["it"] for row in report] == steps, report
    assert [row["yin_fine"] for row in report] == fine, report
    assert all(row["yang_fine"] == 0 for row in report), report
    assert all(row["ghosts"] > 0 and row["max_ghost_error"] <= 2e-13 for row in report), report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=lambda value: Path(value).resolve(strict=True))
    parser.add_argument("workdir", type=lambda value: Path(value).resolve())
    parser.add_argument("--mpiexec", default=["mpirun"], type=shlex.split)
    parser.add_argument("--memory-guard", type=lambda value: Path(value).resolve(strict=True))
    args = parser.parse_args()
    if not os.access(args.executable, os.X_OK):
        parser.error("the regression executable must be executable")
    if HERE.parents[2] in args.workdir.parents or HERE.parents[2] == args.workdir:
        parser.error("the run directory must be outside the source tree")
    args.workdir.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    par = (HERE / "amrvac.par").read_text()
    results = {}
    results["dynamic"] = run("dynamic", 4, par)
    check_steps(results["dynamic"], [0, 1, 2, 3], [0, 8, 8, 0])
    assert results["dynamic"][1]["moved_yin_donors"] > 0, results["dynamic"]
    assert results["dynamic"][2]["moved_yin_donors"] == 0, results["dynamic"]
    assert results["dynamic"][3]["moved_yin_donors"] > 0, results["dynamic"]

    results["initial"] = run("initial", 4, par + "\n&amr_test\n initial_refine=T\n/\n")
    check_steps(results["initial"], [0, 1, 2, 3], [8, 8, 8, 0])

    restart = args.workdir / "dynamic" / "amr_{patch}0001.dat"
    results["restart"] = run("restart", 4, par, restart=restart)
    check_steps(results["restart"], [1, 2, 3], [8, 8, 0])

    results["fixed"] = run("fixed", 4, par.replace("refine_max_level=2", "refine_max_level=1"))
    check_steps(results["fixed"], [0], [0])

    results["native"] = run("native", 2, par + "\n&amr_test\n native=T\n/\n")
    assert results["native"] == [dict(it=it, fine=fine) for it, fine in enumerate([0, 8, 8, 0])], results["native"]

    run("above-native-level", 2, par.replace("refine_max_level=2", "refine_max_level=21"),
        expected="Reset nlevelshi and recompile!")
    run("forced-io", 2, par.replace(" autoconvert=", " level_io=2\n autoconvert="),
        expected="Yin-Yang AMR does not support forced output refinement levels")
    ct = par.replace("typedivbfix='glm'", "typedivbfix='ct'\n type_ct='uct_hll'")
    ct = ct.replace("&methodlist", "&methodlist\n ct_stage_split=T\n yy_ct_experimental=T")
    run("ct", 2, ct, expected="Yin-Yang CT requires a fixed level-one mesh")

    (args.workdir / "results.json").write_text(json.dumps(results, indent=2) + "\n")
