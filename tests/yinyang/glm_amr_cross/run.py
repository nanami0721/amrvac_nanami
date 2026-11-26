#!/usr/bin/env python3
"""Run checked cross-patch GLM AMR reconstructions and restriction tests."""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
FIELDS = (
    "it", "level", "ghosts", "coarse", "same", "fine",
    "mixed_corners", "special_callbacks", "closure_corners", "migrated_donor_sites",
)


def run(executable, workdir, launcher, guard, name, ranks, text, restart=None):
    directory = workdir / name
    directory.mkdir()
    (directory / "amrvac.par").write_text(text)
    command = [*launcher, "-np", str(ranks), str(executable)]
    if restart:
        command += ["-if", str(restart)]
    logs = directory
    if guard:
        logs = directory / "memory"
        command = [
            sys.executable, str(guard), "--out", str(logs),
            "--limit-mib", "1024", "--timeout-seconds", "120", "--", *command,
        ]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    with (directory / "command.log").open("w") as out, (directory / "command.err").open("w") as err:
        process = subprocess.run(command, cwd=directory, env=env, stdout=out, stderr=err, timeout=150)
    checks = []
    visits = 0
    with (logs / "command.log").open() as stdout:
        for line in stdout:
            if line.startswith("CROSS AMR "):
                *counts, error = line.split()[2:]
                checks.append(dict(zip(FIELDS, map(int, counts)), max_error=float(error)))
            if "Finished AMRVAC" in line:
                visits += 1
    if process.returncode or visits != 2 or (logs / "command.err").stat().st_size:
        raise RuntimeError(f"{name} did not finish cleanly. See {logs}")
    if not checks or not all(row["max_error"] <= 2e-12 for row in checks):
        raise RuntimeError(f"{name} failed the ghost comparison. See {logs}")
    print(f"{name} passed", flush=True)
    return checks


def require_coverage(rows, highest, special=False):
    assert max(row["level"] for row in rows) == highest, rows
    for field in ("coarse", "same", "fine", "mixed_corners"):
        assert sum(row[field] for row in rows) > 0, (field, rows)
    if special:
        assert sum(row["special_callbacks"] for row in rows) > 0, rows
        assert sum(row["closure_corners"] for row in rows) > 0, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=lambda value: Path(value).resolve(strict=True))
    parser.add_argument("workdir", type=lambda value: Path(value).resolve())
    parser.add_argument("--mpiexec", default=["mpirun"], type=shlex.split)
    parser.add_argument("--memory-guard", type=lambda value: Path(value).resolve(strict=True))
    args = parser.parse_args()
    source = HERE.parents[2]
    if args.workdir == source or source in args.workdir.parents:
        parser.error("the run directory must be outside the source tree")
    args.workdir.mkdir(parents=True)
    par = (HERE / "amrvac.par").read_text()
    case_args = (args.executable, args.workdir, args.mpiexec, args.memory_guard)
    report = {}

    report["nonuniform"] = run(*case_args, "nonuniform", 4, par)
    require_coverage(report["nonuniform"], 4)
    assert report["nonuniform"][0]["level"] == 1, report["nonuniform"]
    assert report["nonuniform"][-1]["level"] < 4, report["nonuniform"]
    assert sum(row["migrated_donor_sites"] for row in report["nonuniform"]) > 0, report["nonuniform"]

    constant = par.replace("refine_max_level=4", "refine_max_level=3")
    constant += "\n&cross_test\n uniform_state=T\n initial_refine=T\n/\n"
    report["constant"] = run(*case_args, "constant", 2, constant)
    require_coverage(report["constant"], 3)
    assert any(row["it"] == 0 and row["level"] == 3 for row in report["constant"])

    stretched = par.replace("&meshlist", "&meshlist\n stretch_dim(1)='uni'\n qstretch_baselevel(1)=1.04d0")
    report["stretched"] = run(*case_args, "stretched", 4, stretched)
    require_coverage(report["stretched"], 4)

    symmetric = stretched.replace("stretch_dim(1)='uni'", "stretch_dim(1)='symm'\n nstretchedblocks_baselevel(1)=2")
    report["symmetric"] = run(*case_args, "symmetric", 4, symmetric)
    require_coverage(report["symmetric"], 4)

    special = stretched.replace("typeboundary_min1=8*'cont'", "typeboundary_min1=8*'special'")
    special = special.replace("&mhd_list", "&mhd_list\n mhd_internal_e=T")
    special += "\n&cross_test\n special_inner=T\n/\n"
    report["special"] = run(*case_args, "special", 4, special)
    require_coverage(report["special"], 4, special=True)

    restart = args.workdir / "nonuniform" / "cross_{patch}0003.dat"
    report["restart"] = run(*case_args, "restart", 4, par, restart)
    require_coverage(report["restart"], 4)
    assert report["restart"][0]["it"] == 3 and report["restart"][0]["level"] == 4
    (args.workdir / "results.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
