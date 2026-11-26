#!/usr/bin/env python3
"""Build and run the Yin-Yang CT-AMR curl commutator outside the checkout."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
FIELDS = (
    "area", "commutator_rms", "edge_then_curl_rms", "curl_then_face_rms",
    "native_curl_rms", "reference_rms", "commutator_peak", "edge_then_curl_peak",
    "curl_then_face_peak", "native_curl_peak", "rate_divb_peak", "edge_change_peak",
    "common_reference_rms",
)
DIRECTIONS = ("fine_to_coarse", "coarse_to_fine")
BLOCK_WIDTHS = (4, 8, 16)
BASE_MESHES = {
    # Block shape, base-patch cell counts, lower theta, and outer phi in turns.
    "wide-overlap": ((4, 4, 4), (8, 12, 32), 0.09765625, 0.40234375),
    "standard": ((4, 4, 16), (8, 16, 48), 0.109375, 0.390625),
}


def guarded(guard, root, name, command, directory, env, limit, *, check_stderr=False):
    memory = root / f"{name}-memory"
    with (root / f"{name}.log").open("w") as log:
        status = subprocess.run(
            [sys.executable, str(guard), "--out", str(memory), "--limit-mib", str(limit),
             "--timeout-seconds", "180", "--", *map(str, command)],
            cwd=directory, env=env, stdout=log, stderr=log,
        )
    if status.returncode or (check_stderr and (memory / "command.err").stat().st_size):
        raise RuntimeError(f"{name} did not finish cleanly. See {memory}")
    return memory


def parameters(width, geometry):
    blocks, domain, theta_min, phi_max = BASE_MESHES[geometry]
    settings = {
        "xprobmin2": f"{theta_min:.8f}d0", "xprobmax2": f"{0.5 - theta_min:.8f}d0",
        "xprobmin3": f"{-phi_max:.8f}d0", "xprobmax3": f"{phi_max:.8f}d0",
    }
    for prefix, counts in (("block_nx", blocks), ("domain_nx", domain)):
        for dim, count in enumerate(counts, start=1):
            settings[f"{prefix}{dim}"] = width * count // 4
    lines = (HERE / "amrvac.par").read_text().splitlines()
    changed = set()
    for index, line in enumerate(lines):
        name, _, _ = line.strip().partition("=")
        if name in settings:
            lines[index] = f" {name}={settings[name]}"
            changed.add(name)
    if changed != settings.keys():
        raise ValueError(f"Missing expected parameters {settings.keys() - changed}")
    return "\n".join(lines) + "\n"


def run_case(args, executable, env, width, ranks):
    name = f"block{width}-np{ranks}"
    directory = args.workdir / name
    directory.mkdir()
    (directory / "amrvac.par").write_text(parameters(width, args.geometry))
    memory = guarded(
        args.guard, args.workdir, name, [*args.mpiexec, "-np", ranks, executable],
        directory, env, args.run_memory_mib, check_stderr=True,
    )

    result = {}
    finished = 0
    with (memory / "command.log").open() as stdout:
        for line in stdout:
            if "Finished AMRVAC" in line:
                finished += 1
            if not line.startswith("CT_COMMUTATOR "):
                continue
            _, patch, count, *values = line.split()
            direction = DIRECTIONS[int(patch)]
            if direction in result or len(values) != len(FIELDS):
                raise RuntimeError(f"Unexpected CT commutator report in {memory}")
            row = dict(zip(FIELDS, map(float, values)))
            row["count"] = int(count)
            result[direction] = row
    if finished != 2 or set(result) != set(DIRECTIONS):
        raise RuntimeError(f"Incomplete directional result in {memory}")
    for direction, row in result.items():
        values = tuple(row.values())
        if not all(math.isfinite(value) and value >= 0 for value in values):
            raise RuntimeError(f"Degenerate {direction} check in {memory}")
        if (row["count"] <= 0 or row["area"] <= 0 or row["commutator_rms"] < 1e-12
                or row["common_reference_rms"] < 1e-2 or row["edge_change_peak"] < 1e-10):
            raise RuntimeError(f"Unresolved manufactured signal in {memory}")
        if row["rate_divb_peak"] > 1e-10:
            raise RuntimeError(f"Non-solenoidal native curl in {memory}")
        row["relative_commutator_rms"] = row["commutator_rms"] / row["common_reference_rms"]
    print(name, "completed", flush=True)
    return result


def compare_ranks(reference, repeated):
    for direction in DIRECTIONS:
        for field, value in reference[direction].items():
            other = repeated[direction][field]
            if not math.isclose(value, other, rel_tol=1e-11, abs_tol=1e-12):
                raise RuntimeError(f"MPI layout changed {direction} {field}: {value} vs {other}")


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=lambda value: Path(value).resolve(strict=True),
                        help="external AMRVAC candidate source copy")
    parser.add_argument("workdir", type=lambda value: Path(value).resolve(),
                        help="new external result directory")
    parser.add_argument("--guard", type=lambda value: Path(value).resolve(strict=True))
    parser.add_argument("--mpiexec", type=shlex.split, default=["mpirun", "--oversubscribe"])
    parser.add_argument("--arch", default="debug")
    parser.add_argument("--geometry", choices=BASE_MESHES, default="wide-overlap")
    parser.add_argument("--run-memory-mib", type=int, default=2048)
    parser.add_argument("--repeat-ranks", type=int, default=4)
    args = parser.parse_args()
    checkout = next((parent for parent in HERE.parents if (parent / ".git").exists()), None)
    for path in (args.source, args.workdir):
        if (path / ".git").exists() or (checkout and (path == checkout or checkout in path.parents)):
            parser.error("use an external source copy and an external run directory")
    if args.guard is None:
        args.guard = args.source / "tests/yinyang/interface_transport/convergence/memory_guard.py"
    if args.repeat_ranks and (args.repeat_ranks <= 2 or args.repeat_ranks % 2):
        parser.error("--repeat-ranks must be zero or an even integer greater than two")
    args.workdir.mkdir(parents=True, exist_ok=False)
    build = args.workdir / "build"
    build.mkdir()
    shutil.copy2(HERE / "mod_usr.t", build / "mod_usr.t")
    env = dict(os.environ, AMRVAC_DIR=str(args.source), OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    guarded(args.guard, args.workdir, "setup", [args.source / "setup.pl", "-d=3",
            f"-arch={args.arch}"], build, env, 1024)
    guarded(args.guard, args.workdir, "build", ["make", "-j1",
            "F90=mpif90 --param ggc-min-expand=5 --param ggc-min-heapsize=4096"],
            build, env, 4096)
    executable = build / "amrvac"
    manifest = {
        "source": str(args.source), "architecture": args.arch, "geometry": args.geometry,
        "executable_sha256": file_hash(executable),
        "source_sha256": {
            str(path.relative_to(args.source)): file_hash(path)
            for path in (*sorted((args.source / "src").rglob("*.t")),
                         args.source / "src/makefile", args.source / f"arch/{args.arch}.defs")
        },
        "fixture_sha256": {name: file_hash(HERE / name)
                           for name in ("amrvac.par", "mod_usr.t", "run.py")},
    }
    (args.workdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    report_path = args.workdir / "results.json"
    report = {"cases": {}, "orders": {}}
    for width in BLOCK_WIDTHS:
        report["cases"][f"block{width}-np2"] = run_case(args, executable, env, width, 2)
        report_path.write_text(json.dumps(report, indent=2) + "\n")
    if args.repeat_ranks:
        middle = BLOCK_WIDTHS[1]
        repeated = f"block{middle}-np{args.repeat_ranks}"
        report["cases"][repeated] = run_case(args, executable, env, middle, args.repeat_ranks)
        report_path.write_text(json.dumps(report, indent=2) + "\n")
        compare_ranks(report["cases"][f"block{middle}-np2"], report["cases"][repeated])
    ratios = [b / a for a, b in zip(BLOCK_WIDTHS, BLOCK_WIDTHS[1:])]
    report["grid_spacing_ratios"] = ratios
    for direction in DIRECTIONS:
        rows = [report["cases"][f"block{width}-np2"][direction] for width in BLOCK_WIDTHS]
        errors = [row["commutator_rms"] for row in rows]
        report["orders"][direction] = [
            math.log(a / b) / math.log(ratio)
            for a, b, ratio in zip(errors, errors[1:], ratios)
        ]
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    if any(order <= 0 for orders in report["orders"].values() for order in orders):
        raise RuntimeError(f"Commutator RMS did not decrease at both refinements. See {report_path}")
    print("CT AMR commutator RMS decreased in both directions", flush=True)


if __name__ == "__main__":
    main()
