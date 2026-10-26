#!/usr/bin/env python3
"""Run Cartesian elasticity BDDC studies; retain each log and aggregate CSV results."""

import argparse
import csv
import itertools
import json
import os
from pathlib import Path
import platform
import re
import shlex
import subprocess


def cases(args):
    for subdomains in args.subdomains:
        ranks = subdomains ** args.dim
        elements = args.elements
        if args.study in ("strong", "jump"):
            if args.global_elements % subdomains:
                raise ValueError("global-elements must be divisible by every subdomain count")
            elements = [args.global_elements // subdomains]
        sigmas = args.sigmas if args.study == "jump" else [1.0]
        for count, sigma, constraints, solver, repeat in itertools.product(
                elements, sigmas, args.constraints, args.solvers, range(args.repeat)):
            if count < 2:
                raise ValueError("each subdomain direction needs at least two elements")
            if args.study == "jump" and subdomains * count % 4:
                raise ValueError("jump meshes must resolve interfaces at 1/4 and 3/4")
            # The direct AIJ baseline has no corner/edge constraints.
            if solver == "aij" and constraints != args.constraints[0]:
                continue
            yield ranks, subdomains, count, sigma, constraints, solver, repeat


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--dim", required=True, type=int, choices=[2, 3])
    parser.add_argument("--study", choices=["weak", "strong", "refine", "jump"], default="weak")
    parser.add_argument("--subdomains", nargs="+", type=int, default=[2, 4],
                        help="subdomains per axis; MPI ranks = this value ** dim")
    parser.add_argument("--elements", nargs="+", type=int, default=[8],
                        help="elements per subdomain direction (weak/refine)")
    parser.add_argument("--global-elements", type=int, default=24,
                        help="fixed global elements per direction (strong/jump)")
    parser.add_argument("--sigmas", nargs="+", type=float, default=[0.001, 0.01, 1, 100, 1000, 10000])
    parser.add_argument("--constraints", nargs="+", choices=["c", "ce"], default=["c", "ce"])
    parser.add_argument("--solvers", nargs="+", choices=["bddc", "aij"], default=["bddc"])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--rtol", type=float, default=1.e-6)
    parser.add_argument("--launcher", default="mpirun", help="launcher and extra arguments, e.g. 'srun --exclusive'")
    parser.add_argument("--petsc-options", default="", help="additional PETSc options, appended to the driver defaults")
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--output", type=Path, default=Path("bddc-results"))
    parser.add_argument("--dry-run", action="store_true", help="print commands without creating files or launching MPI")
    args = parser.parse_args()

    if args.repeat < 1 or any(s < 1 for s in args.subdomains):
        parser.error("repeat and subdomain counts must be positive")
    if args.study == "refine" and len(args.subdomains) != 1:
        parser.error("refine studies use one fixed --subdomains value and multiple --elements")
    if args.study == "weak" and len(args.elements) != 1:
        parser.error("weak studies use one fixed --elements value")
    try:
        configurations = list(cases(args))
    except ValueError as error:
        parser.error(str(error))

    executable = args.exe.resolve()
    output = args.output.resolve()
    if not args.dry_run:
        if not executable.is_file():
            parser.error("benchmark executable does not exist: " + str(executable))
        # Preserve earlier measurements: each study gets its own output directory.
        output.mkdir(parents=True, exist_ok=False)
        metadata = {"arguments": vars(args).copy(), "platform": platform.platform(),
                    "executable": str(executable), "threads": {}}
        metadata["arguments"]["exe"] = str(args.exe)
        metadata["arguments"]["output"] = str(args.output)

    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        environment.setdefault(name, "1")
        if not args.dry_run:
            metadata["threads"][name] = environment[name]
    if not args.dry_run:
        (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    for index, (ranks, subdomains, elements, sigma, constraints, solver, repeat) in enumerate(configurations):
        stem = "{:04d}_{}D_n{}_m{}_{}_{}_sigma{}_r{}".format(
            index, args.dim, ranks, elements, solver, constraints, sigma, repeat)
        result = output / (stem + ".csv")
        case_name = "BDDC_Plane_Stress" if args.dim == 2 else "BDDC_Elasticity_3D"
        command = shlex.split(args.launcher) + ["-n", str(ranks), str(executable), case_name]
        # The executable parses these options into the same globals used by
        # regression cases and the dynamically loaded EXA input libraries.
        command += ["--benchmark", "--elements", str(elements),
                    "--sigma", str(sigma), "--solver", solver,
                    "--use-edges", "true" if constraints == "ce" else "false",
                    "--rtol", str(args.rtol), "--csv", str(result)]
        if args.study == "jump":
            command += ["--jump"]
        if args.petsc_options:
            # Bind the value to its option so Clara accepts a leading PETSc dash.
            command += ["--petsc-options=" + args.petsc_options]
        print(" ".join(shlex.quote(arg) for arg in command), flush=True)
        if args.dry_run:
            continue

        directory = output / stem
        directory.mkdir()
        (output / (stem + ".command.json")).write_text(json.dumps(
            {"command": command, "cwd": str(directory)}, indent=2) + "\n")
        log_path = output / (stem + ".log")
        with log_path.open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=str(directory),
                           env=environment, timeout=args.timeout, check=True)
        with result.open(newline="") as source:
            rows = list(csv.DictReader(source))
        if len(rows) != 1:
            raise RuntimeError("missing workflow measurement in " + str(result))
        log_text = log_path.read_text()
        iterations = re.findall(r"Linear solve converged due to \w+ iterations (\d+)", log_text)
        if not iterations:
            raise RuntimeError("PETSc reported no converged linear solve; see " + str(log_path))
        residuals = re.findall(r"\|\|r\(i\)\|\|/\|\|b\|\|\s+([0-9.eE+-]+)", log_text)
        diagnostics = {"first_linear_iterations": iterations[0],
                       "total_linear_iterations": sum(int(value) for value in iterations),
                       "linear_solves": len(iterations),
                       "last_relative_residual": residuals[-1] if residuals else ""}
        if residuals and float(residuals[-1]) > max(1.e-10, 1.05 * args.rtol):
            raise RuntimeError("PETSc true residual exceeds tolerance; see " + str(log_path))

        # PETSc's ASCII event table reports maximum time across ranks. These are
        # event timings (possibly nested), not disjoint workflow wall-clock stages.
        petsc_log = directory / "petsc.log"
        text = petsc_log.read_text()
        for event, key in [("KSPSolve", "petsc_ksp_seconds"), ("PCSetUp", "petsc_pc_setup_seconds")]:
            match = re.search(r"^" + event + r"\s+(\d+)\s+\S+\s+([0-9.eE+-]+)", text, re.M)
            diagnostics[key] = match.group(2) if match else ""
        row = rows[0]
        row.update(diagnostics)
        row["repeat"] = repeat
        row["log"] = stem + ".log"
        with (output / "results.csv").open("a", newline="") as aggregate:
            writer = csv.DictWriter(aggregate, fieldnames=list(row.keys()))
            if index == 0:
                writer.writeheader()
            writer.writerow(row)


if __name__ == "__main__":
    main()
