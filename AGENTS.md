# Repository Guidelines

## Project Structure & Module Organization

MORIS is a C++ multiphysics and topology-optimization research code. Modules live under `projects/`, including `MRS` (core utilities), `LINALG`, `MTK`, `HMR`, `XTK`, `GEN`, `FEM`, and `SOL`. Modules generally separate implementation in `src/` from tests in `test/`. Application examples and their input files live in `projects/EXA/`. Shared headers are in `include/`, reusable examples in `snippets/`, and build utilities, Spack recipes, installation instructions, and documentation in `share/`.

## Build, Test, and Development Commands

Configure the required third-party libraries using `share/install/` instructions and the applicable Spack environment. Always build outside the source tree:

```bash
cmake -S . -B build-check -DMORIS_USE_EXAMPLES=OFF -DBUILD_CNT=ON -DBUILD_COM=ON
cmake --build build-check --parallel 4
ctest --test-dir build-check -R '^(CNT|COM)-test-1-procs$' --output-on-failure
```

These commands configure, build, and test a reduced core configuration. Add `-DMORIS_HAVE_SHARED=ON` for the shared library. Use `-DMORIS_HAVE_SYMBOLIC=ON` for debug symbols; `MORIS_HAVE_DEBUG` independently controls runtime checks. See `share/cmake/README.rst` for build-system organization.

For CMake changes, run the dependency-light regression checks:

```bash
python3 share/cmake/tests/test_build_configuration.py
python3 share/cmake/tests/test_shared_library.py
```

## Coding Style & Naming Conventions

Use C++17 and the repository's `.clang-format`: four-space indentation and the configured brace/alignment rules. Format changed C++ files with `clang-format -i <file>`. Follow nearby naming conventions: `cl_*.hpp` for classes, `fn_*.hpp` for functions, and `ut_*` or `UT_*` for tests. Existing code commonly prefixes members with `m`, arguments with `a`, and temporaries with `t`. Keep CMake commands lowercase and explain non-obvious dependency or compatibility decisions.

Separate logical steps with blank lines, especially validation, map construction, assembly, metadata setup, and resource cleanup. Keep related declarations and operations together; avoid both dense uninterrupted blocks and blank lines between every statement. Add comments explaining intent, numbering conventions, collective participation, and compatibility workarounds rather than restating individual calls. Use braced, multiline control blocks in new solver code and tests, even though `.clang-format` permits short statements on one line. Preserve nearby method separators and review the result after formatting; automatic formatting does not supply this logical structure.

## EXA Example Setup & Result Checks

Follow neighboring examples under `projects/EXA/` when adding or refactoring an example. Use MORIS modules for mesh generation, model construction, assembly, solvers, and output: typically HMR, GEN/XTK as needed, FEM, SOL/MSI, and VIS. Define the appropriate module parameter-list functions and property callbacks in dynamically loaded input `.cpp` files. Keep the example implementations in those `.cpp` files rather than hiding them in implementation headers or including another `.cpp`. Include the MORIS headers needed by each source file directly.

Define configuration globals in `example_test_case.cpp` and declare the settings needed by the input libraries with `extern`. Set these globals before calling `fn_WRK_Workflow_Main_Interface`. Use globals for case indices, interpolation order, mesh resolution, material settings, and solver variants; do not pass example configuration through environment variables. Prefer the shared EXA main. If performance studies require command-line options, a local main may follow the same MORIS/Catch initialization and populate the same globals. Keep benchmark options separate from fixed regression settings.

Put required solver tuning in the input's SOL parameter lists so the selected problem configuration does not depend on CTest injecting `PETSC_OPTIONS`. Use named preconditioner and linear-algorithm parameters; translate PETSc switches in `cl_DLA_Preconditioner_PETSc.cpp` when no public setter exists. Reserve native PETSc options and `SOL_PETSc_Options` for explicit expert overrides. Check optional factor packages against PETSc's actual capabilities, and pass example-specific backend choices through driver globals.

Inspect active FEM sets, constitutive models, properties, and couplings before classifying a problem. Multiple GEN geometry phases do not necessarily mean multiple active materials. Both checked-in thick-pressure-vessel inputs currently use one active material with Young's modulus and conductivity equal to one. Distinguish those regression cases from a user's separate multi-material input, and document which material contrasts and process counts were actually verified.

Use explicit Catch cases for meaningful configurations, with CTest selecting each case and MPI process count. Add comments to every case explaining the physical setup, solver or constraint variant, behavior exercised, and results checked. Reuse the same physical references across solver variants that solve the same model. Give tests separate working directories when their output filenames could collide. Respect solver lifetime restrictions, particularly PETSc's inability to reinitialize after finalization.

Implement `check_results()` following nearby EXA tests: read VIS Exodus output with `mtk::Exodus_IO_Helper` and compare mesh dimensions, node/element counts, reference-node coordinates, time, displacement, and relevant global IQIs against stored reference values with explicit tolerances. Document how references were obtained and which case and process count each entry represents. Provide a `gPrintReferenceValues` flag for inspecting reference data. Use a designated MPI rank whose output contains the reference node, and distinguish local mesh counts from global totals. VIS reserves nodal fields 0 and 1 for `Node_Id` and `Node_Index`; requested displacement fields follow them. Close an Exodus reader before another checker opens the same file, and keep collective checks outside rank-specific branches.

Use stored numerical references for regression tests rather than rerunning a reference solver inside every case. Arbitrary performance-study meshes may use integral, topology, and finite-field checks instead of references for a different fixed mesh; document that distinction. Treat iteration counts and timings as diagnostics unless the test explicitly establishes them as regression requirements. Update `example_documentation.dox` and any scaling runner when configuration or loading changes, including differences from literature benchmarks.

For solver comparisons, retain the same mesh, loads, material properties, stabilization, and numerical tolerances. Report iterations for each linear solve and their total across Newton corrections; a different preconditioner can change the number of corrections. Separate PETSc setup/solve time from total workflow time and use the same monitoring and thread settings for each comparison. Krylov residual estimates can differ substantially from explicitly computed residuals in ill-conditioned cut-cell systems; retain nonlinear residual corrections and physical field/IQI checks. Clear inherited `PETSC_OPTIONS` when verifying that managed settings are sufficient, and prevent stale output from satisfying a failed workflow's checks.

## External Problem Inputs

The shell wrappers use `share/scripts/moris_input_build.sh` to configure and build the standalone CMake project in `share/cmake/input/`. Keep this workflow free of Python build helpers and source-tree `tmp/input.so` links. Each input has a `.moris-input-build/` cache beside its source, keyed by the selected MORIS build/configuration, and its shared object is produced beside the input source. External inputs must supply their own configuration defaults rather than assume that EXA test-driver globals exist in the normal MORIS executable.

When external-input settings are missing, first verify the `moris.sh` selected by `PATH`, `MORISROOT`, and `MORISBUILDDBG`/`MORISBUILDOPT`. The selected build must contain `input_build/<configuration>/MorisInputConfig.cmake`; reconfigure that existing build to generate the export. Do not accidentally use a different checkout or introduce a placeholder-link fallback.

## PETSc Matrices & BDDC

PETSc matrix storage is selected before matrix construction through the solver warehouse parameter `SOL_PETSc_Mat_Type`. The default is `"aij"`; BDDC (`PCType = "bddc"`) requires `"is"` (MATIS). Preserve the existing AIJ path when changing BDDC support. Configuration and limitations are documented in `projects/SOL/DLA/doc/DLA_BDDC.dox`.

`Matrix_PETSc` lives in `projects/SOL/SOL_CORE/src/`; BDDC setup is in `projects/SOL/DLA/src/cl_DLA_Preconditioner_PETSc.cpp`. Distinguish MORIS DOF IDs, PETSc global IDs, and MATIS subdomain indices. MATIS uses owned-and-shared DOFs for local assembly and the same global vector ownership as AIJ. Each physical element contribution must be assembled once, and MORIS element values use column-major storage.

Overlap maps can contain unused ghost DOFs. BDDC compacts a separate preconditioning matrix while retaining the original matrix and map for reassembly. Keep compaction collective on all ranks when any rank needs it, and preserve both row and column support for nonsymmetric matrices. Retained strongly constrained DOFs receive unit rows and columns and require zero Newton-correction RHS entries; constraint flags must reach every contributing subdomain. Optional `bddc_use_dof_splitting` identifies fields from the requested physical DOF types.

Choose solver symmetry from the assembled coupled operator. The pressure-vessel thermoelastic operator is nonsymmetric because temperature drives displacement through thermal expansion; retain FGMRES and `bddc_symmetric = false`. Use disconnected-subdomain detection and separate physical-field constraints for these cut meshes. Cartesian subdomain-corner metadata is intended for supported nodal box meshes; do not enable it indiscriminately on enriched cut meshes. Ghost-stabilization matrices must reach the same MATIS assembly path as bulk contributions, with each ghost pair contributed once.

Cut cells and ghost penalties can produce unequal interface stiffness even without material jumps. Test `bddc_use_stiffness_scaling = true` before relying on multiplicity averaging. Deluxe scaling supersedes diagonal stiffness weights and has separate sub-Schur assumptions: for the nonsymmetric vessel, set `use_bddc_defaults = false`, `bddc_use_deluxe_scaling = true`, `bddc_sub_schurs_symmetric = false`, and `bddc_sub_schurs_posdef = false`, together with `bddc_symmetric = false`. PETSc 3.23.2 otherwise assumes SPD sub-Schur problems independently of coarse-basis symmetry; enabling deluxe alone caused a Cholesky zero pivot. The native sub-Schur option names use the PC prefix directly, without a `pc_bddc_` infix. Do not apply SPD-only adaptive spectral-constraint recipes to the full nonsymmetric operator. Preserve native defaults for other inputs and verify serial empty-coarse-problem handling separately from parallel MUMPS configurations.

## Testing Guidelines

Unit tests use Catch via `include/catch.hpp` and are registered with CTest. Add regression cases beside the affected module; register new sources in its test `CMakeLists.txt`. Test names commonly follow `<MODULE>-test-<N>-procs`. Exercise MPI process counts when changing communication behavior. No repository-wide numerical coverage threshold is established.

For PETSc/BDDC changes, the dedicated `DLA-BDDC-test-{1,2,4}-procs` cases cover AIJ/MATIS assembly, repeated assembly, constraints, unused ghosts, multivectors, and BDDC solves. PETSc cannot be reinitialized after finalization, so the BDDC Catch cases run in separate processes. With a configured PETSc build, run:

```bash
ctest --test-dir build_dbg -R '^DLA(-BDDC)?-test-(1|2|4)-procs$' --output-on-failure
ctest --test-dir build_dbg -R '^Pressure_Vessel_3D(_Immersed)?-BDDC-(1|4)-procs$' --output-on-failure
```

The pressure-vessel BDDC tests retain the original numerical references and tolerances. Their four-rank configuration requires PETSc built with MUMPS; CMake checks `PETSC_HAVE_MUMPS`, since `MORIS_USE_MUMPS` alone does not establish PETSc capability. Keep serial and parallel PETSc options separate: forcing MUMPS on the serial BDDC coarse problem can fail when that problem is empty.

## Commit & Pull Request Guidelines

History uses short descriptive subjects, sometimes prefixed with module tags such as `[GEN]`, `[OPT]`, or `[GEN/UT]`. Follow that pattern and keep changes focused. Describe the problem, resulting behavior, affected options, and exact validation commands/results in the PR. Link relevant issues and identify untested configurations. Include numerical comparisons or visual results when they substantiate simulation changes.
