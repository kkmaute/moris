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

## Testing Guidelines

Unit tests use Catch via `include/catch.hpp` and are registered with CTest. Add regression cases beside the affected module; register new sources in its test `CMakeLists.txt`. Test names commonly follow `<MODULE>-test-<N>-procs`. Exercise MPI process counts when changing communication behavior. No repository-wide numerical coverage threshold is established.

## Commit & Pull Request Guidelines

History uses short descriptive subjects, sometimes prefixed with module tags such as `[GEN]`, `[OPT]`, or `[GEN/UT]`. Follow that pattern and keep changes focused. Describe the problem, resulting behavior, affected options, and exact validation commands/results in the PR. Link relevant issues and identify untested configurations. Include numerical comparisons or visual results when they substantiate simulation changes.
