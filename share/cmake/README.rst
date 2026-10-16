MORIS CMake build configuration
==============================

The top-level ``CMakeLists.txt`` coordinates the build. Configuration is split
by responsibility so changes to compiler policy do not require editing package
selection or installation rules:

* ``moris_options.cmake`` declares the public cache options and their defaults.
* ``moris_paths.cmake`` defines the established build and installation layout.
* ``moris_features.cmake`` validates numerical backends and defines feature macros.
* ``moris_compiler.cmake`` handles C++17, warnings, PIC, IPO/LTO, and instrumentation.
* ``moris_runtime.cmake`` locates the runtime libraries needed by existing TPLs.
* ``moris_packages.cmake`` selects modules, tests, examples, and Valgrind commands.
* ``moris_modules.cmake`` lists module selectors and documented ALL exclusions.
* ``moris_install.cmake`` exports the top-level package configuration.
* ``utilities/moris_package.cmake`` generates every installed component package
  from ``morisComponentConfig.cmake.in``.
* ``utilities/moris_shared.cmake`` assembles the optional monolithic shared library.
* ``utilities/moris_testing.cmake`` registers unit tests and their MPI process counts.

Compatibility
-------------

CMake 3.17 remains the minimum. The policy range is updated through CMake 4.4.
Existing ``MORIS_*`` and ``BUILD_*`` option names remain supported.
``BUILD_ALL`` and ``BUILD_NONE`` retain their existing one-shot cache behavior.

Module selection
----------------

``BUILD_ALL=ON`` enables maintained modules (including MIG), the main executable,
and examples. GUI requires Qt6 and remains opt-in; an explicit ``BUILD_GUI=ON``
is preserved. Tutorials and the standalone SDF executable also remain opt-in.
MOD is disabled in this checkout, and no standalone STK module exists, so ALL
disables those two selectors.
``BUILD_STK`` does not control the STK dependency supplied by Trilinos.

``BUILD_NONE=ON`` disables all module executable/test selectors, the main
executable, the SDF executable, tutorials, and examples. Required core libraries
are still built.
Both controls reset to OFF after configuration; enabling both is an error.
Individual ``BUILD_<module>`` selectors can be changed on later configurations.
``MORIS_USE_TESTS`` remains the master switch for unit-test compilation.

``BUILD_SDF_EXE=ON`` builds and installs the standalone ``sdf`` executable and
selects its required libraries. It defaults to OFF, independently of
``BUILD_SDF`` (SDF unit tests) and ``BUILD_ALL``. The SDF library remains available
to dependent modules when the executable is disabled.

``MORIS_USE_EXAMPLES`` is the canonical examples option. ``BUILD_EXA`` is a
synchronized legacy alias: it initializes a fresh cache if the canonical option
is absent, and subsequent alias edits are accepted when the canonical value has
not changed. A changed canonical option takes precedence.

Installed packages
------------------

All component packages use the same generation helper, including header-only
components and XTK. It handles component dependencies, static dependency cycles,
OpenMP discovery, and the shared-library import while preserving package names
and installation directories. Module files retain target and header installation.

``find_package(moris REQUIRED COMPONENTS XTK)`` selects installed library
components; omitting COMPONENTS imports all installed library packages. Missing
required components fail discovery, while missing OPTIONAL_COMPONENTS are allowed.
Repeated lookups may add components. Executables and examples are not advertised
as installed library components.

Build compatibility
-------------------

Select compilers with a toolchain file, ``CMAKE_C_COMPILER`` /
``CMAKE_CXX_COMPILER``, or ``CC`` / ``CXX`` before the first configuration.
The legacy ``MORIS_C_COMPILER`` / ``MORIS_CXX_COMPILER`` aliases are accepted
when the corresponding standard CMake compiler variable is unset. Compiler
selection is never changed by dependency discovery after ``project()``.

Standard CMake flags are preserved. ``MORIS_CXX_FLAGS`` supplies additional
project flags. Single-configuration builds default to Release unless a build
type is supplied. ``MORIS_HAVE_SYMBOLIC`` retains its Debug selection;
``MORIS_HAVE_DEBUG`` independently controls MORIS runtime checks.

``MORIS_USE_LTO`` uses CMake's IPO capability check and compiler-specific tools
for optimized configurations. Static libraries remain position-independent
because input-file plugins can link them. MPI remains a build dependency in
serial mode because current logger headers use MPI directly; the
``MORIS_HAVE_PARALLEL`` definition still controls MORIS parallel behavior.

Static archives keep their per-package output directories. Shared mode produces
``lib/libmoris.so`` (or the platform equivalent) while preserving public component
target names. Installation retains the existing package and header directories.
External libraries are not bundled into the installation.

Set the environment variable ``ZLIB_DIR`` to the zlib installation prefix used
by the dependency stack (for example, the Spack zlib-ng prefix). Runtime-library
discovery searches that prefix and its ``lib`` / ``lib64`` directories for
``libz``. A nonempty ``ZLIB_DIR`` overrides the cached ``MORIS_LIB_Z`` on each
configuration; an invalid prefix is an error rather than a fallback to system
zlib. With the variable unset, the existing cache and normal search are retained.

Runtime requirements are carried by the exported ``moris_runtime`` target and
linked explicitly by MORIS libraries. Unrelated targets no longer inherit them
through directory-wide linking. ``OPENSSL_DIR`` selects SSL and crypto from a
single prefix; ``DL_DIR``, ``RT_DIR``, and ``INTEL_RUNTIME_DIR`` provide optional
overrides for the loader, realtime, and Intel runtimes. Explicit prefixes refresh cached library
locations in both compiler branches. The canonical crypto cache entry is
``MORIS_LIB_Crypto``; the old mixed-case spelling is migrated automatically.
Fortran discovery preserves explicit overrides and otherwise asks the active
compiler for its runtime before searching system locations.

``MORIS_USE_GPERFTOOLS`` controls profiling throughout the build and installed
packages. The old ``USE_GPERFTOOLS`` option initializes it only when the canonical
option is absent. The ``moris::gperftools`` target propagates headers, allocation
compiler flags, and the combined ``tcmalloc_and_profiler`` library. Its finder
and wrapper are installed so downstream consumers can rediscover the dependency.

Unit tests always retain one-process coverage when
``MORIS_HAVE_PARALLEL_TESTS=OFF``. Enabling it restores each module's established
process counts, including HMR's 12-process case. Tests declare their processor
requirements to CTest. This applies to unit tests; application examples retain
their existing registration policy.

Formatting and file layout
--------------------------

The repository's ``.cmake-format.py`` defines four-space indentation, lowercase
commands, and a 100-column layout. Use cmake-format 0.6.13 for reproducible
formatting, including package templates::

    cmake-format -i projects/COM/CMakeLists.txt projects/COM/src/CMakeLists.txt
    cmake-format --check share/cmake/morisComponentConfig.cmake.in

Keep module entry points in source, tests, then tutorials order. Within library
files, group headers and sources, dependencies, target setup, package
configuration, and installation. Test and example files should group inputs,
dependencies, targets, then test registration. Preserve ordering where CMake
scope or dependency initialization requires it.

Use short section comments and explain non-obvious requirements beside the
relevant command. New repository-owned files use the copyright and MIT license
header from neighboring files; retain third-party copyright and license notices.
Formatting must not reorder source lists, dependencies, or commands.

Validation
----------

The focused tests require Python 3, CMake, and a C++17 compiler. They do not
require the MORIS TPL stack::

    python3 share/cmake/tests/test_build_configuration.py
    python3 share/cmake/tests/test_shared_library.py
    python3 share/cmake/tests/test_module_selection.py
    python3 share/cmake/tests/test_package_configuration.py
    python3 share/cmake/tests/test_runtime_configuration.py
    python3 share/cmake/tests/test_test_registration.py

The first test also accepts a CMake executable as its positional argument.
It exercises user flags, Debug/Release and symbolic modes, PIC, IPO, plugin
includes, and Ninja Multi-Config when Ninja is available. The shared-library
test checks cyclic component dependencies, runtime linking, relocated exports,
and the static build path.
The module-selection test covers ALL/NONE, exclusions, cache transitions, and
the examples alias. The package-configuration test builds and runs relocated
static/shared consumers with OpenMP off and on, including cyclic dependencies,
repeated lookups, and missing required/optional components. Its OpenMP cases
require a compiler with OpenMP support.
It also checks profiling requirements using a small stand-in library by default;
set ``MORIS_TEST_GPERFTOOLS_PREFIX`` to an installed gperftools prefix to exercise
the real profiler. The runtime test covers explicit prefixes, cache changes,
Intel/non-Intel discovery, and absence of directory-wide link requirements.
The registration test covers serial/MPI transitions and CTest processor counts.

``.github/workflows/cmake-regression.yml`` runs these six tests and
``python3 share/cmake/tests/check_format.py`` on a runner without the MORIS TPL
stack. The integration workflows explicitly build debug/shared and
optimized/static configurations on the configured self-hosted runners.

With the MORIS dependencies available, a reduced integration build can be run
without the expensive application examples::

    cmake -S . -B build-check -DMORIS_USE_EXAMPLES=OFF -DBUILD_CNT=ON -DBUILD_COM=ON
    cmake --build build-check --parallel
    ctest --test-dir build-check -R '^(CNT|COM)-test-1-procs$' --output-on-failure

Repeat with ``-DMORIS_HAVE_SHARED=ON`` to exercise the shared installation.
