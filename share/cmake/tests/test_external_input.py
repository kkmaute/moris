#!/usr/bin/env python3
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
"""Exercise isolated input builds and legacy script entry points without MORIS TPLs."""
import ctypes
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import textwrap

REPOSITORY = Path(__file__).resolve().parents[3]
SCRIPTS = REPOSITORY / "share/scripts"


def run(*args, **kwargs):
    subprocess.run([str(arg) for arg in args], check=True, **kwargs)


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text))


def main():
    placeholder = REPOSITORY / "projects/mains/input_file.cpp"
    placeholder_state = (placeholder.read_bytes(), placeholder.stat().st_mtime_ns, placeholder.is_symlink())
    with tempfile.TemporaryDirectory(prefix="moris external input ") as temporary:
        root = Path(temporary)
        fixture = root / "fixture"
        write(fixture / "CMakeLists.txt", r'''
            cmake_minimum_required(VERSION 3.17...4.4)
            project(InputBuildFixture LANGUAGES CXX)
            set(MORIS moris)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSER_FLAG=11")
            set(CMAKE_CXX_FLAGS_DEBUG "-g -DFIXTURE_DEBUG")
            set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DFIXTURE_RELEASE")
            add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-DFIXTURE_OPTION=13>")
            add_compile_definitions(FIXTURE_DEFINITION=17)
            include_directories("${CMAKE_CURRENT_SOURCE_DIR}/headers")
            # A legacy nonexistent include must not break the imported target.
            include_directories("${CMAKE_CURRENT_SOURCE_DIR}/obsolete")
            add_library(shared_object_file SHARED placeholder.cpp)
            add_executable(moris main.cpp)
            # Exercise the exported target location, not a guessed mains path.
            set_target_properties(moris PROPERTIES
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/runtime files")
            if(SHARED_FIXTURE)
                add_library(moris_shared SHARED payload.cpp)
                target_compile_definitions(shared_object_file PRIVATE SHARED_FIXTURE)
                target_link_libraries(shared_object_file PRIVATE moris_shared)
            endif()
            include("${MORIS_CMAKE_DIR}/utilities/moris_input_build.cmake")
            moris_export_input_build(shared_object_file)
        ''')
        write(fixture / "placeholder.cpp", "// must remain untouched\n")
        write(fixture / "main.cpp", r'''
            #include <cstdlib>
            #include <fstream>
            int main()
            {
                std::ofstream output(std::getenv("FIXTURE_EXECUTION"));
            #ifdef FIXTURE_DEBUG
                output << "Debug";
            #else
                output << "Release";
            #endif
                return output.good() ? 0 : 1;
            }
        ''')
        write(fixture / "payload.cpp", 'extern "C" int payload() { return 19; }\n')
        write(fixture / "headers/fixture.hpp", "#define HEADER_VALUE 23\n")
        problem = root / "problem files"
        source = problem / "test.problem.cpp"
        write(source, r'''
            #include "fixture.hpp"
            #include <type_traits>
            #if __cplusplus < 201703L
            #error Expected C++17
            #endif
            #if !defined(FIXTURE_DEBUG) && !defined(FIXTURE_RELEASE)
            #error Missing configuration flags
            #endif
            #ifdef SHARED_FIXTURE
            extern "C" int payload();
            #endif
            extern "C" int input_value()
            {
                static_assert(std::is_same_v<int, int>);
                int value = USER_FLAG + FIXTURE_OPTION + FIXTURE_DEFINITION + HEADER_VALUE;
            #ifdef SHARED_FIXTURE
                value += payload();
            #endif
                return value;
            }
        ''')
        environment = os.environ.copy()
        environment["MORISROOT"] = str(REPOSITORY)
        environment.pop("MORISBUILDOPT", None)
        environment.pop("MORIS_INPUT_CONFIGURATION", None)
        for configuration, shared in [("Debug", False), ("Release", True)]:
            build = root / ("selected " + configuration)
            run("cmake", "-S", fixture, "-B", build,
                "-DCMAKE_BUILD_TYPE=" + configuration,
                "-DSHARED_FIXTURE=" + ("ON" if shared else "OFF"),
                "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=" + ("ON" if configuration == "Debug" else "OFF"),
                "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION_DEBUG=OFF",
                "-DMORIS_CMAKE_DIR=" + str(REPOSITORY / "share/cmake"))
            targets = ["moris", "moris_shared"] if shared else ["moris"]
            run("cmake", "--build", build, "--target", *targets)
            marker = build / "lib/input_file.so"
            write(marker, "untouched legacy output")
            selected_state = {str(f): (f.stat().st_size, f.stat().st_mtime_ns)
                              for f in build.rglob("*") if f.is_file()}
            run(SCRIPTS / "create_shared_object.sh", problem, build, source.name, env=environment)
            assert selected_state == {str(f): (f.stat().st_size, f.stat().st_mtime_ns)
                                      for f in build.rglob("*") if f.is_file()}
            # Use a distinct pathname to avoid dlopen reusing an earlier handle.
            loaded = root / (configuration + ".so")
            shutil.copyfile(source.with_suffix(".so"), loaded)
            library = ctypes.CDLL(str(loaded))
            assert library.input_value() == 64 + (19 if shared else 0)
            # Preserve the original Python helper's cache naming across migration.
            identity = hashlib.sha256((str(build.resolve()) + "\n" + configuration).encode()).hexdigest()[:16]
            assert (problem / ".moris-input-build/test.problem" / identity / "CMakeCache.txt").is_file()
            assert marker.read_text() == "untouched legacy output"
            assert (fixture / "placeholder.cpp").read_text() == "// must remain untouched\n"
            assert not (fixture / "input_file.locked").exists()

        # Inspect the actual external target's initialized IPO properties. A
        # Debug override must survive even when the global IPO default is ON.
        probe = root / "ipo probe"
        write(probe / "CMakeLists.txt", '''
            cmake_minimum_required(VERSION 3.17...4.4)
            include("${MORIS_INPUT_CONFIG}")
            project(InputIPOProbe LANGUAGES CXX)
            add_subdirectory("${MORIS_EXTERNAL_PROJECT}" input)
            get_target_property(global_ipo moris_external_input INTERPROCEDURAL_OPTIMIZATION)
            get_target_property(debug_ipo moris_external_input INTERPROCEDURAL_OPTIMIZATION_DEBUG)
            file(WRITE "${CMAKE_BINARY_DIR}/observed.txt" "${global_ipo};${debug_ipo}")
        ''')
        run("cmake", "-S", probe, "-B", root / "ipo probe build",
            "-DCMAKE_BUILD_TYPE=Debug",
            "-DMORIS_INPUT_CONFIG=" + str(root / "selected Debug/input_build/Debug/MorisInputConfig.cmake"),
            "-DMORIS_EXTERNAL_PROJECT=" + str(REPOSITORY / "share/cmake/input"),
            "-DMORIS_INPUT_SOURCE=" + str(source),
            "-DMORIS_INPUT_OUTPUT_DIRECTORY=" + str(problem))
        assert (root / "ipo probe build/observed.txt").read_text() == "ON;OFF"

        assert len(list((problem / ".moris-input-build/test.problem").iterdir())) == 2
        # The main wrapper must resolve its siblings without PATH entries and
        # pass absolute problem/build paths and the optional prefix to MPI.
        debug_build = root / "selected Debug"
        # Missing exports must identify the actual build/configuration and give
        # a usable reconfiguration command, rather than placeholder paths.
        config = debug_build / "input_build/Debug/MorisInputConfig.cmake"
        saved_config = config.with_suffix(".saved")
        config.rename(saved_config)
        try:
            missing = subprocess.run([str(SCRIPTS / "create_shared_object.sh"), str(problem),
                                      str(debug_build), source.name], env=environment,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            assert missing.returncode != 0
            diagnostic = missing.stdout.decode()
            assert str(config) in diagnostic
            assert 'cmake -S "{}" -B "{}"'.format(fixture, debug_build) in diagnostic
        finally:
            saved_config.rename(config)
        fake_bin = root / "fake-bin"
        write(fake_bin / "mpirun", r'''            #!/bin/sh
            printf '%s\0' "$@" > "$MPI_ARGUMENTS"
            if [ "${EXECUTE_FIXTURE-}" = 1 ]; then
                shift 2
                exec "$@"
            fi
        ''')
        (fake_bin / "mpirun").chmod(0o755)
        # Make any runtime Python dependency fail explicitly. The test harness
        # still runs in Python, but input compilation and execution must not.
        write(fake_bin / "python3", '''            #!/bin/sh
            touch "$PYTHON_INVOCATION"
            exit 99
        ''')
        (fake_bin / "python3").chmod(0o755)
        environment["PATH"] = str(fake_bin) + os.pathsep + environment["PATH"]
        environment["MORISBUILDDBG"] = str(debug_build)
        environment["MPI_ARGUMENTS"] = str(root / "mpi-arguments.bin")
        environment["PYTHON_INVOCATION"] = str(root / "python-was-invoked")
        run(SCRIPTS / "moris.sh", "dbg", "3", source.name, "valgrind --quiet", cwd=str(problem), env=environment)
        arguments = (root / "mpi-arguments.bin").read_bytes().decode().split("\0")[:-1]
        assert arguments == ["-np", "3", "valgrind", "--quiet",
                             str(debug_build / "runtime files/moris"), str(source.with_suffix(".so"))]
        # Keep the compile-only mode wrapper and the legacy debug argument working.
        run(SCRIPTS / "moris_create_shared_object.sh", "dbg", source.stem, "dbg", cwd=str(problem), env=environment)
        assert not (root / "python-was-invoked").exists()
        # Independent inputs can compile concurrently against the same MORIS
        # build; neither a shared placeholder nor a global lock is involved.
        jobs = []
        for name, value in [("first", 31), ("second", 37)]:
            extra = problem / (name + ".cpp")
            write(extra, 'extern "C" int input_value() { return ' + str(value) + '; }\n')
            job = subprocess.Popen([str(SCRIPTS / "create_shared_object.sh"), str(problem),
                                    str(debug_build), extra.name], env=environment,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            jobs.append((job, extra, value))
        for job, extra, value in jobs:
            output = job.communicate()[0]
            assert job.returncode == 0, output.decode()
            assert ctypes.CDLL(str(extra.with_suffix(".so"))).input_value() == value

        # A compiler error must propagate and must never launch a stale .so.
        (root / "mpi-arguments.bin").unlink()
        source.write_text("#error intentional compiler failure\n")
        failure = subprocess.run([str(SCRIPTS / "moris.sh"), "dbg", "3", source.name],
                                 cwd=str(problem), env=environment, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert failure.returncode != 0
        assert not (root / "mpi-arguments.bin").exists()
        assert source.with_suffix(".so").is_file()

        # Compile and actually launch both configurations through moris.sh.
        # Fake MPI executes the real fixture binary, so a missing or wrong
        # configuration path cannot pass by merely recording its arguments.
        if shutil.which("ninja"):
            build = root / "multi-config"
            run("cmake", "-S", fixture, "-B", build, "-G", "Ninja Multi-Config",
                "-DMORIS_CMAKE_DIR=" + str(REPOSITORY / "share/cmake"))
            source.write_text('extern "C" int input_value() { return 29; }\n')
            environment["MORISBUILDDBG"] = str(build)
            environment["EXECUTE_FIXTURE"] = "1"
            environment["FIXTURE_EXECUTION"] = str(root / "executed-configuration.txt")
            for configuration in ["Debug", "Release"]:
                run("cmake", "--build", build, "--target", "moris", "--config", configuration)
                environment["MORIS_INPUT_CONFIGURATION"] = configuration
                run(SCRIPTS / "moris.sh", "dbg", "2", source.name,
                    cwd=str(problem), env=environment)
                arguments = (root / "mpi-arguments.bin").read_bytes().decode().split("\0")[:-1]
                assert arguments == ["-np", "2", str(build / "runtime files" / configuration / "moris"),
                                     str(source.with_suffix(".so"))]
                assert (root / "executed-configuration.txt").read_text() == configuration
                loaded = root / ("multi-" + configuration + ".so")
                shutil.copyfile(source.with_suffix(".so"), loaded)
                assert ctypes.CDLL(str(loaded)).input_value() == 29

        assert not (root / "python-was-invoked").exists()

    assert (placeholder.read_bytes(), placeholder.stat().st_mtime_ns, placeholder.is_symlink()) == placeholder_state
    print("External inputs: isolated static/shared builds, Debug IPO overrides, multi-config compilation/execution, paths with spaces, concurrent inputs, shell-only wrapper arguments, and failure propagation passed.")


if __name__ == "__main__":
    main()
