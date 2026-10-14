#!/usr/bin/env python3
"""Exercise shared MORIS target assembly without third-party dependencies.

Run with: python3 share/cmake/tests/test_shared_library.py
"""
from pathlib import Path
import subprocess
import tempfile
import textwrap


def run(*args):
    subprocess.run(args, check=True)


def write(path, contents):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(contents))


def main():
    cmake_dir = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="moris-shared-test-") as temporary:
        root = Path(temporary)
        source, build, prefix = root / "source", root / "build", root / "install"
        write(source / "CMakeLists.txt", r'''
            cmake_minimum_required(VERSION 3.17)
            project(SharedRegression VERSION 1.0.0 LANGUAGES CXX)
            set(MORIS moris)
            option(MORIS_HAVE_SHARED "Build the monolithic shared library" ON)
            set(MORIS_VERSION 1.0.0)
            set(MORIS_VERSION_MAJOR 1)
            set(CONFIG_BUILD_DIR "${CMAKE_BINARY_DIR}/cmake")
            set(CONFIG_INSTALL_DIR share)
            set(HEADER_INSTALL_DIR include)
            set(LIB_INSTALL_DIR lib)
            set(EXE_INSTALL_DIR bin)
            set(MORIS_SOURCE_DIRS A B)
            set(A_TPL_DEPENDENCIES mock)
            set(MORIS_TPL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tpls")
            file(WRITE "${CMAKE_BINARY_DIR}/generated/paths.hpp" "// generated\n")
            include("${MORIS_CMAKE_DIR}/utilities/moris_shared.cmake")
            moris_initialize_shared_library()
            if(NOT MORIS_HAVE_SHARED)
                include("${MORIS_TPL_DIR}/mock_new.cmake")
            endif()
            add_library(LINALG-lib INTERFACE)
            target_include_directories(LINALG-lib INTERFACE
                "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/headers>"
                "$<INSTALL_INTERFACE:include>")
            target_compile_definitions(LINALG-lib INTERFACE HEADER_VALUE=MOCK_VALUE)
            add_subdirectory(a)
            add_subdirectory(b)
            target_link_libraries(LINALG-lib INTERFACE A-lib)
            moris_finalize_shared_library()
            install(TARGETS A-lib B-lib LINALG-lib EXPORT Components)
            install(EXPORT Components DESTINATION share/moris)
            install(FILES headers/value.hpp DESTINATION include)
            install(FILES tpls/mock_new.cmake DESTINATION share/find_wrappers)
            install(FILES "${MORIS_CMAKE_DIR}/utilities/moris_macros.cmake"
                DESTINATION share/utilities)
            add_executable(check main.cpp)
            target_link_libraries(check PRIVATE A-lib)
        ''')
        write(source / "a/CMakeLists.txt", '''
            moris_add_module(A-lib A_TARGET a.cpp)
            target_link_libraries(${A_TARGET} PUBLIC B-lib LINALG-lib moris::mock)
        ''')
        write(source / "b/CMakeLists.txt", '''
            moris_add_module(B-lib B_TARGET b.cpp)
            target_link_libraries(${B_TARGET} PUBLIC A-lib)
        ''')
        write(source / "tpls/mock_new.cmake", '''
            if(NOT TARGET moris::mock)
                add_library(moris::mock INTERFACE IMPORTED)
                set_property(TARGET moris::mock PROPERTY INTERFACE_COMPILE_DEFINITIONS MOCK_VALUE=40)
            endif()
        ''')
        write(source / "headers/value.hpp", '''
            #pragma once
            #ifndef HEADER_VALUE
            #error Missing propagated compile definition
            #endif
            int a();
            int b();
        ''')
        write(source / "a/a.cpp", '#include "value.hpp"\nint a() { return b() + 2; }\n')
        write(source / "b/b.cpp", '#include "value.hpp"\nint b() { return HEADER_VALUE; }\n')
        write(source / "main.cpp", '#include "value.hpp"\nint main() { return a() == 42 ? 0 : 1; }\n')
        run("cmake", "-S", str(source), "-B", str(build),
            f"-DMORIS_CMAKE_DIR={cmake_dir}", f"-DCMAKE_INSTALL_PREFIX={prefix}")
        run("cmake", "--build", str(build), "--parallel", "2")
        run(str(build / "check"))
        run("cmake", "--install", str(build))
        assert not list(prefix.rglob("*.a")), "Shared mode installed a static archive"
        assert not list(prefix.rglob("*.o")), "Internal objects leaked into the installation"
        # Load exports from a relocated installation, not from the build tree.
        moved = root / "relocated"
        prefix.rename(moved)
        consumer = root / "consumer"
        write(consumer / "CMakeLists.txt", '''
            cmake_minimum_required(VERSION 3.17)
            project(Consumer LANGUAGES CXX)
            set(PACKAGE_PREFIX_DIR "${MORIS_TEST_PREFIX}")
            include("${MORIS_TEST_PREFIX}/share/moris/morisSharedConfig.cmake")
            include("${MORIS_TEST_PREFIX}/share/moris/Components.cmake")
            add_executable(consumer main.cpp)
            target_link_libraries(consumer PRIVATE B-lib)
        ''')
        write(consumer / "main.cpp", (source / "main.cpp").read_text())
        run("cmake", "-S", str(consumer), "-B", str(root / "consumer-build"), f"-DMORIS_TEST_PREFIX={moved}")
        run("cmake", "--build", str(root / "consumer-build"))
        run(str(root / "consumer-build/consumer"))
        # The same cyclic component graph must retain the original static path.
        static_build = root / "static-build"
        run("cmake", "-S", str(source), "-B", str(static_build),
            f"-DMORIS_CMAKE_DIR={cmake_dir}", "-DMORIS_HAVE_SHARED=OFF")
        run("cmake", "--build", str(static_build), "--parallel", "2")
        run(str(static_build / "check"))
        assert list(static_build.rglob("*.a")), "Static mode did not build archives"
    print("Shared-library cycle, runtime, relocated export, and static-mode checks passed.")


if __name__ == "__main__":
    main()
