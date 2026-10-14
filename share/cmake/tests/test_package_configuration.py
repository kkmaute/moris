#!/usr/bin/env python3
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
"""Build and consume relocated component packages in static/shared and OpenMP modes."""
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap


def run(*args):
    subprocess.run(args, check=True)


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text))


def main():
    cmake_dir = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="moris-package-test-") as temporary:
        root = Path(temporary)
        # CI uses a small API-compatible stand-in; a local prefix can exercise
        # the real profiler without making the basic regression depend on it.
        profiling_prefix = os.environ.get("MORIS_TEST_GPERFTOOLS_PREFIX")
        if not profiling_prefix:
            mock = root / "profiler"
            write(mock / "CMakeLists.txt", '''
                cmake_minimum_required(VERSION 3.17)
                project(Profiler LANGUAGES CXX)
                add_library(tcmalloc_and_profiler SHARED profiler.cpp)
                install(TARGETS tcmalloc_and_profiler LIBRARY DESTINATION lib)
                install(FILES profiler.h DESTINATION include/gperftools)
            ''')
            write(mock / "profiler.h", 'extern "C" int ProfilerStart(const char*);\nextern "C" void ProfilerStop();\n')
            write(mock / "profiler.cpp", 'extern "C" int ProfilerStart(const char*) { return 1; }\nextern "C" void ProfilerStop() {}\n')
            profiling_prefix = str(root / "profiler-install")
            run("cmake", "-S", str(mock), "-B", str(root / "profiler-build"),
                "-DCMAKE_INSTALL_PREFIX=" + profiling_prefix)
            run("cmake", "--build", str(root / "profiler-build"))
            run("cmake", "--install", str(root / "profiler-build"))
        source = root / "source"
        write(source / "CMakeLists.txt", '''
            cmake_minimum_required(VERSION 3.17)
            project(Packages VERSION 1.0.0 LANGUAGES CXX)
            set(CMAKE_CXX_STANDARD 17)
            set(MORIS moris)
            set(MORIS_VERSION 1.0.0)
            set(MORIS_VERSION_MAJOR 1)
            set(CONFIG_BUILD_DIR "${CMAKE_BINARY_DIR}/cmake")
            set(CONFIG_INSTALL_DIR share)
            set(LIB_INSTALL_DIR lib)
            set(EXE_INSTALL_DIR bin)
            set(MORIS_SOURCE_DIRS XTK CNT)
            include("${MORIS_CMAKE_DIR}/utilities/moris_shared.cmake")
            include("${MORIS_CMAKE_DIR}/utilities/moris_package.cmake")
            if(MORIS_USE_OPENMP)
                find_package(OpenMP REQUIRED COMPONENTS CXX)
            endif()
            set(MORIS_TPL_DIR "${MORIS_CMAKE_DIR}/find_wrappers")
            list(APPEND CMAKE_MODULE_PATH "${MORIS_CMAKE_DIR}/find_modules")
            set(XTK_TPL_DEPENDENCIES gperftools)
            set(CNT_TPL_DEPENDENCIES gperftools)
            include("${MORIS_TPL_DIR}/gperftools_new.cmake")
            moris_initialize_shared_library()
            foreach(component IN ITEMS XTK CNT)
                set(${component}_CONFIG_BUILD_DIR "${CONFIG_BUILD_DIR}/${component}")
                set(${component}_CONFIG_INSTALL_DIR "share/${component}")
                moris_add_module(${component}-lib ${component}_TARGET ${component}.cpp)
                target_link_libraries(${${component}_TARGET} PUBLIC moris::gperftools)
                if(MORIS_USE_OPENMP)
                    target_link_libraries(${${component}_TARGET} PUBLIC OpenMP::OpenMP_CXX)
                endif()
                install(TARGETS ${component}-lib EXPORT ${component}Targets ARCHIVE DESTINATION lib)
                install(EXPORT ${component}Targets DESTINATION share/${component})
            endforeach()
            target_link_libraries(${XTK_TARGET} PUBLIC CNT-lib)
            target_link_libraries(${CNT_TARGET} PUBLIC XTK-lib)
            moris_install_component_package(XTK TARGET ${XTK_TARGET} DEPENDENCIES CNT-lib TPLS gperftools)
            moris_install_component_package(CNT TARGET ${CNT_TARGET} DEPENDENCIES XTK-lib TPLS gperftools)
            moris_finalize_shared_library()
            get_property(MORIS_INSTALLED_COMPONENTS GLOBAL PROPERTY MORIS_INSTALLED_COMPONENTS)
            get_property(MORIS_INSTALLED_COMPONENT_DIRS GLOBAL PROPERTY MORIS_INSTALLED_COMPONENT_DIRS)
            configure_package_config_file("${MORIS_CMAKE_DIR}/morisConfig.cmake.in"
                "${CONFIG_BUILD_DIR}/moris/morisConfig.cmake"
                INSTALL_DESTINATION share/moris PATH_VARS CONFIG_INSTALL_DIR)
            install(FILES "${CONFIG_BUILD_DIR}/moris/morisConfig.cmake" DESTINATION share/moris)
            add_library(moris_general INTERFACE)
            install(TARGETS moris_general EXPORT morisTargets)
            install(EXPORT morisTargets DESTINATION share/moris)
            install(FILES "${MORIS_CMAKE_DIR}/find_modules/FindGperftools.cmake" DESTINATION share/find_modules)
            install(FILES "${MORIS_CMAKE_DIR}/find_wrappers/gperftools_new.cmake" DESTINATION share/find_wrappers)
            install(FILES "${MORIS_CMAKE_DIR}/utilities/moris_macros.cmake" DESTINATION share/utilities)
        ''')
        write(source / "XTK.cpp", 'int cnt();\nint xtk() { return cnt() + 2; }\n')
        write(source / "CNT.cpp", 'int cnt() { return 40; }\n')
        consumer = root / "consumer"
        write(consumer / "CMakeLists.txt", '''
            cmake_minimum_required(VERSION 3.17)
            project(Consumer LANGUAGES CXX)
            if(CHECK_MISSING)
                find_package(moris REQUIRED COMPONENTS Missing)
                message(FATAL_ERROR "Accepted missing required component")
            endif()
            # A standalone XTK import must discover OpenMP and libmoris itself.
            find_package(XTK REQUIRED CONFIG)
            if(EXPECT_OPENMP AND NOT TARGET OpenMP::OpenMP_CXX)
                message(FATAL_ERROR "Missing OpenMP dependency")
            endif()
            find_package(moris REQUIRED COMPONENTS CNT)
            find_package(moris REQUIRED COMPONENTS XTK)
            find_package(moris REQUIRED OPTIONAL_COMPONENTS Missing)
            if(moris_Missing_FOUND)
                message(FATAL_ERROR "Accepted missing optional component")
            endif()
            find_package(moris REQUIRED)
            if(NOT TARGET moris::CNT OR NOT TARGET moris::XTK)
                message(FATAL_ERROR "Repeated package lookup lost components")
            endif()
            add_executable(check main.cpp)
            target_link_libraries(check PRIVATE moris::all_libs)
        ''')
        write(consumer / "main.cpp", '#include <gperftools/profiler.h>\nint xtk();\nint main() { if (!ProfilerStart("consumer.prof")) return 2; ProfilerStop(); return xtk() == 42 ? 0 : 1; }\n')
        for shared in ("OFF", "ON"):
            for openmp in ("OFF", "ON"):
                case = root / ("shared-" + shared + "-openmp-" + openmp)
                build, prefix = case / "build", case / "install"
                run("cmake", "-S", str(source), "-B", str(build),
                    "-DMORIS_CMAKE_DIR=" + str(cmake_dir),
                    "-DCMAKE_INSTALL_PREFIX=" + str(prefix),
                    "-DCMAKE_PREFIX_PATH=" + profiling_prefix,
                    "-DMORIS_HAVE_SHARED=" + shared, "-DMORIS_USE_OPENMP=" + openmp)
                run("cmake", "--build", str(build), "--parallel", "2")
                run("cmake", "--install", str(build))
                relocated = case / "relocated"
                prefix.rename(relocated)
                consumer_build = case / "consumer"
                run("cmake", "-S", str(consumer), "-B", str(consumer_build),
                    "-DCMAKE_PREFIX_PATH=" + str(relocated) + ";" + profiling_prefix, "-DEXPECT_OPENMP=" + openmp)
                run("cmake", "--build", str(consumer_build), "--parallel", "2")
                subprocess.run([str(consumer_build / "check")], cwd=str(consumer_build), check=True)
                result = subprocess.run(
                    ["cmake", "-S", str(consumer), "-B", str(case / "missing"),
                     "-DCMAKE_PREFIX_PATH=" + str(relocated) + ";" + profiling_prefix, "-DCHECK_MISSING=ON"],
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
                assert result.returncode and "Accepted missing" not in result.stdout, result.stdout
    print("Static/shared packages, profiling, OpenMP, relocation, cycles, and component lookup checks passed.")


if __name__ == "__main__":
    main()
