# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.

# Keep the supported minimum while adopting policies tested with CMake 4.4.
cmake_minimum_required(VERSION 3.17...4.4)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not supported. Use cmake -S . -B <build-directory>.")
endif()

# Compiler selection must precede project(). Standard CMake toolchains, CC/CXX,
# and CMAKE_<LANG>_COMPILER remain authoritative; accept the legacy aliases too.
foreach(language IN ITEMS C CXX)
    if(MORIS_${language}_COMPILER AND NOT DEFINED CMAKE_${language}_COMPILER)
        set(CMAKE_${language}_COMPILER "${MORIS_${language}_COMPILER}" CACHE FILEPATH
                                                                             "${language} compiler"
        )
    endif()
endforeach()
project(MORIS VERSION 1.0.0 LANGUAGES C CXX)

set(MORIS_VERSION_RELEASE 0)
if(NOT MORIS_VERSION_RELEASE)
    set(MORIS_VERSION "${MORIS_VERSION}.dbg")
endif()

include("${PROJECT_SOURCE_DIR}/share/cmake/moris_paths.cmake")
include("${MORIS_CMAKE_DIR}/moris_options.cmake")
include("${MORIS_CMAKE_DIR}/utilities/moris_macros.cmake")
include("${MORIS_CMAKE_DIR}/utilities/moris_package.cmake")
include("${MORIS_CMAKE_DIR}/utilities/dynamic_link_input.cmake")
include("${MORIS_CMAKE_DIR}/moris_features.cmake")
include("${MORIS_CMAKE_DIR}/moris_compiler.cmake")
include("${MORIS_CMAKE_DIR}/moris_runtime.cmake")

# Directory requirements cover legacy executables and plugins as well as modules.
# Public module usage requirements continue to be exported by their targets.
list(TRANSFORM MORIS_DEFINITIONS REPLACE "^-D" "")
list(REMOVE_DUPLICATES MORIS_DEFINITIONS)
set(MORIS_COMPILE_DEFINITIONS "${MORIS_DEFINITIONS}")
# Legacy config variables contain command-line -D flags; target/directory APIs
# use the unprefixed copy so older CMake versions handle them consistently.
list(TRANSFORM MORIS_DEFINITIONS PREPEND "-D")
add_compile_definitions(${MORIS_COMPILE_DEFINITIONS})
include_directories(${MORIS_INCDIRS})
set(EXECUTABLE "${MORIS}.exe")
set(EXE_EXT ".exe")

include("${MORIS_CMAKE_DIR}/moris_packages.cmake")
enable_testing()

# Shared mode assembles module objects only after all component usage
# requirements are known. Static mode preserves the component archives.
include("${MORIS_CMAKE_DIR}/utilities/moris_shared.cmake")
moris_initialize_shared_library()
add_subdirectory(projects)
moris_finalize_shared_library()

# External input files reuse this build's requirements without replacing a source
# file or rebuilding a shared input target in the MORIS tree.
if(TARGET shared_object_file)
    include("${MORIS_CMAKE_DIR}/utilities/moris_input_build.cmake")
    moris_export_input_build(shared_object_file)
endif()

add_subdirectory(share/cmake/find_modules)
add_subdirectory(share/cmake/find_wrappers)
add_subdirectory(share/cmake/utilities)
add_subdirectory(share/doc)
include("${MORIS_CMAKE_DIR}/moris_install.cmake")
