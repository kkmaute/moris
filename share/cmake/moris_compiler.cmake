# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Project compilation policy. Do not replace CMAKE_<LANG>_FLAGS: callers and
# toolchain files may supply ABI, architecture, or sanitizer settings there.
include(CheckCXXCompilerFlag)
include(CheckIPOSupported)

set(MORIS_C_COMPILER "${CMAKE_C_COMPILER}")
set(MORIS_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
set(MORIS_C_LINK_EXECUTABLE "${CMAKE_C_COMPILER}")
set(MORIS_CXX_LINK_EXECUTABLE "${CMAKE_CXX_COMPILER}")

# Public headers use C++17. Falling back to C++11/14 only postpones failure.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
# Static archives are also linked into input-file plugins.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Symbolic modes keep their historical Debug selection. Otherwise respect an
# explicit build type and default to Release only for single-config generators.
if(MORIS_HAVE_SYMBOLIC)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build configuration" FORCE)
elseif(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration" FORCE)
endif()
if(DEFINED CACHE{CMAKE_BUILD_TYPE})
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# Each probe needs its own cache entry: reusing HAVE_PEDANTIC for several flags
# incorrectly treats the first successful check as proof that every flag works.
macro(moris_append_cxx_flag flag)
    string(MAKE_C_IDENTIFIER "MORIS_SUPPORTS_${flag}" flag_result)
    check_cxx_compiler_flag("${flag}" ${flag_result})
    if(${flag_result})
        list(APPEND moris_compile_options "${flag}")
    endif()
endmacro()

# Retain the legacy user-supplied MORIS_CXX_FLAGS alongside standard CMake flags.
separate_arguments(moris_compile_options NATIVE_COMMAND "${MORIS_CXX_FLAGS}")
set(moris_link_options ${moris_compile_options})
if(MORIS_HAVE_PEDANTIC)
    foreach(flag IN ITEMS -pedantic -pedantic-errors -Wall -Werror -Wno-unused-parameter)
        moris_append_cxx_flag("${flag}")
    endforeach()
endif()
moris_append_cxx_flag(-m64)
moris_append_cxx_flag(-fextended-identifiers)

if(MORIS_HAVE_SYMBOLIC)
    if(NOT MORIS_HAVE_SYMBOLIC_STRONG)
        moris_append_cxx_flag(-Og)
    endif()
    moris_append_cxx_flag(-g3)
    moris_append_cxx_flag(-ggdb)
elseif(MORIS_USE_INTEL)
    foreach(flag IN ITEMS -xCORE-AVX2 -mkl:sequential -ip)
        moris_append_cxx_flag("${flag}")
        list(APPEND moris_link_options "${flag}")
    endforeach()
endif()

# Let CMake choose the compiler's IPO flags and matching archive tools.
# Only optimized configurations use IPO; unsupported toolchains remain usable.
if(MORIS_USE_LTO AND NOT MORIS_HAVE_SYMBOLIC)
    check_ipo_supported(RESULT moris_ipo_supported OUTPUT moris_ipo_error LANGUAGES CXX)
    if(moris_ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL ON)
    else()
        message(WARNING "MORIS_USE_LTO requested but unavailable: ${moris_ipo_error}")
    endif()
endif()

# Plugins resolve executable-owned symbols (notably the global logger).
# This is a linker requirement, not a C++ compilation flag.
check_cxx_compiler_flag(-rdynamic MORIS_SUPPORTS_RDYNAMIC)
if(MORIS_SUPPORTS_RDYNAMIC)
    list(APPEND moris_link_options -rdynamic)
endif()

if(MORIS_USE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    link_libraries(OpenMP::OpenMP_CXX)
endif()
find_package(Threads REQUIRED)

if(MORIS_USE_CHECK_MEMORY)
    # Only executables own the allocation hooks; do not wrap the shared library.
    foreach(symbol IN ITEMS malloc free posix_memalign)
        add_link_options(
            "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:LINKER:--wrap=${symbol}>"
        )
    endforeach()
endif()

if(MORIS_USE_GPERFTOOLS)
    include("${MORIS_TPL_DIR}/gperftools_new.cmake")

endif()

if(MORIS_USE_MEMORY_DEBUG)
    find_library(MORIS_LIB_MemKind NAMES memkind libmemkind.so.0 HINTS /usr/lib64)
    if(MORIS_LIB_MemKind)
        link_libraries("${MORIS_LIB_MemKind}")
    else()
        message(WARNING "memkind was not found; memory debugging is unavailable.")
    endif()
endif()

if(MORIS_USE_NTOP_LIB)
    set(MORIS_NTOP_LIBRARY "/home/maute/codes/ntop/lib/libntop_core.so.4"
        CACHE FILEPATH "NTOP interface library"
    )
    link_libraries("${MORIS_NTOP_LIBRARY}")
endif()

# USE_SACADO is retained as a compatibility alias for older build scripts.
if(MORIS_USE_SACADO OR USE_SACADO)
    set(MORIS_SACADO_INCLUDE_DIR "$ENV{HOME}/codes/Sacado/Trilinos-master/packages/sacado/src"
        CACHE PATH "Sacado source include directory"
    )
    include_directories(
        "${MORIS_SACADO_INCLUDE_DIR}" "${MORIS_SACADO_INCLUDE_DIR}/mpl"
        "${MORIS_SACADO_INCLUDE_DIR}/parameter"
    )
endif()

# Preserve the narrowly scoped GCC 10.2 workarounds still used on older systems.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION MATCHES "^10\\.2(\\.|$)")
    foreach(flag IN ITEMS -Wno-type-limits -Wno-deprecated-declarations -Wno-maybe-uninitialized)
        moris_append_cxx_flag("${flag}")
    endforeach()
endif()

if(MORIS_USE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES clang-tidy REQUIRED)
    foreach(flag IN
            ITEMS -Wno-unused-command-line-argument -Wno-keyword-macro -Wno-return-type-c-linkage
                  -Wno-deprecated-declarations -Wpessimizing-move
    )
        moris_append_cxx_flag("${flag}")
    endforeach()
    if(DEFINED ENV{MPI_HOME})
        include_directories("$ENV{MPI_HOME}/include")
    endif()
    set(CLANG_TIDY_CHECKS
        "-*"
        readability-redundant-casting
        readability-redundant-control-flow
        readability-duplicate-include
        readability-container-contains
        readability-const-return-type
        modernize-concat-nested-namespaces
        modernize-deprecated-headers
        modernize-use-nullptr
        modernize-use-override
        performance-unnecessary-copy-initialization
        performance-unnecessary-value-param
        performance-move-const-arg
        performance-inefficient-algorithm
        performance-trivially-destructible
        performance-no-automatic-move
        performance-inefficient-vector-operation
        performance-for-range-copy
        performance-faster-string-find
        performance-avoid-endl
    )
    list(JOIN CLANG_TIDY_CHECKS "," CLANG_TIDY_CHECKS)
    # Keep the established fix-in-place behavior and exclude third-party headers.
    set(CMAKE_CXX_CLANG_TIDY
        "${CLANG_TIDY_EXE}" --format-style=file "-header-filter=.+projects.+"
        "--exclude-header-filter=.+spack.+|cl_Parameter\\.hpp" "-checks=${CLANG_TIDY_CHECKS}" -fix
    )
endif()

list(REMOVE_DUPLICATES moris_compile_options)
foreach(flag IN LISTS moris_compile_options)
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${flag}>")
endforeach()
add_link_options(${moris_link_options})
# Legacy package config variables remain available to non-target consumers.
set(moris_legacy_flags ${moris_compile_options} ${moris_link_options}
                       ${CMAKE_CXX17_STANDARD_COMPILE_OPTION} ${CMAKE_CXX_COMPILE_OPTIONS_PIC}
)
list(REMOVE_DUPLICATES moris_legacy_flags)
list(JOIN moris_legacy_flags " " MORIS_CXX_FLAGS)
