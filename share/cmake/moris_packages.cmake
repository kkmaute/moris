# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

include("${MORIS_CMAKE_DIR}/utilities/moris_testing.cmake")
include("${MORIS_CMAKE_DIR}/utilities/moris_select_modules.cmake")
moris_select_modules()
include(${MORIS_DEPENDS_DIR}/MRS_includes.cmake)

if(BUILD_MAIN)
    include(${MORIS_DEPENDS_DIR}/main_includes.cmake)
endif()

# Dependency files establish module ordering as well as transitive dependencies.
foreach(component IN LISTS MORIS_BUILD_MODULES)
    if(BUILD_${component}
       AND NOT component IN_LIST MORIS_CORE_MODULES
       AND NOT component STREQUAL "MOD"
       AND NOT component STREQUAL "STK"
    )
        include("${MORIS_DEPENDS_DIR}/${component}_Depends.cmake")
    endif()
endforeach()
# MORIS valgrind
if(MORIS_USE_VALGRIND)
    set(VALGRIND valgrind)
    set(VALGRIND_OPTIONS --suppressions=${MORIS_SOURCE_DIR}/share/valgrind/moris.supp)
    set(VALGRIND_OPTIONS_EXA --leak-check=full --track-origins=yes
                             --suppressions=${MORIS_SOURCE_DIR}/share/valgrind/moris.supp
    )
else()
    set(VALGRIND "")
    set(VALGRIND_OPTIONS "")
    set(VALGRIND_OPTIONS_EXA "")
endif()

# MORIS unit tests
if(MORIS_USE_TESTS)
    include(${MORIS_DEPENDS_DIR}/test_includes.cmake)
endif()

if(MORIS_HAVE_PARALLEL_TESTS)
    set(MORIS_TEST_PROCS 1 2 4)
    set(MORIS_TEST_PROCS_HMR 1 2 4 12)
endif()

# MORIS examples
if(MORIS_USE_EXAMPLES)
    message(STATUS "Build examples")
    include(${MORIS_DEPENDS_DIR}/EXA_Depends.cmake)
endif()
