# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

if(NOT MORIS_HAVE_PETSC)
    return()
endif()

set(EXAMPLE_MAIN Main_BDDC_Benchmark)
set(EXAMPLE_FILES BDDC_Plane_Stress BDDC_Elasticity_3D)
set(EXAMPLE_SOURCES example_main.cpp example_test_case.cpp)

# Include functionality of paths
include_directories(${CMAKE_BINARY_DIR}/generated/)

# Include directories.
set(EXAMPLE_INCLUDES ${COM} ${IOS} ${WRK} ${LINALG})

foreach(EXAMPLE_INCLUDE ${EXAMPLE_INCLUDES})
    include_directories(../${EXAMPLE_INCLUDE}/src)
endforeach()

# Executable dependencies.
set(EXAMPLE_DEPENDENCIES
    ${HMR}-lib
    ${COM}-lib
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
)

# Example output directories.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Define the executable and its dependencies.
add_executable(${EXAMPLE_MAIN} ${EXAMPLE_SOURCES})
target_compile_definitions(
    ${EXAMPLE_MAIN} PRIVATE BDDC_INPUT_DIRECTORY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
)
target_link_libraries(${EXAMPLE_MAIN} PRIVATE ${EXAMPLE_DEPENDENCIES})
set_target_properties(${EXAMPLE_MAIN} PROPERTIES OUTPUT_NAME ${EXAMPLE_MAIN}${EXE_EXT})

set(SO_TPLS "trilinos" ${ARMADILLO_EIGEN})

# get list of INT subfolders
get_property(INT_SRC_LIST GLOBAL PROPERTY INT_SRC_LIST)
get_property(MTK_SRC_LIST GLOBAL PROPERTY MTK_SRC_LIST)

set(SO_INCLUDES
    ${MORIS_BOOST_INCLUDE_DIRS}
    ${MORIS_TRILINOS_INCLUDE_DIRS}
    ${MORIS_PACKAGE_DIR}/ALG/src
    ${INT_SRC_LIST}
    ${MTK_SRC_LIST}
    ${MORIS_PACKAGE_DIR}/FEM/MSI/src
    ${MORIS_PACKAGE_DIR}/FEM/VIS/src
    ${MORIS_PACKAGE_DIR}/SOL/DLA/src
    ${MORIS_PACKAGE_DIR}/SOL/TSA/src
    ${MORIS_PACKAGE_DIR}/SOL/NLA/src
    ${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src
    ${MORIS_PACKAGE_DIR}/LINALG/src
    ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
    ${MORIS_PACKAGE_DIR}/COM/src
    ${MORIS_PACKAGE_DIR}/PRM/ENM/src
    ${MORIS_PACKAGE_DIR}/MTK/src
    ${MORIS_PACKAGE_DIR}/HMR/src
    ${MORIS_PACKAGE_DIR}/XTK/src
    ${MORIS_PACKAGE_DIR}/PRM/src
    ${MORIS_PACKAGE_DIR}/MRS/COR/src
)

foreach(TPL ${SO_TPLS})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

# create dynamically linked .so files
foreach(EXAMPLE_FILE IN LISTS EXAMPLE_FILES)
    dynamic_link_input(${EXAMPLE_FILE} ${EXAMPLE_FILE} ${EXAMPLE_FILE}.cpp ${SO_INCLUDES})
endforeach()

# The direct AIJ reference and parallel BDDC coarse solve use PETSc's MUMPS.
# Check PETSc itself, independently of MORIS_USE_MUMPS for other libraries.
include(CMakePushCheckState)
include(CheckCXXSourceCompiles)
cmake_push_check_state(RESET)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIR})
if(PETSC_ARCH_INCLUDE_DIR)
    list(APPEND CMAKE_REQUIRED_INCLUDES ${PETSC_ARCH_INCLUDE_DIR})
endif()
check_cxx_source_compiles(
    "#include <petscconf.h>
     #ifndef PETSC_HAVE_MUMPS
     #error PETSc was built without MUMPS
     #endif
     int main() { return 0; }"
    MORIS_BDDC_BENCHMARK_PETSC_HAS_MUMPS
)
cmake_pop_check_state()

foreach(dim 2 3)
    if(dim EQUAL 2)
        set(example_case BDDC_Plane_Stress)
        set(parallel_procs 4)
    else()
        set(example_case BDDC_Elasticity_3D)
        set(parallel_procs 8)
    endif()
    set(process_counts 1)
    if(MORIS_HAVE_PARALLEL_TESTS AND MORIS_BDDC_BENCHMARK_PETSC_HAS_MUMPS)
        list(APPEND process_counts ${parallel_procs})
    endif()

    foreach(procs IN LISTS process_counts)
        set(variants homogeneous jump)
        if(procs GREATER 1)
            list(APPEND variants corners aij)
        endif()
        foreach(variant IN LISTS variants)
            set(test_name EXA-structure-linear-bddc-${dim}D-${variant}-${procs}-procs)
            set(test_directory ${CMAKE_CURRENT_BINARY_DIR}/cases/${test_name})
            file(MAKE_DIRECTORY ${test_directory})
            set(test_case ${example_case})
            if(variant STREQUAL "jump")
                string(APPEND test_case "_Jump")
            elseif(variant STREQUAL "corners")
                string(APPEND test_case "_Corners")
            elseif(variant STREQUAL "aij")
                string(APPEND test_case "_AIJ")
            endif()

            # Follow the other EXA tests: Catch cases, dynamically loaded
            # input libraries, and one complete HMR/GEN/XTK/FEM/SOL/VIS workflow.
            add_test(NAME ${test_name}
                     COMMAND ${MORIS_EXECUTE_COMMAND} -n ${procs} ${VALGRIND}
                             ${VALGRIND_OPTIONS_EXA} $<TARGET_FILE:${EXAMPLE_MAIN}> ${test_case}
            )
            set_tests_properties(
                ${test_name}
                PROPERTIES PROCESSORS
                           ${procs}
                           TIMEOUT
                           600
                           LABELS
                           "EXA;structure;linear;bddc"
                           WORKING_DIRECTORY
                           ${test_directory}
            )
        endforeach()
    endforeach()
endforeach()
