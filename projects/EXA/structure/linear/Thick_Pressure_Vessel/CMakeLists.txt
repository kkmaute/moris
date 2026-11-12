# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

set(EXAMPLE_MAIN Main_Thick_Pressure_Vessel)
set(EXAMPLE_FILES Pressure_Vessel_3D Pressure_Vessel_3D_Immersed)
set(EXAMPLE_PROCS 1 4)

# Sources.
set(EXAMPLE_SOURCES ../../../src/example_main.cpp example_test_case.cpp)

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
target_link_libraries(${EXAMPLE_MAIN} PRIVATE ${EXAMPLE_DEPENDENCIES})
# Inputs stay in the build output directory while each test writes to its own case directory.
target_compile_definitions(
    ${EXAMPLE_MAIN} PRIVATE PRESSURE_VESSEL_INPUT_DIRECTORY="$<TARGET_FILE_DIR:Pressure_Vessel_3D>"
)
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
foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
    dynamic_link_input(${EXAMPLE_FILE} ${EXAMPLE_FILE} ${EXAMPLE_FILE}.cpp ${SO_INCLUDES})
endforeach()

# PETSc cases are launched separately inside each process-count group.
# Check capabilities before adding BDDC variants to the existing vessel tests.
if(MORIS_HAVE_PETSC)
    # Cut-cell subdomains need a pivoting factorization. Check PETSc's actual
    # capabilities rather than MORIS_USE_MUMPS, which also controls other TPLs.
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
        MORIS_PRESSURE_VESSEL_PETSC_HAS_MUMPS
    )
    cmake_pop_check_state()

    if(MORIS_HAVE_PARALLEL_TESTS AND NOT MORIS_PRESSURE_VESSEL_PETSC_HAS_MUMPS)
        message(
            STATUS
                "Parallel pressure-vessel BDDC cases require PETSc with MUMPS; retaining serial BDDC coverage."
        )
    endif()
endif()

set(process_counts 1)
if(MORIS_HAVE_PARALLEL_TESTS)
    set(process_counts ${EXAMPLE_PROCS})
endif()
foreach(procs IN LISTS process_counts)
    set(test_name Thick_Pressure_Vessel-${procs}-procs)
    set(group_cases "")
    set(group_case_count 1)
    set(test_directory ${CMAKE_CURRENT_BINARY_DIR}/cases/${test_name}/amesos)
    file(MAKE_DIRECTORY "${test_directory}")
    moris_append_test_group_case(
        group_cases
        "${test_name}-amesos"
        "${test_directory}"
        600
        ${MORIS_EXECUTE_COMMAND}
        -n
        ${procs}
        ${VALGRIND}
        ${VALGRIND_OPTIONS_EXA}
        $<TARGET_FILE:${EXAMPLE_MAIN}>
    )

    if(MORIS_HAVE_PETSC AND (procs EQUAL 1 OR MORIS_PRESSURE_VESSEL_PETSC_HAS_MUMPS))
        foreach(example IN LISTS EXAMPLE_FILES)
            set(case_name ${example}-BDDC-${procs}-procs)
            set(case_directory ${CMAKE_CURRENT_BINARY_DIR}/cases/${test_name}/${example}-BDDC)
            file(MAKE_DIRECTORY "${case_directory}")
            moris_append_test_group_case(
                group_cases
                "${case_name}"
                "${case_directory}"
                600
                ${MORIS_EXECUTE_COMMAND}
                -n
                ${procs}
                ${VALGRIND}
                ${VALGRIND_OPTIONS_EXA}
                $<TARGET_FILE:${EXAMPLE_MAIN}>
                "${example}_BDDC"
            )
            math(EXPR group_case_count "${group_case_count} + 1")
        endforeach()
    endif()

    math(EXPR group_timeout "600 * ${group_case_count} + 30")
    moris_register_test_group(
        ${test_name}
        "${group_cases}"
        PROCESSORS
        ${procs}
        TIMEOUT
        ${group_timeout}
        LABELS
        "EXA;structure;linear;pressure_vessel"
    )
endforeach()
