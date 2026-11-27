# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

set(EXAMPLE_MAIN Main_Coulomb_Shear)
set(EXAMPLE_FILES Coulomb_Shear)

# Sources.
set(EXAMPLE_SOURCES ../../../src/example_main.cpp example_test_case.cpp)

# Include functionality of paths
include_directories(${CMAKE_BINARY_DIR}/generated/)

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

# The input and shared EXA main use driver globals for fixed regression settings.
target_compile_definitions(Coulomb_Shear PRIVATE MORIS_COULOMB_SHEAR_DRIVER)
target_compile_definitions(
    ${EXAMPLE_MAIN} PRIVATE MORIS_COULOMB_SHEAR_LIBRARY="$<TARGET_FILE:Coulomb_Shear>"
)
add_dependencies(${EXAMPLE_MAIN} Coulomb_Shear)

foreach(CASE Cycle SubdividedLoad RefinedMesh)
    set(CASE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/runs/$<CONFIG>/${CASE}-1")
    file(GENERATE OUTPUT "${CASE_DIRECTORY}/.moris-test-directory" CONTENT "")
    add_test(NAME Coulomb_Shear-${CASE}-1-procs WORKING_DIRECTORY "${CASE_DIRECTORY}"
             COMMAND ${MORIS_EXECUTE_COMMAND} -n 1 $<TARGET_FILE:${EXAMPLE_MAIN}>
                     "[CoulombShear${CASE}]"
    )
endforeach()
