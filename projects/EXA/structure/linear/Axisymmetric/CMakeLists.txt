#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

set(EXAMPLE_MAIN Main_Axisymmetric_Problem)
set(EXAMPLE_FILES Axisymmetric_Problem Axisymmetric_Problem_ConstantThermalField )
set(EXAMPLE_PROCS 1)

# List source files
set(EXAMPLE_SOURCES
    ../../../src/example_main.cpp
    example_test_case.cpp
)

# Include functionality of paths
include_directories( ${CMAKE_BINARY_DIR}/generated/ )

# List include directories
set(EXAMPLE_INCLUDES
    ${COM}
    ${IOS}
    ${WRK}
    ${LINALG}
        )

foreach(EXAMPLE_INCLUDE ${EXAMPLE_INCLUDES})
    include_directories(../${EXAMPLE_INCLUDE}/src)
endforeach()

# List executable dependencies
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
    
# Set the output path for example
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${EXAMPLE_MAIN} ${EXAMPLE_SOURCES})
target_link_libraries(${EXAMPLE_MAIN} PRIVATE ${EXAMPLE_DEPENDENCIES})
set_target_properties(${EXAMPLE_MAIN} PROPERTIES OUTPUT_NAME ${EXAMPLE_MAIN}${EXE_EXT})

# create dynamically linked .so files
set(SO_TPLS
"trilinos"
${ARMADILLO_EIGEN}
)

# get list of INT subfolders
get_property(INT_SRC_LIST GLOBAL PROPERTY INT_SRC_LIST)
get_property(MTK_SRC_LIST GLOBAL PROPERTY MTK_SRC_LIST)

SET(SO_INCLUDES
    ${MORIS_BOOST_INCLUDE_DIRS}
    ${MORIS_TRILINOS_INCLUDE_DIRS}
    ${MORIS_PACKAGE_DIR}/ALG/src
    ${INT_SRC_LIST}
    ${MTK_SRC_LIST}
    ${MORIS_PACKAGE_DIR}/FEM/MSI/src 
    ${MORIS_PACKAGE_DIR}/FEM/VIS/src 
    ${MORIS_PACKAGE_DIR}/GEN/GEN_CORE/src 
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

    ${MORIS_PACKAGE_DIR}/MRS/COR/src)
    
foreach(TPL ${SO_TPLS})
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

# create dynamically linked .so files
foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
    dynamic_link_input(${EXAMPLE_FILE} ${EXAMPLE_FILE} ${EXAMPLE_FILE}.cpp ${SO_INCLUDES})
endforeach()

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    foreach(PROCS ${EXAMPLE_PROCS})
        add_test(NAME Axisymmetric_Problem-${PROCS}-procs
        COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS_EXA} $<TARGET_FILE:${EXAMPLE_MAIN}>
        WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    endforeach()
endif()
