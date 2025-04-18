#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Parameters Tests ----------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES 
    test_main.cpp
    UT_WRK_Workflow.cpp
    UT_WRK_Remeshing_Mini_Performer.cpp )

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${WRK}-lib
    ${XTK}-lib
    ${HMR}-lib
    )

set(GEN_INCLUDES
	"ADV"
	"ADV/fields"
	"PDV"
    )

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()
include_directories(${MORIS_PACKAGE_DIR}/FEM/INT/src/IP)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(WRK_TEST ${WRK}-test)

# Create executable
add_executable(${WRK_TEST} ${TEST_SOURCES})
target_link_libraries(${WRK_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${WRK_TEST} PROPERTIES OUTPUT_NAME ${WRK_TEST}${EXE_EXT})

## create dynamically linked .so files
set(WRK_TEST_INPUT_FILE_BASE WRK_Input_1)

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
    ${MORIS_PACKAGE_DIR}/FEM/INT/src
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
    ${MORIS_PACKAGE_DIR}/MTK/src
    ${MORIS_PACKAGE_DIR}/PRM/src
    ${MORIS_PACKAGE_DIR}/MRS/COR/src
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/IP
    ${MORIS_PACKAGE_DIR}/PRM/ENM/src
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/CORE
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/IWG
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/CM
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/MM
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/SP
    ${MORIS_PACKAGE_DIR}/FEM/INT/src/ELEM
)
    
foreach(TPL ${SO_TPLS})
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

foreach(base_name ${WRK_TEST_INPUT_FILE_BASE})
    dynamic_link_input(${base_name} ${base_name} ${base_name}.cpp ${SO_INCLUDES})
endforeach()

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(WRK_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${WRK_TEST_PROCS})
        set(WRK_PARALLEL_TEST ${WRK_TEST}-${PROCS}-procs)
        add_test(NAME ${WRK_PARALLEL_TEST} 
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${WRK_TEST}>)
    endforeach()
endif()

