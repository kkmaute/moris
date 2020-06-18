# Parameters Tests ----------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES 
    test_main.cpp
    UT_WRK_Performer_Manager.cpp
    UT_WRK_Workflow.cpp
    UT_WRK_L2.cpp )

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${WRK}-lib
    ${XTK}-lib
    ${HMR}-lib
    )
    
    set(XTK_INCLUDES
    "assert"
    "containers"
    "core"
    "geomeng"
    "geometry"
    "ios"
    "linalg"
    "mesh"
    "tools"
    "topology"
    "xtk"                                                                                                                                                                                                                                                                      
    )  
    
    # Include depedency directories                                                                                                                                                                                                                                                
foreach(XTK_INCLUDE ${XTK_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})                                                                                                                                                                                                           
endforeach()  

set(GEN_INCLUDES                                                                                                                                                                                                                                                               
    "additional"
    "geomeng"
    "geometry"
    "property"
    )

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
endforeach()

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(WRK_TEST ${WRK}-test)

# Create executable
add_executable(${WRK_TEST} ${TEST_SOURCES})
target_link_libraries(${WRK_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${WRK_TEST} PROPERTIES OUTPUT_NAME ${WRK_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(WRK_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${WRK_TEST_PROCS})
        set(WRK_PARALLEL_TEST ${WRK_TEST}-${PROCS}-procs)
        add_test(NAME ${WRK_PARALLEL_TEST} 
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${WRK_TEST}>)
    endforeach()
endif()


## create dynamically linked .so files
set(GEN_TEST_INPUT_FILE_BASE WRK_Input_1)

set(SO_TPLS
"trilinos"
${ARMADILLO_EIGEN}
)

SET(SO_INCLUDES 
	${Boost_INCLUDE_DIRS}
	${Trilinos_INCLUDE_DIRS}
	${MORIS_PACKAGE_DIR}/ALG/src
	${MORIS_PACKAGE_DIR}/FEM/INT/src 
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
    ${MORIS_PACKAGE_DIR}/MRS/COR/src)
    
    
foreach(TPL ${SO_TPLS})
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

foreach(base_name ${GEN_TEST_INPUT_FILE_BASE})
    dynamic_link_input(${base_name} ${base_name} ${base_name}.cpp ${SO_INCLUDES})
endforeach()
