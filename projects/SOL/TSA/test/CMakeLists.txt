# Time Solver Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

    include_directories(TSA_Test_Proxy)
    
# List source files
set(TEST_SOURCES
    cl_TSA_Time_Solver_Tests.cpp
    cl_TSA_Time_Solver_Staggered_Tests.cpp
    TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy.cpp
    TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy2.cpp
    test_main.cpp)

# List test includes
set(TEST_INCLUDES
    #${SOL}/${DLA}
    #${SOL}/${SOL_CORE}
    #${LINALG}
    #${COM}
    #${SOL}/${NLA}
    #${SOL}/${TSA}
    #${ALG}
    #${FEM}/${MSI}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
	${ALG}-lib
	${NLA}-lib
	${MSI}-lib
    #${DLA}-lib
    #${SOL_CORE}-lib
    #${COM}-lib
    #${NLA}-lib
    #${TSA}-lib
    #${MSI}-lib
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(TSA_TEST ${TSA}-test)

# Create executable
add_executable(${TSA_TEST} ${TEST_SOURCES})
target_link_libraries(${TSA_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${TSA_TEST} PROPERTIES OUTPUT_NAME ${TSA_TEST}${EXE_EXT})

foreach(TPL ${TSA_TEST_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${TSA_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${TSA_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${TSA_TEST}
    	PRIVATE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Mark as test
add_test(NAME ${TSA_TEST}
    COMMAND $<TARGET_FILE:${TSA_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(TSA_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${TSA_TEST_PROCS})
        set(TSA_PARALLEL_TEST ${TSA_TEST}-${PROCS}-procs)
        add_test(NAME ${TSA_PARALLEL_TEST} 
            COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${TSA_TEST}>)
    endforeach()
endif()

