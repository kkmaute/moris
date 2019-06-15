# Optimization Tests ------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    opt_test.cpp )

# List test includes
set(TEST_INCLUDES
    #${LINALG}
    #${COM}
    #${OPT}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${OPT}-lib
    ${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    #${MORIS_GCMMA_LIBS}
    #${MORIS_SNOPT_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(OPT_TEST ${OPT}-test)

# Create executable
add_executable(${OPT_TEST} ${TEST_SOURCES})
target_link_libraries(${OPT_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${MOPT_TEST} PROPERTIES OUTPUT_NAME ${OPT_TEST}${EXE_EXT})

foreach(TPL ${OPT_TEST_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${OPT_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${OPT_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${OPT_TEST}
    	PRIVATE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Mark as test
add_test(NAME ${OPT_TEST}
    COMMAND $<TARGET_FILE:${OPT_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(OPT_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${OPT_TEST_PROCS})
		set(OPT_PARALLEL_TEST ${OPT_TEST}-${PROCS}-procs)
		add_test(NAME ${OPT_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${OPT_TEST}>)
	endforeach()
endif()

