# Model Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Bspline.cpp
    EgnManager.cpp )

# List test includes
set(TEST_INCLUDES
    #${COM}
    #${MOD}
    #${ALG}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    #${COM}-lib
    #${MOD}-lib
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MOD_TEST ${MOD}-test)

# Create executable
add_executable(${MOD_TEST} ${TEST_SOURCES})
target_link_libraries(${MOD_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${MOD_TEST} PROPERTIES OUTPUT_NAME ${MTK_TEST}${EXE_EXT})

foreach(TPL ${MOD_TEST_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MOD_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${MOD_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${MOD_TEST}
    	PRIVATE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Mark as test
add_test(NAME ${MOD_TEST}
    COMMAND $<TARGET_FILE:${MOD_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MOD_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MOD_TEST_PROCS})
		set(MOD_PARALLEL_TEST ${MOD_TEST}-${PROCS}-procs)
		add_test(NAME ${MOD_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MOD_TEST}>)
	endforeach()
endif()

