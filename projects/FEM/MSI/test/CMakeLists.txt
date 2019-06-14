# MSI ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_MSI_Dof_Manager_Test.cpp
    cl_MSI_Equation_Object_Test.cpp
    cl_MSI_Pdof_Host_Test.cpp
    cl_MSI_Model_Solver_Interface_Test.cpp
    cl_MSI_Multigrid.cpp
    cl_MSI_SpaceTime_Test.cpp
    )

# List test includes
set(TEST_INCLUDES
    #${SOL}/${DLA}
    #${SOL}/${NLA}
    #${COM}
    #${LINALG}
    #${ALG}
    #${MTK}
    #${FEM}/${INT}
    #${FEM}/${MSI}
    #${FEM}/${MDL}
    #${MTK}/${MAP}
    #${HMR}
)

#include_directories(${MORIS_INCDIRS})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies

set(TEST_DEPENDENCIES
	${LINALG}-lib
	${ALG}-lib
	${COM}-lib
	${NLA}-lib
	${MSI}-lib
	${MDL}-lib
    #${HMR}-lib
    #${MAP}-lib
    #${MTK}-lib
    #${MDL}-lib # for L2
    #${INT}-lib # for L2
    #${MSI}-lib # for L2
    #${NLA}-lib # for L2
#    ${SDF}-lib # for SDF
    #${COM}-lib
    ${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_TRILINOS_LIBS}
    #${MORIS_ARPACK_LIBS}
    #${MORIS_BOOST_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)
    
# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MSI_TEST ${MSI}-test)

# Create executable
add_executable(${MSI_TEST} ${TEST_SOURCES})
target_link_libraries(${MSI_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MSI_TEST} PROPERTIES OUTPUT_NAME ${MSI_TEST}${EXE_EXT})

foreach(TPL ${MSI_TEST_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MSI_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${MSI_TEST}
    	PRIVATE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${MSI_TEST}
    	PRIVATE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Mark as test
add_test(NAME ${MSI_TEST}
    COMMAND $<TARGET_FILE:${MSI_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MSI_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MSI_TEST_PROCS})
		set(MSI_PARALLEL_TEST ${MSI_TEST}-${PROCS}-procs)
		add_test(NAME ${MSI_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MSI_TEST}>)
	endforeach()
endif()




