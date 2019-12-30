# Geometry Engine Tests ---------------------------------------------------
# -------------------------------------------------------------------------
    include_directories(GE_Test_Proxy)

# List source files
set(TEST_SOURCES
    test_main.cpp
  #  cl_Fiber_Problem.cpp
    cl_GeSDFCore.cpp
    cl_GeSDFGen.cpp
    cl_GeTriangle.cpp
    cl_GE_test.cpp
    UT_Bounding_Box.cpp
  #  UT_Finite_Difference_Checks.cpp
    UT_GE_HMR_Interaction.cpp
    UT_GE_Properties.cpp
    
   # geometry_classes/cl_GE_test_analytic.cpp
   # geometry_classes/cl_GE_test_discrete.cpp
   # geometry_classes/cl_GE_test_sdf.cpp

	new_impl_test/cl_General_Functionalities.cpp
   # new_impl_test/cl_Parameter_Sweep.cpp
    new_impl_test/cl_Pdv_Tests.cpp
    )

#FIXME: Remove this once xtk has removed subdirectories from its package
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


    
# List header dependencies                                                                                                                                                                                                                                                     
set(GEN_INCLUDES                                                                                                                                                                                                                                                               
    "additional"
    "field"
    "geomeng"
    "geometry"
    "property"
    )

# Include depedency directories
foreach(GEN_INCLUDE ${GEN_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDE})
endforeach()
     

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${ALG}-lib
    ${GEN}-lib
    ${INT}-lib
  # ${MTK}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(GEN_TEST ${GEN}-test)

# Create executable
add_executable(${GEN_TEST} ${TEST_SOURCES})
target_link_libraries(${GEN_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${GEN_TEST} PROPERTIES OUTPUT_NAME ${GEN_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(GEN_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${GEN_TEST_PROCS})
		set(GEN_PARALLEL_TEST ${GEN_TEST}-${PROCS}-procs)
		add_test(NAME ${GEN_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} $<TARGET_FILE:${GEN_TEST}>)
	endforeach()
endif()

