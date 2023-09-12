#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Visualization Tests --------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    UT_VIS_Visualization_Mesh.cpp
    UT_VIS_Output_Manager.cpp
    UT_VIS_Output_Dof.cpp
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

# GEN directories
set(GEN_INCLUDES
	"field"
	"field/geometry"
	"field/property"
	"pdv"
	)

foreach(GEN_INCLUDES ${GEN_INCLUDES})
	include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()
    
# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)
    
# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${ALG}-lib
    ${VIS}-lib
    ${MTK}-lib
    ${HMR}-lib
    ${XTK}-lib
    ${GEN}-lib
    ${MORIS_BASE_LIBS}
    )    
    
# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(VIS_TEST ${VIS}-test)

# Create executable
add_executable(${VIS_TEST} ${TEST_SOURCES})
target_link_libraries(${VIS_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${VIS_TEST} PROPERTIES OUTPUT_NAME ${VIS_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(VIS_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${VIS_TEST_PROCS})
		set(VIS_PARALLEL_TEST ${VIS_TEST}-${PROCS}-procs)
		add_test(NAME ${VIS_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${VIS_TEST}>)
	endforeach()
endif()





