# Integration and Interpolation Tests -------------------------------------
# -------------------------------------------------------------------------
include_directories(DLA_Test_Proxys)

# List source files
set(TEST_SOURCES
    cl_MDL_DiffusionElement_Mesh_Interaction_Test.cpp
    cl_MDL_DiffusionElement_Cut_Mesh.cpp
    UT_MDL_Mesh_Model_Helper.cpp
    UT_MDL_XTK_DiffusionElement.cpp
    UT_MDL_XTK_HMR_DiffusionElement.cpp
    test_main.cpp
    #MDL_Test_Proxys/cl_MDL_Element_Proxy.hpp
    #MDL_Test_Proxys/cl_MDL_Node_Proxy.hpp
    )

# List test includes
set(TEST_INCLUDES
    #${LINALG}
    #${IOS}
    #${COM}
    #${MTK}
    #${HMR}
    #${ALG}
    #${SOL}/${DLA}
    #${SOL}/${NLA}
    #${SOL}/${TSA}
    #${SOL}/${SOL_CORE}
    #${FEM}/${INT}
    #${FEM}/${MSI}
    #${FEM}/${MDL}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

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

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
	${MTK}-lib
	${ALG}-lib
    #${COM}-lib
    #${MDL}-lib
    #${MSI}-lib
    ${INT}-lib
    #${MTK}-lib
    #${TOL}-lib
    ${HMR}-lib
    ${SDF}-lib
    ${MAP}-lib
    ${XTK}-lib
    #${DLA}-lib
    #${NLA}-lib
    #${TSA}-lib
    #${SOL_CORE}-lib
    ${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MDL_TEST ${MDL}-test)

# Create executable
add_executable(${MDL_TEST} ${TEST_SOURCES})
target_link_libraries(${MDL_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MDL_TEST} PROPERTIES OUTPUT_NAME ${MDL_TEST}${EXE_EXT})

# Mark as test
add_test(NAME ${MDL_TEST}
    COMMAND $<TARGET_FILE:${MDL_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(MDL_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${MDL_TEST_PROCS})
        set(MDL_PARALLEL_TEST ${MDL_TEST}-${PROCS}-procs)
        add_test(NAME ${MDL_PARALLEL_TEST} 
            COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MDL_TEST}>)
    endforeach()
endif()

