# Geometry Engine Tests ---------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
	test_main.cpp
	ut_GEN_Pdv.cpp
	ut_GEN_Geometry.cpp
	ut_GEN_Property.cpp
	ut_GEN_Intersection.cpp
	fn_GEN_create_simple_mesh.cpp
    )
    
# List header dependencies                                                                                                                                                                                                                                                     
set(GEN_MAIN_INCLUDES                                                                                                                                                                                                                                                               
    "additional"
    "pdv"
    "geometry"
    "property"
    )

# Include depedency directories
foreach(GEN_MAIN_INCLUDE ${GEN_MAIN_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_MAIN_INCLUDE})
endforeach()
     
# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${GEN_CORE}-lib
    ${COR}-lib
    ${ALG}-lib
    ${GEN_MAIN}-lib
    ${INT}-lib
    ${SOL_CORE}-lib
	${HMR}-lib
	${MORIS_BASE_LIBS}
	${PRM}-lib
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
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${GEN_TEST}>)
	endforeach()
endif()

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
