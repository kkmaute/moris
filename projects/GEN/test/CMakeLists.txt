# Geometry Engine Tests ---------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_GeometryEngineFactory.cpp
    cl_GeometryObject.cpp
    cl_GeSDFCore.cpp
    cl_GeSDFGen.cpp
    cl_GeTriangle.cpp )

# List test includes
set(TEST_INCLUDES
	${LINALG}
    ${COM})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${GEN}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${GEN}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${GEN}-test
    COMMAND ${GEN}-test${EXE_EXT} )

