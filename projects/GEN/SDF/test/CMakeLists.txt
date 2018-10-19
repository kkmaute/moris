# SDF Tests ---------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_SDF_Triangle.cpp 
    cl_SDF_Core.cpp)

# List test includes
set(TEST_INCLUDES
    ${GEN}/${SDF}
    ${MTK}
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
    ${SDF}-lib
    ${MTK}-lib
    ${COM}-lib
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_BASE_LIBS} )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${SDF}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${SDF}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${SDF}-test
    COMMAND ${SDF}-test${EXE_EXT} )

