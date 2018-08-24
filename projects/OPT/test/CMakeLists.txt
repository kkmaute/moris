# Optimization Tests ------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    opt_test.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${COM}
    ${OPT} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${OPT}-lib
    ${MORIS_LDLIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_GCMMA_LIBS} )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${OPT}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${OPT}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${OPT}-test
    COMMAND ${OPT}-test${EXE_EXT} )

