# Integration and Interpolation Tests -------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${IOS}
    ${COM}
    ${MTK}
    ${HMR}
    ${FEM}/${INT}
    ${FEM}/${MSI}
    ${FEM}/${MDL}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MDL}-lib
    ${MSI}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${MDL}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MDL}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MDL}-test
    COMMAND ${MDL}-test${EXE_EXT} )

