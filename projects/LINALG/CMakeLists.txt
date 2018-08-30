# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------
set(TEST_SOURCES
    test_main.cpp
    fn_det.cpp)
    

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${ALG})


# Third party library dependencies
list(REMOVE_DUPLICATES LINALG_TPL_DEPENDENCIES)
foreach(LINALG_TPL ${LINALG_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${LINALG_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_LINALG)
    add_subdirectory(test)
endif()
