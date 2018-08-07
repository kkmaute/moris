# HMR ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
list(REMOVE_DUPLICATES HMR_TPL_DEPENDENCIES)
foreach(HMR_TPL ${HMR_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${HMR_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_HMR)
    add_subdirectory(test)
endif()
