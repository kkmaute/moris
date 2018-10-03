# Interpolation and Integration -------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
list(REMOVE_DUPLICATES SDF_TPL_DEPENDENCIES)
foreach(SDF_TPL ${SDF_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${SDF_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_SDF)
#    add_subdirectory(test)
endif()