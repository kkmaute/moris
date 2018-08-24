# Interpolation and Integration -------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
list(REMOVE_DUPLICATES MDL_TPL_DEPENDENCIES)
foreach(MDL_TPL ${MDL_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${MDL_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MDL)
    add_subdirectory(test)
endif()