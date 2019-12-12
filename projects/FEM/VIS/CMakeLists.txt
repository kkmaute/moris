# VIS -----------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES VIS_TPL_DEPENDENCIES)
#foreach(VIS_TPL ${VIS_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${VIS_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_VIS)
    add_subdirectory(test)
endif()
