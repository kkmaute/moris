# MIG ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES MTK_TPL_DEPENDENCIES)
#foreach(MIG_TPL ${MIG_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${MIG_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MIG)
    add_subdirectory(test)
endif()
