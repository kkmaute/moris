# Distributed Linear Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
list(REMOVE_DUPLICATES DLA_TPL_DEPENDENCIES)
foreach(DLA_TPL ${DLA_TPL_DEPENDENCIES})
    include(${MORIS_TPL_DIR}/${DLA_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_DLA)
    add_subdirectory(test)
endif()
