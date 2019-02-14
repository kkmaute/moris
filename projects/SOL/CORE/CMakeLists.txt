# Core Solver Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
list(REMOVE_DUPLICATES CORE_TPL_DEPENDENCIES)
foreach(CORE_TPL ${CORE_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${CORE_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_CORE)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here
