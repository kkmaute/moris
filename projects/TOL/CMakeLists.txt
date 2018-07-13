# Tools -------------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_TOL)
    add_subdirectory(test)
endif()
