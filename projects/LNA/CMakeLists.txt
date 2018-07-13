# Linear Algebra ----------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_LNA)
    add_subdirectory(test)
endif()
