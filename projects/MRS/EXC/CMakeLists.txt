# Exceptions --------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_EXC)
    add_subdirectory(test)
endif()
