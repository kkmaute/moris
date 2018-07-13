# Distributed Linear Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_COM)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here
