# Distributed Linear Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/PETSc.cmake)
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES DLA_TPL_DEPENDENCIES)
foreach(DLA_TPL ${DLA_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${DLA_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_DLA)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here
