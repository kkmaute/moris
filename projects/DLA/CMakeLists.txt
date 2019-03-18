# Distributed Linear Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_CMAKE_DIR}/PETSc.cmake)
#include(${MORIS_CMAKE_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES DLA_TPL_DEPENDENCIES)
foreach(DLA_TPL ${DLA_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${DLA_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_DLA)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here
