# Non Linear Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_CMAKE_DIR}/PETSc.cmake)
#include(${MORIS_CMAKE_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES NLA_TPL_DEPENDENCIES)
foreach(NLA_TPL ${NLA_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${NLA_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_NLA)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here
