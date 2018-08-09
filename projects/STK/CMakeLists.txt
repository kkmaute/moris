# STK ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_CMAKE_DIR}/PETSc.cmake)
#include(${MORIS_CMAKE_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES STK_TPL_DEPENDENCIES)
foreach(STK_TPL ${STK_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${STK_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_STK)
    add_subdirectory(test)
endif()
