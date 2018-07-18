# STK ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/PETSc.cmake)
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES STK_TPL_DEPENDENCIES)
foreach(STK_TPL ${STK_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${STK_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_STK)
    add_subdirectory(test)
endif()
