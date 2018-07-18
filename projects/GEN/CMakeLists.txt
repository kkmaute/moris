# Geometry Engine ---------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/PETSc.cmake)
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES GEN_TPL_DEPENDENCIES)
foreach(GEN_TPL ${GEN_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${GEN_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_GEN)
    add_subdirectory(test)
endif()
