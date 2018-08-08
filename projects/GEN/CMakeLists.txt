# Geometry Engine ---------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_CMAKE_DIR}/PETSc.cmake)
#include(${MORIS_CMAKE_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES GEN_TPL_DEPENDENCIES)
foreach(GEN_TPL ${GEN_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${GEN_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_GEN)
    add_subdirectory(test)
endif()
