# Linear Algebra ----------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_CMAKE_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES LNA_TPL_DEPENDENCIES)
foreach(LNA_TPL ${LNA_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${LNA_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_LNA)
    add_subdirectory(test)
endif()
