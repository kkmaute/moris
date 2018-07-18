# Linear Algebra ----------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES LNA_TPL_DEPENDENCIES)
foreach(LNA_TPL ${LNA_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${LNA_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_LNA)
    add_subdirectory(test)
endif()
