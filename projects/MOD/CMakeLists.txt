# Model -------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES MOD_TPL_DEPENDENCIES)
foreach(MOD_TPL ${MOD_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${MOD_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MOD)
    add_subdirectory(test)
endif()
