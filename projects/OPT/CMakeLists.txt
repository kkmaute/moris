# Optimization ------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES OPT_TPL_DEPENDENCIES)
foreach(OPT_TPL ${OPT_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${OPT_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_OPT)
    add_subdirectory(test)
endif()
