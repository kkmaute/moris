# Algorithms --------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES ALG_TPL_DEPENDENCIES)
foreach(ALG_TPL ${ALG_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${ALG_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_ALG)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here
