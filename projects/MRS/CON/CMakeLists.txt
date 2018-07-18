# Containers --------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES CON_TPL_DEPENDENCIES)
foreach(CON_TPL ${CON_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${CON_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_CON)
    add_subdirectory(test)
endif()
