# Containers --------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES CON_TPL_DEPENDENCIES)
foreach(CON_TPL ${CON_TPL_DEPENDENCIES})
    include(${MORIS_TPL_DIR}/${CON_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_CON)
    add_subdirectory(test)
endif()
