# IOS ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_CMAKE_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES IOS_TPL_DEPENDENCIES)
foreach(IOS_TPL ${IOS_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${IOS_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_IOS)
    add_subdirectory(test)
endif()
