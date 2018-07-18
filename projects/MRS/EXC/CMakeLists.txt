# Exceptions --------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES EXC_TPL_DEPENDENCIES)
foreach(EXC_TPL ${EXC_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${EXC_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_EXC)
    add_subdirectory(test)
endif()
