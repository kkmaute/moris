# Assert ------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/viennacl.cmake)
list(REMOVE_DUPLICATES ASR_TPL_DEPENDENCIES)
foreach(ASR_TPL ${ASR_TPL_DEPENDENCIES})
    include(${MORIS_TPL_DIR}/${ASR_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_ASR)
    add_subdirectory(test)
endif()
