# Communications ----------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/viennacl.cmake)
#list(REMOVE_DUPLICATES COM_TPL_DEPENDENCIES)
#foreach(COM_TPL ${COM_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${COM_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_COM)
    add_subdirectory(test)
endif()
