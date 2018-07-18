# Communications ----------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${PROJECT_SOURCE_DIR}/share/cmake/viennacl.cmake)
list(REMOVE_DUPLICATES COM_TPL_DEPENDENCIES)
foreach(COM_TPL ${COM_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${COM_TPL}.cmake)
endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_COM)
    add_subdirectory(test)
endif()
