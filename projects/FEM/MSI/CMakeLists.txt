# MSI -----------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#set(MSI_TPL_DEPENDENCIES     "armadillo"     "viennacl"     "PETSc"    #"trilinos"    ${MATH_LIB}    "mpi"    ${MATRIX_LIB})
list(REMOVE_DUPLICATES MSI_TPL_DEPENDENCIES)
foreach(MSI_TPL ${MSI_TPL_DEPENDENCIES})
    include(${PROJECT_SOURCE_DIR}/share/cmake/${MSI_TPL}.cmake)
endforeach()

#include(${PROJECT_SOURCE_DIR}/share/cmake/MatthewCMake/DLA_Depends.cmake)

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MSI)
    add_subdirectory(test)
endif()
