# MSI -----------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#set(MSI_TPL_DEPENDENCIES     "armadillo"     "viennacl"     "PETSc"    #"trilinos"    ${MATH_LIB}    "mpi"    ${MATRIX_LIB})
list(REMOVE_DUPLICATES MSI_TPL_DEPENDENCIES)
foreach(MSI_TPL ${MSI_TPL_DEPENDENCIES})
    include(${MORIS_CMAKE_DIR}/${MSI_TPL}.cmake)
endforeach()

#include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MSI)
    add_subdirectory(test)
endif()
