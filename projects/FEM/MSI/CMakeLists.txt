#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MSI -----------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#set(MSI_TPL_DEPENDENCIES     "armadillo"     "viennacl"     "PETSc"    #"trilinos"    ${ACML_LAPACK_MKL_OPENBLAS}    "mpi"    ${ARMADILLO_EIGEN})
#list(REMOVE_DUPLICATES MSI_TPL_DEPENDENCIES)
#foreach(MSI_TPL ${MSI_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${MSI_TPL}.cmake)
#endforeach()

#include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MSI)
    add_subdirectory(test)
endif()

