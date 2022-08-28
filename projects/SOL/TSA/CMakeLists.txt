#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/PETSc.cmake)
#include(${MORIS_TPL_DIR}/viennacl.cmake)
#list(REMOVE_DUPLICATES TSA_TPL_DEPENDENCIES)
#foreach(TSA_TPL ${TSA_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${TSA_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_TSA)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here

