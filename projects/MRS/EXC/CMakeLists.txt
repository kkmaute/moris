#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Exceptions --------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/viennacl.cmake)
#list(REMOVE_DUPLICATES EXC_TPL_DEPENDENCIES)
#foreach(EXC_TPL ${EXC_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${EXC_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_EXC)
    add_subdirectory(test)
endif()

