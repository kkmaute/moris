#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Optimization ------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/viennacl.cmake)
#list(REMOVE_DUPLICATES OPT_TPL_DEPENDENCIES)
#foreach(OPT_TPL ${OPT_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${OPT_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_OPT)
    add_subdirectory(test)
endif()

