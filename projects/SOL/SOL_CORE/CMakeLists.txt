#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# SOL_CORE Solver Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES SOL_CORE_TPL_DEPENDENCIES)
#foreach(SOL_CORE_TPL ${SOL_CORE_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${SOL_CORE_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

## Add test directory # maybe add core tests when in debug
#if(MORIS_USE_TESTS AND BUILD_SOL_CORE)
#    add_subdirectory(test)
#endif()

