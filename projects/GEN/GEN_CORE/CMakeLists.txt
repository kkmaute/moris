#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# GEN_CORE Geometry Engine ------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES GEN_CORE_TPL_DEPENDENCIES)
#foreach(GEN_CORE_TPL ${GEN_CORE_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${GEN_CORE_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

## Add test directory # maybe add core tests when in debug
#if(MORIS_USE_TESTS AND BUILD_SOL_CORE)
#    add_subdirectory(test)
#endif()

