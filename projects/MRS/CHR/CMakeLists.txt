#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Chronos -----------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES CHR_TPL_DEPENDENCIES)
#foreach(CHR_TPL ${CHR_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${CHR_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_CHR)
    add_subdirectory(test)
endif()

