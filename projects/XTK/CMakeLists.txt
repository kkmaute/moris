#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# XTK ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Based on MORIS Compiler Options, set XTK Options

SET(XTK_DEFINITIONS)

#if(MORIS_USE_EIGEN)
#        LIST(APPEND XTK_DEFINITIONS "-DXTK_USE_EIGEN")
#else()
#        message(FATAL_ERROR "XTK only uses Eigen" )
#endif()


ADD_DEFINITIONS(${XTK_DEFINITIONS})

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_XTK)
    add_subdirectory(test)
endif()

if( BUILD_TUTORIALS AND BUILD_XTK )
# Add source directory
add_subdirectory(tutorials)
endif()


