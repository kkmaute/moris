# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# XTK

# Based on MORIS Compiler Options, set XTK Options

set(XTK_DEFINITIONS)

#if(MORIS_USE_EIGEN)
#        LIST(APPEND XTK_DEFINITIONS "-DXTK_USE_EIGEN")
#else()
#        message(FATAL_ERROR "XTK only uses Eigen" )
#endif()

add_definitions(${XTK_DEFINITIONS})

# Build the module library.
add_subdirectory(src)

# Build tests when this module is selected.
if(MORIS_USE_TESTS AND BUILD_XTK)
    add_subdirectory(test)
endif()

if(BUILD_TUTORIALS AND BUILD_XTK)
    # Build the module library.
    add_subdirectory(tutorials)
endif()
