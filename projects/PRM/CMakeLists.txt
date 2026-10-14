# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Parameters

# Build the module library.
add_subdirectory(src)
#add_subdirectory(ENM)

# Build tests when this module is selected.
if(MORIS_USE_TESTS AND BUILD_PRM)
    add_subdirectory(test)
endif()
