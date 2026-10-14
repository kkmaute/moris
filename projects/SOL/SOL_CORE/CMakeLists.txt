# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# SOL_CORE Solver Algebra

# Build the module library.
add_subdirectory(src)

## Add test directory # maybe add core tests when in debug
#if(MORIS_USE_TESTS AND BUILD_SOL_CORE)
#    add_subdirectory(test)
#endif()
