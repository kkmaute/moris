# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Compatibility entry point; module selection is centralized.
include("${CMAKE_CURRENT_LIST_DIR}/moris_select_modules.cmake")
moris_set_build_option(BUILD_ALL ON)
moris_select_modules()
