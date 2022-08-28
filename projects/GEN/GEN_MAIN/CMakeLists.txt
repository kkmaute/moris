#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Geometry Engine ---------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#include(${MORIS_TPL_DIR}/PETSc.cmake)
#include(${MORIS_TPL_DIR}/viennacl.cmake)
#list(REMOVE_DUPLICATES GEN_TPL_DEPENDENCIES)
#foreach(GEN_TPL ${GEN_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${GEN_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

#if(BUILD_SDF)
#	add_subdirectory(SDF)
#endif()

# Add test directory
if(MORIS_USE_TESTS AND BUILD_GEN)
    add_subdirectory(test)
endif()

