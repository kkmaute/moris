#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# ViennaCL libraries ------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT VIENNACL_FOUND_ONCE)
    find_package(ViennaCL)
    mark_as_advanced(ViennaCL_DIR)
    if(ViennaCL_FOUND)
        set(VIENNACL_FOUND_ONCE TRUE CACHE INTERNAL 
            "ViennaCL was found." FORCE)
        set(MORIS_VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIRS} CACHE INTERNAL "ViennaCL include directories" FORCE)
    endif()
    message(STATUS "VIENNACL_LIBRARIES: ${VIENNACL_LIBRARIES}")
endif()

include_directories(${MORIS_VIENNACL_INCLUDE_DIRS})


