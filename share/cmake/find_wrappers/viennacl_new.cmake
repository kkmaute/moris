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
    
    if(ViennaCL_FOUND)
        set(VIENNACL_FOUND_ONCE TRUE CACHE INTERNAL 
            "ViennaCL was found." FORCE)
        
        set(MORIS_VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIRS} 
        	CACHE PATH "ViennaCL include directories")
        
        mark_as_advanced(ViennaCL_DIR
        	MORIS_VIENNACL_INCLUDE_DIRS )
    endif()
    message(STATUS "VIENNACL_LIBRARIES: ${VIENNACL_LIBRARIES}")
endif()

if(NOT TARGET viennacl)
	add_library(viennacl INTERFACE IMPORTED GLOBAL)
endif()

#include_directories(${VIENNACL_INCLUDE_DIRS})


