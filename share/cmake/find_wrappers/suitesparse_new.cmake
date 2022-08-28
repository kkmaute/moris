#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# SuiteSparse libraries ------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT SUITESPARSE_FOUND_ONCE)
    find_package(SuiteSparse)
    
    if(SuiteSparse_FOUND)
        set(SUITESPARSE_FOUND_ONCE TRUE CACHE INTERNAL 
            "SuiteSparse was found." FORCE)
        
        set(MORIS_SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} 
        	CACHE INTERNAL "SuiteSparse include directories" FORCE )
        
        mark_as_advanced(SuiteSparse_DIR
        	MORIS_SUITESPARSE_INCLUDE_DIRS
        	)
    endif()
    message(STATUS "SUITESPARSE_LIBRARIES: ${SUITESPARSE_LIBRARIES}")
endif()

if(NOT TARGET suitesparse)
	add_library(suitesparse INTERFACE IMPORTED GLOBAL)
endif()

#include_directories(${SUITESPARSE_INCLUDE_DIRS})

