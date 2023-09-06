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
         set(SUITESPARSE_FOUND_ONCE TRUE 
            CACHE INTERNAL "SuiteSparse was found.")

        set(MORIS_SUITESPARSE_INCLUDE_DIRS ${SUITESPARSE_INCLUDE_DIRS} 
            CACHE INTERNAL "SuiteSparse include directories" )

        set(MORIS_SUITESPARSE_LIBRARIES ${SUITESPARSE_LIBRARIES}
            CACHE INTERNAL "SuiteSparse libraries.")

        mark_as_advanced(SuiteSparse_DIR
            MORIS_SUITESPARSE_INCLUDE_DIRS
            MORIS_SUITESPARSE_LIBRARIES
            )
    endif()
    message(STATUS "SUITESPARSE_INCLUDE_DIRS: ${SUITESPARSE_INCLUDE_DIRS}")
    message(STATUS "SUITESPARSE_LIBRARIES:    ${SUITESPARSE_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::suitesparse)
    add_library(${MORIS}::suitesparse STATIC IMPORTED GLOBAL)
    set_target_properties(${MORIS}::suitesparse
        PROPERTIES IMPORTED_LOCATION ${MORIS_SUITESPARSE_LIBRARIES})
endif()

#include_directories(${SUITESPARSE_INCLUDE_DIRS})

