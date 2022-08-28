#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# ARPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ARPACK_FOUND_ONCE)
    find_package(ARPACK)
    
    if(ARPACK_FOUND)
        set(ARPACK_FOUND_ONCE TRUE CACHE INTERNAL "ARPACK was found.")
        
        set(MORIS_ARPACK_LIBRARIES ${ARPACK_LIBRARIES}
        	CACHE INTERNAL "ARPACK libraries.")
        
        mark_as_advanced(MORIS_ARPACK_LIBRARIES)
    endif()
    
    message(STATUS "ARPACK_LIBRARIES: ${ARPACK_LIBRARIES}")
endif()

set(MORIS_ARPACK_LIBS ${ARPACK_LIBRARIES})

