#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# LBFGSB libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LBFGSB_FOUND_ONCE)
    find_package(LBFGSB)
    
    if(LBFGSB_FOUND)
        #set(LBFGSB_FOUND_ONCE TRUE CACHE INTERNAL "LBFGSB was found.")
        set(LBFGSB_FOUND_ONCE TRUE)
        
        #set(MORIS_LBFGSB_INCLUDE_DIRS ${LBFGSB_LIBRARY_DIRS}
        #	CACHE PATH "LBFGSB library directories." )
        #set(MORIS_LBFGSB_LIBRARIES ${LBFGSB_LIBRARIES}
        #	CACHE PATH "LBFGSB libraries." )
        set(MORIS_LBFGSB_INCLUDE_DIRS ${LBFGSB_LIBRARY_DIRS})
        set(MORIS_LBFGSB_LIBRARIES ${LBFGSB_LIBRARIES})
        
        mark_as_advanced(MORIS_LBFGSB_INCLUDE_DIRS
        	MORIS_LBFGSB_LIBRARIES )
    endif()
    message(STATUS "LBFGSB_LIBRARIES: ${LBFGSB_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::lbfgsb)
	_import_libraries(LBFGSB_LIBRARY_TARGETS ${LBFGSB_LIBRARIES})

	add_library(${MORIS}::lbfgsb INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::lbfgsb INTERFACE  ${LBFGSB_LIBRARY_TARGETS} ${MORIS_LIB_GFortran} )
endif()

