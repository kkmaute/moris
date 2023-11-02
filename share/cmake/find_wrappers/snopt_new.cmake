#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# SNOPT libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT SNOPT_FOUND_ONCE)
    find_package(SNOPT)
    
    if(SNOPT_FOUND)
        #set(SNOPT_FOUND_ONCE TRUE CACHE INTERNAL "SNOPT was found.")
        set(SNOPT_FOUND_ONCE TRUE)
        
        #set(MORIS_SNOPT_INCLUDE_DIRS ${SNOPT_LIBRARY_DIRS}
        #	CACHE PATH "SNOPT library directories." )
        #set(MORIS_SNOPT_LIBRARIES ${SNOPT_LIBRARIES}
        #	CACHE PATH "SNOPT libraries." )
        set(MORIS_SNOPT_INCLUDE_DIRS ${SNOPT_LIBRARY_DIRS})
        set(MORIS_SNOPT_LIBRARIES ${SNOPT_LIBRARIES})
        
        mark_as_advanced(MORIS_SNOPT_INCLUDE_DIRS
        	MORIS_SNOPT_LIBRARIES )
    endif()
    message(STATUS "SNOPT_LIBRARIES: ${SNOPT_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::snopt)
	_import_libraries(SNOPT_LIBRARY_TARGETS ${SNOPT_LIBRARIES})

	add_library(${MORIS}::snopt INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::snopt INTERFACE ${SNOPT_LIBRARY_TARGETS} ${MORIS_LIB_GFortran} )
endif()

