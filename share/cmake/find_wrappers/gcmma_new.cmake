#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# GCMMA libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT GCMMA_FOUND_ONCE)
    find_package(GCMMA)

	if(GCMMA_FOUND)
        #set(GCMMA_FOUND_ONCE TRUE CACHE INTERNAL "GCMMA was found.")
        set(GCMMA_FOUND_ONCE TRUE)
        
        #set(MORIS_GCMMA_INCLUDE_DIRS ${GCMMA_INCLUDE_DIRS}
        #	CACHE PATH "GCMMA include directories." )
        #set(MORIS_GCMMA_LIBRARY_DIRS ${GCMMA_LIBRARY_DIRS}
        #	CACHE PATH "GCMMA library directories." )
        #set(MORIS_GCMMA_LIBRARIES ${GCMMA_LIBRARIES}
        #	CACHE PATH "GCMMA libraries." )
        set(MORIS_GCMMA_INCLUDE_DIRS ${GCMMA_INCLUDE_DIRS})
        set(MORIS_GCMMA_LIBRARY_DIRS ${GCMMA_LIBRARY_DIRS})
        set(MORIS_GCMMA_LIBRARIES ${GCMMA_LIBRARIES})
        
        mark_as_advanced(MORIS_GCMMA_INCLUDE_DIRS
        	MORIS_GCMMA_LIBRARY_DIRS
        	MORIS_GCMMA_LIBRARIES )
    endif()

    message(STATUS "GCMMA_LIBRARIES: ${GCMMA_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::gcmma)
	_import_libraries(GCMMA_LIBRARY_TARGETS ${GCMMA_LIBRARIES})
	
	add_library(${MORIS}::gcmma INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::gcmma INTERFACE ${GCMMA_LIBRARY_TARGETS})
endif()

