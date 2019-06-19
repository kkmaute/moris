# -------------------------------------------------------------------------
# GCMMA libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT GCMMA_FOUND_ONCE)
    find_package(GCMMA)

	if(GCMMA_FOUND)
        set(GCMMA_FOUND_ONCE TRUE CACHE INTERNAL "GCMMA was found.")
        
        set(MORIS_GCMMA_INCLUDE_DIRS ${GCMMA_INCLUDE_DIRS}
        	CACHE PATH "GCMMA include directories." )
        set(MORIS_GCMMA_LIBRARY_DIRS ${GCMMA_LIBRARY_DIRS}
        	CACHE PATH "GCMMA library directories." )
        set(MORIS_GCMMA_LIBRARIES ${GCMMA_LIBRARIES}
        	CACHE PATH "GCMMA libraries." )
        
        mark_as_advanced(MORIS_GCMMA_INCLUDE_DIRS
        	MORIS_GCMMA_LIBRARY_DIRS
        	MORIS_GCMMA_LIBRARIES )
    endif()

    message(STATUS "GCMMA_LIBRARIES: ${GCMMA_LIBRARIES}")
endif()

if(NOT TARGET gcmma)
	add_library(gcmma STATIC IMPORTED GLOBAL)
	set_target_properties(gcmma PROPERTIES
		IMPORTED_LOCATION ${MORIS_GCMMA_LIBRARIES} )
endif()

#include_directories(${GCMMA_INCLUDE_DIRS})
#link_directories(${GCMMA_LIBRARY_DIRS})
#set(MORIS_GCMMA_LIBS ${GCMMA_LIBRARIES})
