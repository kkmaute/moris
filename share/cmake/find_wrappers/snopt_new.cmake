# -------------------------------------------------------------------------
# SNOPT libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT SNOPT_FOUND_ONCE)
    find_package(SNOPT)
    
    if(SNOPT_FOUND)
        set(SNOPT_FOUND_ONCE TRUE CACHE INTERNAL "SNOPT was found.")
        
        set(MORIS_SNOPT_INCLUDE_DIRS ${SNOPT_LIBRARY_DIRS}
        	CACHE PATH "SNOPT library directories." )
        set(MORIS_SNOPT_LIBRARIES ${SNOPT_LIBRARIES}
        	CACHE PATH "SNOPT libraries." )
        
        mark_as_advanced(MORIS_SNOPT_INCLUDE_DIRS
        	MORIS_SNOPT_LIBRARIES )
    endif()
    message(STATUS "SNOPT_LIBRARIES: ${SNOPT_LIBRARIES}")
endif()

if(NOT TARGET snopt)
	add_library(snopt STATIC IMPORTED GLOBAL)
	set_target_properties(snopt PROPERTIES
		IMPORTED_LOCATION ${MORIS_SNOPT_LIBRARIES} )
endif()

#link_directories(${SNOPT_LIBRARY_DIRS})
#set(MORIS_SNOPT_LIBS ${SNOPT_LIBRARIES})
