# -----------------------------------------------------------------------------
# SuperLU libraries and includes ----------------------------------------------
# -----------------------------------------------------------------------------

if(NOT SUPERLU_FOUND_ONCE)
    find_package(SuperLU)
    if(SUPERLU_FOUND) # SuperLU_FOUND should be used here*
        set(SUPERLU_FOUND_ONCE TRUE
        	CACHE INTERNAL "SuperLU was found.")
        
        set(MORIS_SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDES}
        	CACHE PATH "SuperLU include directories.")
        set(MORIS_SUPERLU_LIBRARIES ${SUPERLU_LIBRARIES}
        	CACHE INTERNAL "SuperLU libraries.")
        
        mark_as_advanced(MORIS_SUPERLU_INCLUDE_DIRS
        	MORIS_SUPERLU_LIBRARIES )
        
        #add_library(superlu STATIC IMPORTED GLOBAL)
        #set_target_properties(superlu 
        #	PROPERTIES IMPORTED_LOCATION ${MORIS_SUPERLU_LIBRARIES} )
    endif()
    message(STATUS "SUPERLU_INCLUDES: ${SUPERLU_INCLUDES}")
	message(STATUS "SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")
endif()

if(NOT TARGET superlu)
	add_library(superlu STATIC IMPORTED GLOBAL)
	set_target_properties(superlu
		PROPERTIES IMPORTED_LOCATION ${MORIS_SUPERLU_LIBRARIES})
endif()

#include_directories(${SUPERLU_INCLUDES})
#set(MORIS_SUPERLU_LIBS ${SUPERLU_LIBRARIES})

# -------------------------------------------------------------------------
# SuperLU_DIST

if(NOT SUPERLU_DIST_FOUND_ONCE)
    find_package(SuperLU_DIST)
    if(SUPERLU_DIST_FOUND) # SuperLU_DIST_FOUND should be used here*
        set(SUPERLU_DIST_FOUND_ONCE TRUE
            CACHE INTERNAL "SuperLU_DIST was found." )
        
        #set(MORIS_SUPERLU_LIBRARIES ${MORIS_SUPERLU_LIBRARIES} ${SUPERLU_DIST_LIBRARIES}
        #	CACHE INTERNAL "SuperLU_DIST libraries." FORCE)
        set(MORIS_SUPERLU_DIST_LIBRARIES ${SUPERLU_DIST_LIBRARIES}
        	CACHE INTERNAL "SuperLU_DIST libraries." FORCE)
        
        mark_as_advanced(MORIS_SUPERLU_LIBRARIES)
        
        #add_library(superlu_dist STATIC IMPORTED GLOBAL)
        #set_target_properties(superlu_dist
       # 	PROPERTIES IMPORTED_LOCATION ${SUPERLU_DIST_LIBRARIES} )
    endif()
    message(STATUS "SUPERLU_DIST_LIBRARIES: ${SUPERLU_DIST_LIBRARIES}")
endif()

if(NOT TARGET superlu_dist)
	add_library(superlu_dist STATIC IMPORTED GLOBAL)
	set_target_properties(superlu_dist
		PROPERTIES IMPORTED_LOCATION ${MORIS_SUPERLU_DIST_LIBRARIES} )
	
	if(NOT TARGET superlu)
		add_library(superlu INTERFACE IMPORTED GLOBAL)
	endif()
	target_link_libraries(superlu INTERFACE superlu_dist)
endif()

#list(APPEND MORIS_SUPERLU_LIBS ${SUPERLU_DIST_LIBRARIES})

# *variable only gets set for one or the other for some reason... 
# Using the librairies variable as a workaround
