#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# OPENBLAS libraries ----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT OPENBLAS_FOUND_ONCE)
    find_package(OPENBLAS)
	
    if(OPENBLAS_FOUND)
        #set(OPENBLAS_FOUND_ONCE TRUE CACHE INTERNAL "OPENBLAS was found." FORCE)
        set(OPENBLAS_FOUND_ONCE TRUE)
        
        #set(MORIS_OPENBLAS_LIBRARIES ${OPENBLAS_LIBRARIES}
        #	CACHE INTERNAL "OPENBLAS libraries.")
        #set(MORIS_OPENBLAS_DEFINITIONS "-DMORIS_HAVE_OPENBLAS"
        #	CACHE INTERNAL "Moris preprocessor definitions for OPENBLAS.")
        set(MORIS_OPENBLAS_DEFINITIONS "-DMORIS_HAVE_OPENBLAS")
        set(MORIS_OPENBLAS_LIBRARIES OPENBLAS::openblas)
        
        mark_as_advanced(MORIS_OPENBLAS_LIBRARIES
        	MORIS_OPENBLAS_DEFINITIONS )
    endif()
    
    message(STATUS "OPENBLAS_LIBRARIES: ${OPENBLAS_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::openblas)
	add_library(${MORIS}::openblas INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::openblas INTERFACE ${MORIS_OPENBLAS_LIBRARIES})
endif()

#add_definitions("-DMORIS_HAVE_OPENBLAS")
#set(MORIS_OPENBLAS_LIBS ${LAPACK_LIBRARIES})
#set(MORIS_OPENBLAS_LAPACK_MKL_LIBS ${OPENBLAS_LIBRARIES})

