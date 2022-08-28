#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# ACML libraries ----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ACML_FOUND_ONCE)
    find_package(ACML)
    
    if(ACML_FOUND)
        set(ACML_FOUND_ONCE TRUE CACHE INTERNAL "ACML was found." FORCE)
        
        set(MORIS_ACML_LIBRARIES ${ACML_LIBRARIES}
        	CACHE INTERNAL "ACML libraries.")
        set(MORIS_ACML_DEFINITIONS "-DMORIS_HAVE_ACML"
        	CACHE INTERNAL "Moris preprocessor definitions for ACML.")
        
        mark_as_advanced(MORIS_ACML_LIBRARIES
        	MORIS_ACML_DEFINITIONS )
    endif()
    
    message(STATUS "ACML_LIBRARIES: ${ACML_LIBRARIES}")
endif()

add_definitions("-DMORIS_HAVE_ACML")
#set(MORIS_ACML_LIBS ${LAPACK_LIBRARIES})
set(MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS ${ACML_LIBRARIES})

