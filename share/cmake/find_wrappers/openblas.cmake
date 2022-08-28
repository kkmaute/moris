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
    find_package(OPENBLAS)
    
    if(OPENBLAS_FOUND)
        set(OPENBLAS_FOUND_ONCE TRUE CACHE INTERNAL "OPENBLAS was found." FORCE)
        
        set(MORIS_ACML_LIBRARIES ${OPENBLAS_LIBRARIES}
        	CACHE INTERNAL "OPENBLAS libraries.")
        set(MORIS_ACML_DEFINITIONS "-DMORIS_HAVE_OPENBLAS"
        	CACHE INTERNAL "Moris preprocessor definitions for OPENBLAS.")
        
        mark_as_advanced(MORIS_OPENBLAS_LIBRARIES
        	MORIS_OPENBLAS_DEFINITIONS )
    endif()
    
    message(STATUS "OPENBLAS_LIBRARIES: ${OPENBLAS_LIBRARIES}")
endif()

add_definitions("-DMORIS_HAVE_OPENBLAS")
set(MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS ${OPENBLAS_LIBRARIES})

