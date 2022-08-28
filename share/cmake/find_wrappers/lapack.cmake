#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# LAPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LAPACK_FOUND_ONCE)
    find_package(LAPACK)
    if(LAPACK_FOUND)
        set(LAPACK_FOUND_ONCE TRUE CACHE INTERNAL "LAPACK was found.")
        
        set(MORIS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES}
        	CACHE INTERNAL "LAPACK libraries.")
        set(MORIS_LAPACK_DEFINTIONS "-DMORIS_HAVE_LAPACK"
        	CACHE INTERNAL "Moris preprocessor definitions for LAPACK.")
        
        mark_as_advanced(MORIS_LAPACK_LIBRARIES
        	MORIS_LAPACK_DEFINTIONS )
    endif()
    message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
endif()

add_definitions("-DMORIS_HAVE_LAPACK")
set(MORIS_LAPACK_LIBS ${LAPACK_LIBRARIES})
set(MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS ${LAPACK_LIBRARIES})

