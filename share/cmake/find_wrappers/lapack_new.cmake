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
set(LAPACK_ENV_VARS
    $ENV{LAPACKDIR}
    $ENV{LAPACK_DIR}
    $ENV{lapack_DIR}
    $ENV{LAPACK_ROOT}
    $ENV{lapack_ROOT}
    $ENV{LAPACK_PATH}
    $ENV{lapack_PATH} )
    
    #find_package(LAPACK REQUIRED PATHS ${LAPACK_ENV_VARS})
    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND)
        set(LAPACK_FOUND_ONCE TRUE)
        
        #set(MORIS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES}
        #	CACHE INTERNAL "LAPACK libraries.")
        #set(MORIS_LAPACK_DEFINTIONS "-DMORIS_HAVE_LAPACK"
        #	CACHE INTERNAL "Moris preprocessor definitions for LAPACK.")
        set(MORIS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
        set(MORIS_LAPACK_DEFINTIONS "-DMORIS_HAVE_LAPACK")
        
        mark_as_advanced(MORIS_LAPACK_LIBRARIES
        	MORIS_LAPACK_DEFINTIONS )
    endif()
    message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::lapack)
	_import_libraries(LAPACK_LIBRARY_TARGETS "${MORIS_LAPACK_LIBRARIES}")

	add_library(${MORIS}::lapack INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::lapack INTERFACE ${LAPACK_LIBRARY_TARGETS})
endif()

#add_definitions("-DMORIS_HAVE_LAPACK")
#set(MORIS_LAPACK_LIBS ${LAPACK_LIBRARIES})
#set(MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS ${LAPACK_LIBRARIES})

