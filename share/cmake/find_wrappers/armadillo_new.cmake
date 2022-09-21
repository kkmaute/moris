#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# Armadillo libraries -----------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ARMADILLO_FOUND_ONCE)
    set(ARMADILLO_ENV_VARS
        $ENV{ARMADILLODIR}
        $ENV{ARMADILLO_DIR}
        $ENV{Armadillo_DIR}
        $ENV{ARMADILLO_ROOT}
        $ENV{Armadillo_ROOT}
        $ENV{ARMADILLO_PATH}
        $ENV{Armadillo_PATH} )

    find_package(Armadillo REQUIRED HINTS ${ARMADILLO_ENV_VARS})
    find_library(ARMADILLO_LIBRARIES NAMES libarmadillo HINTS ${ARMADILLO_LIBRARY_DIRS} )
    
    #set(MORIS_ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS}
    #    CACHE PATH "Armadillo include directories." )
    ##set(MORIS_ARMADILLO_LIBRARY_DIRS ${ARMADILLO_LIBRARY_DIRS}
    ##    CACHE PATH "Armadillo library directories." )
    ##set(MORIS_ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES}
    ##    CACHE FILEPATH "Armadillo libraries." )
    #set(MORIS_ARMADILLO_LIBRARIES "${ARMADILLO_LIBRARY_DIRS}/libarmadillo.so"
    #    CACHE FILEPATH "Armadillo libraries." )
    #set(MORIS_ARMADILLO_DEFINITIONS "-DMORIS_USE_ARMA"
    #	CACHE INTERNAL "Moris preprocessor definitions for Armadillo.")
    set(MORIS_ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS})
    set(MORIS_ARMADILLO_DEFINITIONS "-DMORIS_USE_ARMA")
    set(MORIS_ARMADILLO_LIBRARIES "${ARMADILLO_LIBRARIES}")
    
    mark_as_advanced(Armadillo_DIR
        MORIS_ARMADILLO_INCLUDE_DIRS
        MORIS_ARMADILLO_LIBRARIES
        MORIS_ARMADILLO_DEFINITIONS )
    
    if(Armadillo_FOUND)
        #set(ARMADILLO_FOUND_ONCE TRUE CACHE INTERNAL "Armadillo was found.")
        set(ARMADILLO_FOUND_ONCE TRUE)
    endif()
    
    message(STATUS "ARMADILLO_INCLUDE_DIRS: ${ARMADILLO_INCLUDE_DIRS}")
	message(STATUS "ARMADILLO_LIBRARY_DIRS: ${ARMADILLO_LIBRARY_DIRS}")
	message(STATUS "ARMADILLO_LIBRARIES: ${ARMADILLO_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::armadillo)
	set(MORIS_ARMADILLO_TPLS
		${ACML_LAPACK_MKL}
		"arpack"
		)
	
	add_library(${MORIS}::armadillo STATIC IMPORTED GLOBAL)
	
	foreach(TPL ${MORIS_ARMADILLO_TPLS})
		include(${TPL}_new)
		list(APPEND MORIS_ARMADILLO_TPL_TARGETS ${MORIS}::${TPL})
	endforeach()
	
	set_target_properties(${MORIS}::armadillo PROPERTIES
		IMPORTED_LOCATION ${MORIS_ARMADILLO_LIBRARIES})
	target_link_libraries(${MORIS}::armadillo INTERFACE ${MORIS_ARMADILLO_TPL_TARGETS})
endif()

#add_definitions("-DMORIS_USE_ARMA")
#include_directories(${MORIS_ARMADILLO_INCLUDE_DIRS})
#link_directories(${MORIS_ARMADILLO_LIBRARY_DIRS})
#list(APPEND MORIS_ARMADILLO_EIGEN_LIBS "-larmadillo")
set(ARMADILLO_FOUND_ONCE TRUE)

