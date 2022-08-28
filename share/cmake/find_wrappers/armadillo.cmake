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
    
    set(ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS}
        CACHE PATH "Armadillo include directories." )
    set(ARMADILLO_LIBRARY_DIRS ${ARMADILLO_LIBRARY_DIRS}
        CACHE PATH "Armadillo library directories." )
    set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARIES}
        CACHE FILEPATH "Armadillo libraries." )
    
    mark_as_advanced(Armadillo_DIR
        ARMADILLO_INCLUDE_DIRS
        ARMADILLO_LIBRARY_DIRS
        ARMADILLO_LIBRARIES )
    
    if(Armadillo_FOUND)
        set(ARMADILLO_FOUND_ONCE TRUE CACHE INTERNAL "Armadillo was found.")
    endif()
    
    message(STATUS "ARMADILLO_INCLUDE_DIRS: ${ARMADILLO_INCLUDE_DIRS}")
	message(STATUS "ARMADILLO_LIBRARY_DIRS: ${ARMADILLO_LIBRARY_DIRS}")
	message(STATUS "ARMADILLO_LIBRARIES: ${ARMADILLO_LIBRARIES}")
endif()

add_definitions("-DMORIS_USE_ARMA")
include_directories(${ARMADILLO_INCLUDE_DIRS})
link_directories(${ARMADILLO_LIBRARY_DIRS})
#list(APPEND MORIS_ARMADILLO_EIGEN_LIBS "-larmadillo")
list(APPEND MORIS_ARMADILLO_EIGEN_LIBS ${MORIS_ARMADILLO_LIBRARIES})

