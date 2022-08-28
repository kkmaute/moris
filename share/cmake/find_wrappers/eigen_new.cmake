#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# Eigen libraries ---------------------------------------------------------
# -------------------------------------------------------------------------
if(NOT EIGEN_FOUND_ONCE)
    set(EIGEN_ENV_VARS
        $ENV{EIGENDIR}
        $ENV{EIGEN_DIR}
        $ENV{Eigen_DIR}
        $ENV{Eigen3_DIR}
        $ENV{EIGEN_ROOT}
        $ENV{Eigen_ROOT}
        $ENV{Eigen3_ROOT}
        $ENV{EIGEN_PATH}
        $ENV{Eigen_PATH}
        $ENV{Eigen3_PATH})

    find_package(Eigen
        REQUIRED
        NAMES Eigen3
        HINTS ${EIGEN_ENV_VARS})
    
    set(MORIS_EIGEN_INCLUDE_DIRS "${EIGEN3_ROOT_DIR}/include"
        CACHE PATH "Eigen include directories." )
    set(MORIS_EIGEN_TARGETS "Eigen3::Eigen"
    	CACHE INTERNAL "Eigen library.")
    set(MORIS_EIGEN_TARGET_FILE "${EIGEN3_ROOT_DIR}/share/Eigen3Targets.cmake"
        CACHE FILEPATH "Eigen targets file.")
    set(MORIS_EIGEN_DEFINITIONS "-DMORIS_USE_EIGEN"
    	CACHE INTERNAL "Moris preprocessor definitions for Eigen.")
    
    mark_as_advanced(MORIS_EIGEN_INCLUDE_DIRS
    	MORIS_EIGEN_TARGETS
    	MORIS_EIGEN_TARGET_FILE
    	MORIS_EIGEN_DEFINITIONS
        Eigen_DIR
        )
    
    if(EIGEN3_FOUND)
        set(EIGEN_FOUND_ONCE TRUE CACHE INTERNAL "Eigen was found.")
    endif()
    
    message(STATUS "EIGEN3_ROOT_DIR: ${EIGEN3_ROOT_DIR}")
else()
    include(${MORIS_EIGEN_TARGET_FILE})
endif()

if(NOT TARGET ${MORIS}::eigen)
	add_library(${MORIS}::eigen INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::eigen INTERFACE ${MORIS_EIGEN_TARGETS} )
endif()

#add_definitions("-DMORIS_USE_EIGEN")
#include_directories("${MORIS_EIGEN_INCLUDE_DIRS}")
#list(APPEND MORIS_ARMADILLO_EIGEN_LIBS "Eigen3::Eigen")

