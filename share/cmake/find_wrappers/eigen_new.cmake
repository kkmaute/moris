# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Eigen libraries
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
        $ENV{Eigen3_PATH}
    )

    find_package(Eigen REQUIRED NAMES Eigen3 HINTS ${EIGEN_ENV_VARS})

    set(MORIS_EIGEN_TARGETS "Eigen3::Eigen" CACHE INTERNAL "Eigen library.")
    set(MORIS_EIGEN_TARGET_FILE "${EIGEN3_ROOT_DIR}/share/eigen3/cmake/Eigen3Targets.cmake"
        CACHE FILEPATH "Eigen targets file."
    )
    set(MORIS_EIGEN_DEFINITIONS "-DMORIS_USE_EIGEN"
        CACHE INTERNAL "Moris preprocessor definitions for Eigen."
    )

    include(${MORIS_EIGEN_TARGET_FILE})

    mark_as_advanced(
        MORIS_EIGEN_INCLUDE_DIRS MORIS_EIGEN_TARGETS MORIS_EIGEN_TARGET_FILE
        MORIS_EIGEN_DEFINITIONS Eigen_DIR
    )

    if(EIGEN3_FOUND)
        set(EIGEN_FOUND_ONCE TRUE CACHE INTERNAL "Eigen was found.")
    endif()

    message(STATUS "EIGEN3_ROOT_DIR:          ${EIGEN3_ROOT_DIR}")
    message(STATUS "MORIS_EIGEN_TARGET_FILE:  ${MORIS_EIGEN_TARGET_FILE}")
    message(STATUS "MORIS_EIGEN_DEFINITIONS:  ${MORIS_EIGEN_DEFINITIONS}")
else()
    include(${MORIS_EIGEN_TARGET_FILE})
endif()

# Input plugins consume this variable without linking Eigen3::Eigen. Use the
# package's actual include layout and repair values cached by older wrappers.
get_target_property(MORIS_EIGEN_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
set(MORIS_EIGEN_INCLUDE_DIRS "${MORIS_EIGEN_INCLUDE_DIRS}" CACHE PATH "Eigen include directories."
                                                                 FORCE
)

if(NOT TARGET ${MORIS}::eigen)
    add_library(${MORIS}::eigen INTERFACE IMPORTED GLOBAL)
    target_link_libraries(${MORIS}::eigen INTERFACE ${MORIS_EIGEN_TARGETS})
endif()

#add_definitions("-DMORIS_USE_EIGEN")
#include_directories("${MORIS_EIGEN_INCLUDE_DIRS}")
#list(APPEND MORIS_ARMADILLO_EIGEN_LIBS "Eigen3::Eigen")
