#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

IF (PETSC_INCLUDE_DIRS)
  SET(PETSC_FIND_QUIETLY TRUE)
ENDIF (PETSC_INCLUDE_DIRS)

SET (PETSC_HEADER "petsc.h")
SET (PETSC_LIB_OPT "petsc")
SET (PETSC_LIB_DBG "${PETSC_LIB_OPT}")

set(PETSC_ENV_VARS
  $ENV{PETSCDIR}
  $ENV{PETSC_DIR}
  $ENV{PETSc_DIR}
  $ENV{PETSC_ROOT}
  $ENV{PETSc_ROOT}
  $ENV{PETSC_PATH}
  $ENV{PETSc_PATH} )

FIND_PATH(PETSC_DIR NAMES include/${PETSC_HEADER}
    PATHS
    ${PETSC_ENV_VARS}
    /usr/lib/petscdir
)

FIND_PATH(PETSC_INCLUDE_DIR NAMES ${PETSC_HEADER}
    PATHS
    ${PETSC_DIR}/include
    ${PETSC_INCLUDE_PATH}
    /usr/local/include/petsc
    /usr/local/include
    /usr/include/petsc
    /usr/include
)

SET(PETSC_ARCH $ENV{PETSC_ARCH} CACHE PATH "The PETSC architecture.")
FIND_PATH(PETSC_ARCH_INCLUDE_DIR NAMES petscconf.h
    PATHS
    ${PETSC_DIR}/${PETSC_ARCH}/include
)

IF (EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib)
    SET(PETSC_LIBRARY_DIR ${PETSC_DIR}/${PETSC_ARCH}/lib
        CACHE PATH "The directory where the petsc library resides.")
ELSE()
    SET(PETSC_LIBRARY_DIR NOTFOUND
        CACHE PATH "The directory where the petsc library resides.")
ENDIF()
FIND_LIBRARY( PETSC_LIBRARY_DEBUG
    NAMES ${PETSC_LIB_DBG}
    HINTS
    ${PETSC_LIBRARY_DIR}
    "${CMAKE_SOURCE_DIR}/MacOS/Libs/gurobi40"
    )
FIND_LIBRARY( PETSC_LIBRARY_RELEASE
    NAMES ${PETSC_LIB_OPT}
    HINTS
    ${PETSC_LIBRARY_DIR}
    "${CMAKE_SOURCE_DIR}/MacOS/Libs/gurobi40"
    )
INCLUDE (FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(Petsc DEFAULT_MSG PETSC_LIBRARY_DEBUG PETSC_LIBRARY_RELEASE PETSC_INCLUDE_DIR)

IF (PETSC_FOUND)
    SET(PETSC_LIBRARIES "${PETSC_LIBRARY}")
    SET(PETSC_INCLUDE_DIRS "${PETSC_INCLUDE_DIR}")
    IF (EXISTS ${PETSC_ARCH_INCLUDE_DIR})
        LIST(APPEND PETSC_INCLUDE_DIRS "${PETSC_ARCH_INCLUDE_DIR}")
    ENDIF()
    SET( PETSC_LIBRARY
        debug ${PETSC_LIBRARY_DEBUG}
        optimized ${PETSC_LIBRARY_RELEASE} )
ENDIF (PETSC_FOUND)

mark_as_advanced(PETSC_DIR
    PETSC_INCLUDE_DIR
    PETSC_ARCH
    PETSC_ARCH_INCLUDE_DIR
    PETSC_LIBRARY_DIR
    PETSC_LIBRARY_DEBUG
    PETSC_LIBRARY_RELEASE )

