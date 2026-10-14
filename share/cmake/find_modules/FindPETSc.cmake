# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

if(PETSC_INCLUDE_DIRS)
    set(PETSC_FIND_QUIETLY TRUE)
endif(PETSC_INCLUDE_DIRS)

set(PETSC_HEADER "petsc.h")
set(PETSC_LIB_OPT "petsc")
set(PETSC_LIB_DBG "${PETSC_LIB_OPT}")

set(PETSC_ENV_VARS
    $ENV{PETSCDIR}
    $ENV{PETSC_DIR}
    $ENV{PETSc_DIR}
    $ENV{PETSC_ROOT}
    $ENV{PETSc_ROOT}
    $ENV{PETSC_PATH}
    $ENV{PETSc_PATH}
)

find_path(PETSC_DIR NAMES include/${PETSC_HEADER} PATHS ${PETSC_ENV_VARS} /usr/lib/petscdir)

find_path(PETSC_INCLUDE_DIR NAMES ${PETSC_HEADER}
          PATHS ${PETSC_DIR}/include ${PETSC_INCLUDE_PATH} /usr/local/include/petsc
                /usr/local/include /usr/include/petsc /usr/include
)

set(PETSC_ARCH $ENV{PETSC_ARCH} CACHE PATH "The PETSC architecture.")
find_path(PETSC_ARCH_INCLUDE_DIR NAMES petscconf.h PATHS ${PETSC_DIR}/${PETSC_ARCH}/include)

if(EXISTS ${PETSC_DIR}/${PETSC_ARCH}/lib)
    set(PETSC_LIBRARY_DIR ${PETSC_DIR}/${PETSC_ARCH}/lib
        CACHE PATH "The directory where the petsc library resides."
    )
else()
    set(PETSC_LIBRARY_DIR NOTFOUND CACHE PATH "The directory where the petsc library resides.")
endif()
find_library(
    PETSC_LIBRARY_DEBUG NAMES ${PETSC_LIB_DBG} HINTS ${PETSC_LIBRARY_DIR}
                                                     "${CMAKE_SOURCE_DIR}/MacOS/Libs/gurobi40"
)
find_library(
    PETSC_LIBRARY_RELEASE NAMES ${PETSC_LIB_OPT} HINTS ${PETSC_LIBRARY_DIR}
                                                       "${CMAKE_SOURCE_DIR}/MacOS/Libs/gurobi40"
)
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
    PETSc DEFAULT_MSG PETSC_LIBRARY_DEBUG PETSC_LIBRARY_RELEASE PETSC_INCLUDE_DIR
)

if(PETSC_FOUND)
    set(PETSC_LIBRARIES "${PETSC_LIBRARY}")
    set(PETSC_INCLUDE_DIRS "${PETSC_INCLUDE_DIR}")
    if(EXISTS ${PETSC_ARCH_INCLUDE_DIR})
        list(APPEND PETSC_INCLUDE_DIRS "${PETSC_ARCH_INCLUDE_DIR}")
    endif()
    set(PETSC_LIBRARY debug ${PETSC_LIBRARY_DEBUG} optimized ${PETSC_LIBRARY_RELEASE})
endif(PETSC_FOUND)

mark_as_advanced(
    PETSC_DIR
    PETSC_INCLUDE_DIR
    PETSC_ARCH
    PETSC_ARCH_INCLUDE_DIR
    PETSC_LIBRARY_DIR
    PETSC_LIBRARY_DEBUG
    PETSC_LIBRARY_RELEASE
)
