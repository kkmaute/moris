#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# OPENBLAS Find Module --------------------------------------------------------
# -------------------------------------------------------------------------

if (OPENBLAS_LIBRARIES)
  set(OPENBLAS_FIND_QUIETLY TRUE)
endif (OPENBLAS_LIBRARIES)

set(OPENBLAS_ENV_VARS
  $ENV{OPENBLASDIR}
  $ENV{OPENBLAS_DIR}
  $ENV{OPENBLAS_DIR}
  $ENV{OPENBLAS_ROOT}
  $ENV{openblas_ROOT}
  $ENV{OPENBLAS_PATH}
  $ENV{openblas_PATH} )

find_library(OPENBLAS_LIBRARIES
  NAMES
  openblas
  PATHS
  ${OPENBLAS_ENV_VARS}
  /usr/lib64/
  PATH_SUFFIXES
  lib64
)

find_file(OPENBLAS_LIBRARIES
  NAMES
  libopenblas.so
  PATHS
  /usr/lib64
  ${OPENBLAS_ENV_VARS}
  ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib64
)

if(NOT OPENBLAS_LIBRARIES)
    message(STATUS "Multi-threaded library not found, looking for single-threaded")
    find_library(OPENBLAS_LIBRARIES
        NAMES
        openblas
        PATHS
        ${OPENBLAS_ENV_VARS}
        ${LIB_INSTALL_DIR}
        PATH_SUFFIXES
        lib64
        )
    find_file(OPENBLAS_LIBRARIES
        libopenblas.so
        PATHS
        /usr/lib64
        ${OPENBLAS_ENV_VARS}
        ${LIB_INSTALL_DIR}
        PATH_SUFFIXES
        lib64
        )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OPENBLAS DEFAULT_MSG OPENBLAS_LIBRARIES)

mark_as_advanced(OPENBLAS_LIBRARIES)

_import_libraries(OPENBLAS_LIBRARY_TARGETS ${OPENBLAS_LIBRARIES})

add_library(OPENBLAS::openblas INTERFACE IMPORTED GLOBAL)
target_link_libraries(OPENBLAS::openblas INTERFACE ${OPENBLAS_LIBRARY_TARGETS})

add_library(OPENBLAS::all_libs INTERFACE IMPORTED)
target_link_libraries(OPENBLAS::all_libs INTERFACE ${OPENBLAS_LIBRARY_TARGETS})
