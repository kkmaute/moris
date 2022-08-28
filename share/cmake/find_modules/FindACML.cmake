#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# ACML Find Module --------------------------------------------------------
# -------------------------------------------------------------------------

if (ACML_LIBRARIES)
  set(ACML_FIND_QUIETLY TRUE)
endif (ACML_LIBRARIES)

set(ACML_ENV_VARS
  $ENV{ACMLDIR}
  $ENV{ACML_DIR}
  $ENV{acml_DIR}
  $ENV{ACML_ROOT}
  $ENV{acml_ROOT}
  $ENV{ACML_PATH}
  $ENV{acml_PATH} )

find_library(ACML_LIBRARIES
  NAMES
  acml_mp acml_mv
  PATHS
  ${ACML_ENV_VARS}
  /usr/lib/acml
  PATH_SUFFIXES
  lib
)

find_file(ACML_LIBRARIES
  NAMES
  libacml_mp.so
  PATHS
  /usr/lib
  ${ACML_ENV_VARS}
  ${LIB_INSTALL_DIR}
  PATH_SUFFIXES
  lib
)

if(NOT ACML_LIBRARIES)
    message(STATUS "Multi-threaded library not found, looking for single-threaded")
    find_library(ACML_LIBRARIES
        NAMES
        acml acml_mv
        PATHS
        ${ACML_ENV_VARS}
        ${LIB_INSTALL_DIR}
        PATH_SUFFIXES
        lib
        )
    find_file(ACML_LIBRARIES
        libacml.so libacml_mv.so
        PATHS
        /usr/lib
        ${ACML_ENV_VARS}
        ${LIB_INSTALL_DIR}
        PATH_SUFFIXES
        lib
        )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACML DEFAULT_MSG ACML_LIBRARIES)

mark_as_advanced(ACML_LIBRARIES)

_import_libraries(ACML_LIBRARY_TARGETS ${ACML_LIBRARIES})

add_library(ACML::acml INTERFACE IMPORTED GLOBAL)
target_link_libraries(ACML::acml INTERFACE ${ACML_LIBRARY_TARGETS})


