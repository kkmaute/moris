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
  $ENV{AMDBLIS_DIR}
  $ENV{AMDLIBFLAME_DIR}
)

find_library(AMBLIS_LIBRARIES
  NAMES
  blis
  PATHS
  ${ACML_ENV_VARS}
  PATH_SUFFIXES
  lib
  lib64
)

find_library(AMDFLAME_LIBRARIES
  NAMES
  flame 
  PATHS
  ${ACML_ENV_VARS}
  PATH_SUFFIXES
  lib
  lib64
)
find_library(AMDAOCLDTL_LIBRARIES
  NAMES
  aocldtl 
  PATHS
  ${ACML_ENV_VARS}
  PATH_SUFFIXES
  lib
  lib64
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ACML DEFAULT_MSG AMBLIS_LIBRARIES)
find_package_handle_standard_args(ACML DEFAULT_MSG AMDFLAME_LIBRARIES)
find_package_handle_standard_args(ACML DEFAULT_MSG AMDAOCLDTL_LIBRARIES)

SET(ACML_LIBRARIES "${AMBLIS_LIBRARIES}")
LIST(APPEND ACML_LIBRARIES "${AMDFLAME_LIBRARIES}")
LIST(APPEND ACML_LIBRARIES "${AMDAOCLDTL_LIBRARIES}")

mark_as_advanced(AMBLIS_LIBRARIES)
mark_as_advanced(AMDFLAME_LIBRARIES)
mark_as_advanced(AMDAOCLDTL_LIBRARIES)

_import_libraries(AMBLIS_LIBRARY_TARGETS ${AMBLIS_LIBRARIES})
_import_libraries(AMDFLAME_LIBRARY_TARGETS ${AMDFLAME_LIBRARIES})
_import_libraries(AMDAOCLDTL_LIBRARY_TARGETS ${AMDAOCLDTL_LIBRARIES})

if (NOT TARGET ACML::acml )
    add_library(ACML::acml INTERFACE IMPORTED GLOBAL)
    target_link_libraries(ACML::acml INTERFACE 
         ${AMBLIS_LIBRARY_TARGETS} 
         ${AMDFLAME_LIBRARY_TARGETS} 
         ${AMDAOCLDTL_LIBRARY_TARGETS} )
endif()

if (NOT TARGET ACML::all_libs )
    add_library(ACML::all_libs INTERFACE IMPORTED)
    target_link_libraries(ACML::all_libs INTERFACE 
         ${AMBLIS_LIBRARY_TARGETS} 
         ${AMDFLAME_LIBRARY_TARGETS} 
         ${AMDAOCLDTL_LIBRARY_TARGETS} )
endif()

