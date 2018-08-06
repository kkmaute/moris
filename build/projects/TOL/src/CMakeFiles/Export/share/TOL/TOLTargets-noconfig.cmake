#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TOL-lib" for configuration ""
set_property(TARGET TOL-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(TOL-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/TOL/libTOL.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS TOL-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_TOL-lib "${_IMPORT_PREFIX}/lib/TOL/libTOL.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
