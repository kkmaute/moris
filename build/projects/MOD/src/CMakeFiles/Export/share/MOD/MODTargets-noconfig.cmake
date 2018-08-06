#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MOD-lib" for configuration ""
set_property(TARGET MOD-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(MOD-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/MOD/libMOD.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MOD-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_MOD-lib "${_IMPORT_PREFIX}/lib/MOD/libMOD.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
