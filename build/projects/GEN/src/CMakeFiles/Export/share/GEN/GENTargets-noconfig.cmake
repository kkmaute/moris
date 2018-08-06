#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GEN-lib" for configuration ""
set_property(TARGET GEN-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(GEN-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/GEN/libGEN.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS GEN-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_GEN-lib "${_IMPORT_PREFIX}/lib/GEN/libGEN.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
