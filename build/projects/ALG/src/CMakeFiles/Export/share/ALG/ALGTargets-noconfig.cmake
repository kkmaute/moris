#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ALG-lib" for configuration ""
set_property(TARGET ALG-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(ALG-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/ALG/libALG.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS ALG-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_ALG-lib "${_IMPORT_PREFIX}/lib/ALG/libALG.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
