#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "DLA-lib" for configuration ""
set_property(TARGET DLA-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(DLA-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/DLA/libDLA.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS DLA-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_DLA-lib "${_IMPORT_PREFIX}/lib/DLA/libDLA.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
