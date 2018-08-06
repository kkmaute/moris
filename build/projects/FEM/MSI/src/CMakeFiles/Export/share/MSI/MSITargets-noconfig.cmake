#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "MSI-lib" for configuration ""
set_property(TARGET MSI-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(MSI-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/MSI/libMSI.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS MSI-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_MSI-lib "${_IMPORT_PREFIX}/lib/MSI/libMSI.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
