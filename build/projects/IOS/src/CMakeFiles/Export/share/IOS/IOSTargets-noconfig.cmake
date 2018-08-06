#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "IOS-lib" for configuration ""
set_property(TARGET IOS-lib APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(IOS-lib PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/IOS/libIOS.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS IOS-lib )
list(APPEND _IMPORT_CHECK_FILES_FOR_IOS-lib "${_IMPORT_PREFIX}/lib/IOS/libIOS.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
