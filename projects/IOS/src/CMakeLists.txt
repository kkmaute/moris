# IOS --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${IOS}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    ios.hpp
    cl_Logger.hpp
    cl_Logger_Sink.hpp
    fn_iostream.hpp
    fn_Logger.hpp
    fn_to_stdio.hpp )

# List header dependencies
set(IOS_INCLUDES)

# Include depedency directories
foreach(IOS_INCLUDE ${IOS_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${IOS_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Library

# List source files
set(LIB_SOURCES
    fn_iostream.cpp )

# List library dependencies
set(LIB_DEPENDENCIES "")

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${IOS}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${IOS}-lib ${LIB_DEPENDENCIES})
set_target_properties(${IOS}-lib PROPERTIES OUTPUT_NAME ${IOS})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${IOS}Config.cmake.in
#     ${${IOS}_CONFIG_BUILD_DIR}/${IOS}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${IOS}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${IOS}_CONFIG_BUILD_DIR}/${IOS}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # IOS Installation
# 
# # Install *Config and *ConfigVersion files
# install(FILES ${${IOS}_CONFIG_BUILD_DIR}/${IOS}Config.cmake 
#     ${${IOS}_CONFIG_BUILD_DIR}/${IOS}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${IOS} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${IOS} )

# Install library
install(TARGETS ${IOS}-lib
    DESTINATION ${LIB_INSTALL_DIR}/${IOS}
    EXPORT ${IOS}Targets )

# Export library in to *Targets file
install(EXPORT ${IOS}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${IOS} )
