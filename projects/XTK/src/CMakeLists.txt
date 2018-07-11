# XTK Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${XTK}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Entity_Tracker.hpp
    cl_Mesh_XTK.hpp
    cl_Model.hpp
    cl_XTK_Enums.hpp )

# List header dependencies
set(XTK_INCLUDES
    ${LNA}
    ${IOS}
    ${MTK}
    ${TOL})

# Include depedency directories
foreach(XTK_INCLUDE ${XTK_INCLUDES})
    include_directories(../../${XTK_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Library

# List source files
set(LIB_SOURCES
    cl_XTK_Enums.hpp
    cl_Entity_Tracker.cpp
    cl_Mesh_XTK.cpp
    cl_Model.cpp )

# List library dependencies
set(LIB_DEPENDENCIES "")

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${XTK}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${XTK}-lib ${LIB_DEPENDENCIES})
set_target_properties(${XTK}-lib PROPERTIES OUTPUT_NAME ${XTK})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${XTK}Config.cmake.in
#     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${XTK}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # XTK Installation
# 
# # Install *Config and *ConfigVersion files
# install(FILES ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake 
#     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${XTK} )

# Install library
install(TARGETS ${XTK}-lib
    DESTINATION ${LIB_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )

# Export library in to *Targets file
install(EXPORT ${XTK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )
