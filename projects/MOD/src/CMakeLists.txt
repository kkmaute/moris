# Model Source ------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MOD}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    model.hpp
    cl_Bspline.hpp
    cl_EqnManager.hpp
    cl_Model_Enums.hpp )

# List header dependencies
set(MOD_INCLUDES
    ${ALG}
    ${LNA})

# Include depedency directories
foreach(MOD_INCLUDE ${MOD_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MOD_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Library

# List source files
set(LIB_SOURCES
    cl_Bspline.cpp
    cl_EqnManager.cpp )

# List library dependencies
set(LIB_DEPENDENCIES "")

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MOD}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${MOD}-lib ${LIB_DEPENDENCIES})
set_target_properties(${MOD}-lib PROPERTIES OUTPUT_NAME ${MOD})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MOD}Config.cmake.in
    ${${MOD}_CONFIG_BUILD_DIR}/${MOD}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MOD}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${MOD}_CONFIG_BUILD_DIR}/${MOD}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MOD}_CONFIG_BUILD_DIR}/${MOD}Config.cmake 
    ${${MOD}_CONFIG_BUILD_DIR}/${MOD}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${MOD} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MOD} )

# Install library
install(TARGETS ${MOD}-lib
    DESTINATION ${LIB_INSTALL_DIR}/${MOD}
    EXPORT ${MOD}Targets )

# Export library in to *Targets file
install(EXPORT ${MOD}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${MOD} )
