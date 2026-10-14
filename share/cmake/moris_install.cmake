# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Export the public MORIS configuration and shared definitions.
# Component packages retain their existing names and installation directories.
# Target to catch general moris properties
add_library(${MORIS}_general INTERFACE)
target_compile_definitions(${MORIS}_general INTERFACE ${MORIS_DEFINITIONS})
target_compile_features(${MORIS}_general INTERFACE cxx_std_17)

# Needed for the following macros
include(CMakePackageConfigHelpers)

# Configure *Config file
set(MORIS_HAS_RUNTIME_TARGET TRUE)
# Only library packages belong in the component list; GUI, mains, and examples
# have build directories but do not install component configurations.
get_property(MORIS_INSTALLED_COMPONENTS GLOBAL PROPERTY MORIS_INSTALLED_COMPONENTS)
get_property(MORIS_INSTALLED_COMPONENT_DIRS GLOBAL PROPERTY MORIS_INSTALLED_COMPONENT_DIRS)
configure_package_config_file(
    ${MORIS_CMAKE_DIR}/${MORIS}Config.cmake.in ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS} PATH_VARS CONFIG_INSTALL_DIR
)

# Create *ConfigVersion file
write_basic_package_version_file(
    ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}ConfigVersion.cmake COMPATIBILITY AnyNewerVersion
)

# Install package configuration files.
install(FILES ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}Config.cmake
              ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}ConfigVersion.cmake
        DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS}
)

# Install general moris target
install(TARGETS ${MORIS}_general DESTINATION ${LIB_INSTALL_DIR}/${MORIS} EXPORT ${MORIS}Targets)

install(EXPORT ${MORIS}Targets DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS})

install(FILES "include/catch.hpp" DESTINATION ${HEADER_INSTALL_DIR}/${MORIS})

# Install moris_macros.cmake file
install(FILES ${CMAKE_SOURCE_DIR}/share/cmake/utilities/moris_macros.cmake
        DESTINATION ${CONFIG_INSTALL_DIR}/utilities
)
