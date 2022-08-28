#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Core Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    core.hpp
    banner.hpp
    common.hpp
    typedefs.hpp )
    
configure_file( paths.hpp.in ${CMAKE_BINARY_DIR}/generated/paths.hpp )
include_directories( ${CMAKE_BINARY_DIR}/generated/ ) # Make sure it can be included...
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Core Library
add_library(${COR}-lib INTERFACE)
target_include_directories(${COR}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/>
	$<INSTALL_INTERFACE:${${COR}_HEADER_INSTALL_DIR}> )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${COR_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	list(APPEND COR_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${COR}Config.cmake.in
    ${${COR}_CONFIG_BUILD_DIR}/${COR}Config.cmake
    INSTALL_DESTINATION 
    ${${COR}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${COR}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${COR}_CONFIG_BUILD_DIR}/${COR}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Installation

# Install *Config and *ConfigVersion files
install(FILES ${${COR}_CONFIG_BUILD_DIR}/${COR}Config.cmake 
    ${${COR}_CONFIG_BUILD_DIR}/${COR}ConfigVersion.cmake
    DESTINATION ${${COR}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${COR}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${COR}-lib
    DESTINATION ${${COR}_LIB_INSTALL_DIR}
    EXPORT ${COR}Targets )

# Export library in to *Targets file
install(EXPORT ${COR}Targets
    DESTINATION ${${COR}_CONFIG_INSTALL_DIR} )
