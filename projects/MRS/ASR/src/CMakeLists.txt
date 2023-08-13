#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Assert Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    assert.hpp
    fn_assert.hpp
    fn_log_assert.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Library
set(LIB_DEPENDENCIES
	${IOS}-lib
	${EXC}-lib
	${COR}-lib
	)

add_library(${ASR}-lib INTERFACE)
target_include_directories(${ASR}-lib INTERFACE 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${ASR}_HEADER_INSTALL_DIR}> )
target_link_libraries(${ASR}-lib INTERFACE ${LIB_DEPENDENCIES})

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND ASR_CONFIG_TPL_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND ASR_CONFIG_TPL_DEPENDENCIES
		"find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND ASR_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${ASR}Config.cmake.in
    ${${ASR}_CONFIG_BUILD_DIR}/${ASR}Config.cmake
    INSTALL_DESTINATION 
    ${${ASR}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${ASR}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${ASR}_CONFIG_BUILD_DIR}/${ASR}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert Installation

# Install *Config and *ConfigVersion files
install(FILES ${${ASR}_CONFIG_BUILD_DIR}/${ASR}Config.cmake 
    ${${ASR}_CONFIG_BUILD_DIR}/${ASR}ConfigVersion.cmake
    DESTINATION ${${ASR}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${ASR}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${ASR}-lib
    DESTINATION ${${ASR}_LIB_INSTALL_DIR}
    EXPORT ${ASR}Targets )

# Export library in to *Targets file
install(EXPORT ${ASR}Targets
    DESTINATION ${${ASR}_CONFIG_INSTALL_DIR} )

