#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

#  Model Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MDL}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_MDL_Model.hpp
    #cl_MDL_Mesh_Model_Helper.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Library

# List library source files
set(LIB_SOURCES
    cl_MDL_Model.cpp
    #cl_MDL_Mesh_Model_Helper.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${INT}-lib
	${MTK}-lib
    ${NLA}-lib
    ${TSA}-lib
    ${SDF}-lib
    ${DLA}-lib
    ${PRM}-lib
    ${VIS}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MDL}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${MDL}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${MDL}_HEADER_INSTALL_DIR}> )
target_link_libraries(${MDL}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MDL}-lib PROPERTIES OUTPUT_NAME ${MDL})

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)

foreach(TPL ${MDL_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${MDL}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MDL}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${MDL}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${MDL_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND MDL_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND MDL_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(MDL_COMPILE_DEFINITIONS ${MDL}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MDL}Config.cmake.in
    ${${MDL}_CONFIG_BUILD_DIR}/${MDL}Config.cmake
    INSTALL_DESTINATION 
    ${${MDL}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${MDL}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${MDL}_CONFIG_BUILD_DIR}/${MDL}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MDL}_CONFIG_BUILD_DIR}/${MDL}Config.cmake 
    ${${MDL}_CONFIG_BUILD_DIR}/${MDL}ConfigVersion.cmake
    DESTINATION ${${MDL}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${MDL}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${MDL}-lib
    DESTINATION ${${MDL}_LIB_INSTALL_DIR}
    EXPORT ${MDL}Targets )

# Export library into *Targets file
install(EXPORT ${MDL}Targets
    DESTINATION ${${MDL}_CONFIG_INSTALL_DIR} )

