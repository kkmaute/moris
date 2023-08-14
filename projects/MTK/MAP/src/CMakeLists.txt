#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MAP Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MAP}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	cl_MTK_Mapper_Node.hpp
	cl_MTK_Mapper.hpp
    )
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAP Library

# List library source files
set( LIB_SOURCES
    cl_MTK_Mapper.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${CNT}-lib
    ${MTK}-lib
    ${NLA}-lib
    ${INT}-lib
    ${TSA}-lib
    ${MDL}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MAP}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${MAP}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${MAP}_HEADER_INSTALL_DIR}> )
target_link_libraries(${MAP}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MAP}-lib PROPERTIES OUTPUT_NAME ${MAP})

foreach(TPL ${MAP_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${MAP}-lib
    	PUBLIC ${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MAP}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${MAP}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# map Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${MAP_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND MAP_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND MAP_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(MAP_COMPILE_DEFINITIONS ${MAP}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MAP}Config.cmake.in
    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}Config.cmake
    INSTALL_DESTINATION 
    ${${MAP}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${MAP}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Map Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MAP}_CONFIG_BUILD_DIR}/${MAP}Config.cmake 
    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}ConfigVersion.cmake
    DESTINATION ${${MAP}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${MAP}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${MAP}-lib
    DESTINATION ${${MAP}_LIB_INSTALL_DIR}
    EXPORT ${MAP}Targets )

# Export library into *Targets file
install(EXPORT ${MAP}Targets
    DESTINATION ${${MAP}_CONFIG_INSTALL_DIR} )

