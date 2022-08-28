#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Communications Source ---------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${COM}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Communication_Enums.hpp
    cl_Communication_Manager.hpp
    cl_Communication_Tools.hpp
    cl_MPI_Tools.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Library

# List library source files
set(LIB_SOURCES
    cl_Communication_Manager.cpp
    cl_Communication_Tools.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${COM}-lib STATIC ${LIB_SOURCES})
target_include_directories(${COM}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${COM}_HEADER_INSTALL_DIR}> )
target_link_libraries(${COM}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${COM}-lib PROPERTIES OUTPUT_NAME ${COM})

# Link third party libraries
foreach(TPL ${COM_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${COM}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${COM}-lib 
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${COM}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${COM_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND COM_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND COM_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
    string(APPEND COM_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(COM_COMPILE_DEFINITIONS ${COM}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${COM}Config.cmake.in
    ${${COM}_CONFIG_BUILD_DIR}/${COM}Config.cmake
    INSTALL_DESTINATION
    ${CONFIG_INSTALL_DIR}/${COM}
    PATH_VARS
    ${COM}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${COM}_CONFIG_BUILD_DIR}/${COM}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Installation

# Install *Config and *ConfigVersion files
install(FILES ${${COM}_CONFIG_BUILD_DIR}/${COM}Config.cmake 
    ${${COM}_CONFIG_BUILD_DIR}/${COM}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${COM} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${COM} )

# Install library and executable
install(TARGETS ${COM}-lib
    EXPORT ${COM}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${COM}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${COM}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${COM} )

# Export library in to *Targets file
install(EXPORT ${COM}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${COM} )

