#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# IOS --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    cl_Tracer_Enums.hpp
	cl_GlobalClock.hpp
	cl_Tracer.hpp
    cl_Logger.hpp
    cl_Query.hpp
	cl_XML_Parser.hpp
	fn_Parsing_Tools.hpp
    cl_Library_Enums.hpp
	cl_Library_Factory.hpp
    cl_Library_IO.hpp
    cl_Library_IO_Meshgen.hpp
    cl_Library_IO_Standard.hpp
	fn_Logger.hpp
    fn_stringify.hpp
    fn_stringify_matrix.hpp
	fn_to_stdio.hpp
	HDF5_Tools.hpp
	cl_Ascii.hpp
	IO_Tools.hpp
	ios.hpp
	Log_Constants.hpp
	)

# List library source files
set(LIB_SOURCES
    cl_Library_Enums.cpp
	cl_Library_IO.cpp
    cl_Library_IO_Meshgen.cpp
    cl_Library_IO_Standard.cpp
	cl_GlobalClock.cpp
	cl_Logger.cpp
    cl_Query.cpp
    cl_Query_Table.cpp
    cl_Query_Tree.cpp
    cl_Tracer_Enums.cpp
    cl_XML_Parser.cpp
    cl_Ascii.cpp
    )


# List library dependencies
set(LIB_DEPENDENCIES
    ${COR}-lib
    ${COM}-lib
    ${PRM}-lib
	${GEN_CORE}-lib
    ${MORIS_BASE_LIBS}
     )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# ------------------------------------------------------------------------
# Add git revision
# https://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake

include(${MORIS_CMAKE_DIR}/utilities/GetGitRevisionDescription.cmake)
get_git_head_revision(MORIS_GIT_BRANCH MORIS_GIT_HASH)
message( STATUS "done with GetGitRevisionDescription " ${MORIS_GIT_BRANCH} "  " ${MORIS_GIT_HASH})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cl_Git_info.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/cl_Git_info.cpp" @ONLY)

list(APPEND LIB_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/cl_Git_info.cpp")
#list(APPEND HEADERS moris_git_info.hpp)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Target
#add_library(${IOS}-lib INTERFACE)
add_library(${IOS}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${IOS}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${IOS}_HEADER_INSTALL_DIR}> )

target_link_libraries(${IOS}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${IOS}-lib PROPERTIES OUTPUT_NAME ${IOS})

foreach(TPL ${IOS_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${IOS}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${IOS}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})

    target_compile_definitions(${IOS}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${IOS_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND IOS_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND IOS_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${IOS}Config.cmake.in
    ${${IOS}_CONFIG_BUILD_DIR}/${IOS}Config.cmake
    INSTALL_DESTINATION 
    ${${IOS}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${IOS}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${IOS}_CONFIG_BUILD_DIR}/${IOS}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Installation

# Install *Config and *ConfigVersion files
install(FILES ${${IOS}_CONFIG_BUILD_DIR}/${IOS}Config.cmake 
    ${${IOS}_CONFIG_BUILD_DIR}/${IOS}ConfigVersion.cmake
    DESTINATION ${${IOS}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${IOS}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${IOS}-lib
    DESTINATION ${${IOS}_LIB_INSTALL_DIR}
    EXPORT ${IOS}Targets )

# Export library in to *Targets file
install(EXPORT ${IOS}Targets
    DESTINATION ${${IOS}_CONFIG_INSTALL_DIR} )

