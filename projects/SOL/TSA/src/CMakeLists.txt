#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${TSA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_TSA_Time_Solver_Algorithm.hpp
    cl_TSA_Monolithic_Time_Solver.hpp
    cl_TSA_Staggered_Time_Solver.hpp
    cl_TSA_Time_Solver_Factory.hpp
    cl_TSA_Time_Solver.hpp
 )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Library

# List library source files
set(LIB_SOURCES
    cl_TSA_Time_Solver_Algorithm.cpp
    cl_TSA_Monolithic_Time_Solver.cpp
    cl_TSA_Staggered_Time_Solver.cpp
    cl_TSA_Time_Solver_Factory.cpp
    cl_TSA_Time_Solver.cpp
 )

# List library dependencies
set(LIB_DEPENDENCIES
    ${SOL_CORE}-lib
    ${DLA}-lib
    ${NLA}-lib
    ${MSI}-lib
    ${ENM}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${TSA}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${TSA}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${TSA}_HEADER_INSTALL_DIR}> )
target_link_libraries(${TSA}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${TSA}-lib PROPERTIES OUTPUT_NAME ${TSA})

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)

foreach(TPL ${TSA_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${TSA}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${TSA}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${TSA}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${TSA_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND TSA_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND TSA_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(TSA_COMPILE_DEFINITIONS ${TSA}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${TSA}Config.cmake.in
    ${${TSA}_CONFIG_BUILD_DIR}/${TSA}Config.cmake
    INSTALL_DESTINATION 
    ${${TSA}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${TSA}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${TSA}_CONFIG_BUILD_DIR}/${TSA}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Installation

# Install *Config and *ConfigVersion files
install(FILES ${${TSA}_CONFIG_BUILD_DIR}/${TSA}Config.cmake 
    ${${TSA}_CONFIG_BUILD_DIR}/${TSA}ConfigVersion.cmake
    DESTINATION ${${TSA}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${TSA}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${TSA}-lib
    DESTINATION ${${TSA}_LIB_INSTALL_DIR}
    EXPORT ${TSA}Targets )

# Export library into *Targets file
install(EXPORT ${TSA}Targets
    DESTINATION ${${TSA}_CONFIG_INSTALL_DIR} )

