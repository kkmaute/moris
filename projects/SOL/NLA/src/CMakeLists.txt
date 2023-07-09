#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# NonLinear Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${NLA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	#cl_NLA_Arc_Length.hpp
    cl_NLA_Convergence.hpp
    cl_NLA_Solver_Relaxation.hpp
    cl_NLA_Solver_Pseudo_Time_Control.hpp
    cl_NLA_Solver_Load_Control.hpp
	cl_NLA_Newton_Solver.hpp
	cl_NLA_NLBGS.hpp
	cl_NLA_Nonlinear_Algorithm.hpp
	cl_NLA_Nonlinear_Problem.hpp
	cl_NLA_Nonlinear_Solver_Factory.hpp
	cl_NLA_Nonlinear_Solver.hpp
	)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Library

# List library source files
set(LIB_SOURCES
    #cl_NLA_Arc_Length.cpp   
    cl_NLA_Newton_Solver.cpp
    cl_NLA_Nonlinear_Solver_Factory.cpp
    cl_NLA_Nonlinear_Algorithm.cpp
    cl_NLA_Nonlinear_Solver.cpp
    cl_NLA_Nonlinear_Problem.cpp
    cl_NLA_Convergence.cpp
    cl_NLA_Solver_Relaxation.cpp
    cl_NLA_Solver_Pseudo_Time_Control.cpp
    cl_NLA_Solver_Load_Control.cpp
    cl_NLA_NLBGS.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${DLA}-lib
    ${COM}-lib
    ${MSI}-lib
    ${SOL_CORE}-lib
    ${ENM}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${NLA}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${NLA}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${NLA}_HEADER_INSTALL_DIR}> )
target_link_libraries(${NLA}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${NLA}-lib PROPERTIES OUTPUT_NAME ${NLA})

foreach(TPL ${NLA_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${NLA}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${NLA}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${NLA}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${NLA_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND NLA_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND NLA_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(NLA_COMPILE_DEFINITIONS ${NLA}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${NLA}Config.cmake.in
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}Config.cmake
    INSTALL_DESTINATION 
    ${${NLA}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${NLA}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Installation

# Install *Config and *ConfigVersion files
install(FILES ${${NLA}_CONFIG_BUILD_DIR}/${NLA}Config.cmake 
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}ConfigVersion.cmake
    DESTINATION ${${NLA}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${NLA}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${NLA}-lib
    DESTINATION ${${NLA}_LIB_INSTALL_DIR}
    EXPORT ${NLA}Targets )

# Export library into *Targets file
install(EXPORT ${NLA}Targets
    DESTINATION ${${NLA}_CONFIG_INSTALL_DIR} )

