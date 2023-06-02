#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Optimization Source -----------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${OPT}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_OPT_Algorithm.hpp
    cl_OPT_Algorithm_GCMMA.hpp
    cl_OPT_Algorithm_LBFGS.hpp
    cl_OPT_Algorithm_SQP.hpp
    cl_OPT_Algorithm_Sweep.hpp
    cl_OPT_Manager.hpp
    cl_OPT_Problem.hpp
    cl_OPT_Criteria_Interface.hpp
    cl_OPT_Interface_Manager.hpp
    cl_OPT_Interface_User_Defined.hpp
    fn_OPT_create_problem.hpp
    fn_OPT_create_interface.hpp
    fn_OPT_create_algorithm.hpp)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Library

# List source files
set(LIB_SOURCES
    cl_OPT_Algorithm.cpp
    cl_OPT_Algorithm_GCMMA.cpp
    cl_OPT_Algorithm_LBFGS.cpp
    cl_OPT_Algorithm_SQP.cpp
    cl_OPT_Algorithm_Sweep.cpp
    cl_OPT_Manager.cpp
    cl_OPT_Problem.cpp
    cl_OPT_Criteria_Interface.cpp
    cl_OPT_Interface_Manager.cpp
    cl_OPT_Problem_User_Defined.cpp
    cl_OPT_Interface_User_Defined.cpp
    fn_OPT_create_problem.cpp
    fn_OPT_create_interface.cpp
    fn_OPT_create_algorithm.cpp)

# List library dependencies
set(LIB_DEPENDENCIES
    ${CNT}-lib
    ${LINALG}-lib
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    )
#message(INFO "base lib = ${MORIS_BASE_LIBS}\n\n\n\n\n")
# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${OPT}-lib STATIC ${LIB_SOURCES})
target_include_directories(${OPT}-lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${${OPT}_HEADER_INSTALL_DIR}> )
target_link_libraries(${OPT}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${OPT}-lib PROPERTIES OUTPUT_NAME ${OPT})
target_compile_definitions(${OPT}-lib PUBLIC "-DF77ADD_")

# Link third party libraries
foreach(TPL ${OPT_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${OPT}-lib
        PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${OPT}-lib 
        PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${OPT}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${OPT_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    string(APPEND OPT_CONFIG_TPL_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
    string(APPEND OPT_CONFIG_MORIS_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(OPT_COMPILE_DEFINITIONS ${OPT}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${OPT}Config.cmake.in
    ${${OPT}_CONFIG_BUILD_DIR}/${OPT}Config.cmake
    INSTALL_DESTINATION 
    ${${OPT}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${OPT}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${OPT}_CONFIG_BUILD_DIR}/${OPT}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Installation

# Install *Config and *ConfigVersion files
install(FILES ${${OPT}_CONFIG_BUILD_DIR}/${OPT}Config.cmake 
    ${${OPT}_CONFIG_BUILD_DIR}/${OPT}ConfigVersion.cmake
    DESTINATION ${${OPT}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${OPT}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${OPT}-lib
    DESTINATION ${${OPT}_LIB_INSTALL_DIR}
    EXPORT ${OPT}Targets )

# Export library into *Targets file
install(EXPORT ${OPT}Targets
    DESTINATION ${${OPT}_CONFIG_INSTALL_DIR} )

