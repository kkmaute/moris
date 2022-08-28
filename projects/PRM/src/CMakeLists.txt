#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Parameters Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${PRM}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
        fn_PRM_FEM_Parameters.hpp
        fn_PRM_GEN_Parameters.hpp
        fn_PRM_HMR_Parameters.hpp
        fn_PRM_MSI_Parameters.hpp
        fn_PRM_OPT_Parameters.hpp
        fn_PRM_SOL_Parameters.hpp
        fn_PRM_VIS_Parameters.hpp
        fn_PRM_XTK_Parameters.hpp
        fn_PRM_MORIS_GENERAL_Parameters.hpp
        fn_PRM_MIG_Parameters.hpp
        fn_PRM_WRK_Parameters.hpp
     )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Library

# List library source files
set(LIB_SOURCES
        fn_PRM_FEM_Parameters.cpp)

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${INT}-lib
     )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${PRM}-lib STATIC ${LIB_SOURCES})
target_include_directories(${PRM}-lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${${PRM}_HEADER_INSTALL_DIR}> )
target_link_libraries(${PRM}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${PRM}-lib PROPERTIES OUTPUT_NAME ${PRM})

# Link third party libraries
foreach(TPL ${PRM_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${PRM}-lib
        PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${PRM}-lib 
        PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${PRM}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${PRM_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    string(APPEND PRM_CONFIG_TPL_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
    string(APPEND PRM_CONFIG_MORIS_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be copied into config file
get_target_property(PRM_COMPILE_DEFINITIONS ${PRM}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${PRM}Config.cmake.in
    ${${PRM}_CONFIG_BUILD_DIR}/${PRM}Config.cmake
    INSTALL_DESTINATION
    ${CONFIG_INSTALL_DIR}/${PRM}
    PATH_VARS
    ${PRM}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${PRM}_CONFIG_BUILD_DIR}/${PRM}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Installation

# Install *Config and *ConfigVersion files
install(FILES ${${PRM}_CONFIG_BUILD_DIR}/${PRM}Config.cmake 
    ${${PRM}_CONFIG_BUILD_DIR}/${PRM}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${PRM} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${PRM} )

# Install library and executable
install(TARGETS ${PRM}-lib
    EXPORT ${PRM}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${PRM}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${PRM}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${PRM} )

# Export library in to *Targets file
install(EXPORT ${PRM}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${PRM} )

