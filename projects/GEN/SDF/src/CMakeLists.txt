#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

#  SDF Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${SDF}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_SDF_Arguments.hpp
    cl_SDF_Cell.hpp
    cl_SDF_Core.hpp
    cl_SDF_Facet_Vertex.hpp
    cl_SDF_Facet.hpp
    cl_SDF_Field.hpp
    cl_SDF_Generator.hpp
    cl_SDF_Line.hpp
    cl_SDF_Mesh.hpp
    cl_SDF_Object.hpp
    cl_SDF_Parameters.hpp
    cl_SDF_State.hpp
    cl_SDF_STK.hpp
    cl_SDF_Triangle.hpp
    cl_SDF_Vertex.hpp
    SDF_Tools.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SDF Library

# List library source files
set(LIB_SOURCES
    cl_SDF_Arguments.cpp
    cl_SDF_Cell.cpp
    cl_SDF_Core.cpp
    cl_SDF_Facet_Vertex.cpp
    cl_SDF_Facet.cpp
    cl_SDF_Field.cpp
    cl_SDF_Generator.cpp
    cl_SDF_Line.cpp
    cl_SDF_Mesh.cpp
    cl_SDF_Object.cpp
    cl_SDF_Parameters.cpp
    cl_SDF_STK.cpp
    cl_SDF_Triangle.cpp
    cl_SDF_Vertex.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${COM}-lib
	${CNT}-lib
	${MTK}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${SDF}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${SDF}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${SDF}_HEADER_INSTALL_DIR}> )
target_link_libraries(${SDF}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${SDF}-lib PROPERTIES OUTPUT_NAME ${SDF})

foreach(TPL ${SDF_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${SDF}-lib
        PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${SDF}-lib
        PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${SDF}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SDF Executable

# List executable source files
set(EXE_SOURCES
    main.cpp )

# List executable dependencies
set(EXE_DEPENDENCIES
	${CHR}-lib
    ${SDF}-lib
    ${MTK}-lib
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${SDF}-exe ${EXE_SOURCES})
target_link_libraries(${SDF}-exe PRIVATE ${EXE_DEPENDENCIES})
set_target_properties(${SDF}-exe PROPERTIES OUTPUT_NAME sdf )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SDF Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${SDF_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND SDF_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND SDF_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(SDF_COMPILE_DEFINITIONS ${SDF}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${SDF}Config.cmake.in
    ${${SDF}_CONFIG_BUILD_DIR}/${SDF}Config.cmake
    INSTALL_DESTINATION 
    ${${SDF}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${SDF}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${SDF}_CONFIG_BUILD_DIR}/${SDF}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SDF Installation

# Install *Config and *ConfigVersion files
install(FILES ${${SDF}_CONFIG_BUILD_DIR}/${SDF}Config.cmake 
    ${${SDF}_CONFIG_BUILD_DIR}/${SDF}ConfigVersion.cmake
    DESTINATION ${${SDF}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${SDF}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${SDF}-lib
    DESTINATION ${${SDF}_LIB_INSTALL_DIR}
    EXPORT ${SDF}Targets )

# Install executable
install(TARGETS ${SDF}-exe
	DESTINATION ${${SDF}_EXE_INSTALL_DIR}
	EXPORT ${SDF}Targets )

# Export library into *Targets file
install(EXPORT ${SDF}Targets
    DESTINATION ${${SDF}_CONFIG_INSTALL_DIR} )

