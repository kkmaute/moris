# Visualization -------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${VIS}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_VIS_Output_Manager.hpp
    cl_VIS_Visualization_Mesh.hpp
    cl_VIS_Vertex_Visualization.hpp
    cl_VIS_Cell_Cluster_Visualization.hpp
    cl_VIS_Factory.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Visualization Library

# List library source files
set(LIB_SOURCES
    cl_VIS_Output_Manager.cpp
    cl_VIS_Cell_Cluster_Visualization.cpp
    cl_VIS_Factory.cpp
 )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${COM}-lib
    ${CON}-lib
    ${HMR}-lib
    ${MTK}-lib
    #${MDL}-lib
   )

set(GEN_INCLUDES
	"field"
	"field/geometry"
	"field/property"
	"pdv"
    )

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
endforeach()

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${VIS}-lib STATIC ${LIB_SOURCES})
target_include_directories(${VIS}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${VIS}_HEADER_INSTALL_DIR}> )
target_link_libraries(${VIS}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${VIS}-lib PROPERTIES OUTPUT_NAME ${VIS})

foreach(TPL ${VIS_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${VIS}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${VIS}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${VIS}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Visualization Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${VIS_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND VIS_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND VIS_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(VIS_COMPILE_DEFINITIONS ${VIS}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${VIS}Config.cmake.in
    ${${VIS}_CONFIG_BUILD_DIR}/${VIS}Config.cmake
    INSTALL_DESTINATION 
    ${${VIS}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${VIS}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${VIS}_CONFIG_BUILD_DIR}/${VIS}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Visualization Installation

# Install *Config and *ConfigVersion files
install(FILES ${${VIS}_CONFIG_BUILD_DIR}/${VIS}Config.cmake 
    ${${VIS}_CONFIG_BUILD_DIR}/${VIS}ConfigVersion.cmake
    DESTINATION ${${VIS}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${VIS}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${VIS}-lib
    DESTINATION ${${VIS}_LIB_INSTALL_DIR}
    EXPORT ${VIS}Targets )

# Export library into *Targets file
install(EXPORT ${VIS}Targets
    DESTINATION ${${VIS}_CONFIG_INSTALL_DIR} )
