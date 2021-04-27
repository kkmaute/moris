# GEN Core Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN_CORE}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
        cl_GEN_Pdv_Enums.hpp
		GEN_typedefs.hpp
 )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Library

# List library source files
set(LIB_SOURCES
		cl_GEN_Pdv_Enums.cpp
 )

# List library dependencies
set(LIB_DEPENDENCIES

    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${GEN_CORE}-lib STATIC ${LIB_SOURCES})
target_include_directories(${GEN_CORE}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${GEN_CORE}_HEADER_INSTALL_DIR}> )
target_link_libraries(${GEN_CORE}-lib ${LIB_DEPENDENCIES})
set_target_properties(${GEN_CORE}-lib PROPERTIES OUTPUT_NAME ${GEN_CORE})

foreach(TPL ${GEN_CORE_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${GEN_CORE}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${GEN_CORE}-lib
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${GEN_CORE}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${GEN_CORE_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND GEN_CORE_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	# string(APPEND GEN_CORE_CONFIG_MORIS_DEPENDENCIES
	# 	"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )


        string(APPEND XTK_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )

endforeach()

# Properties to be coppied into config file
get_target_property(GEN_CORE_COMPILE_DEFINITIONS ${GEN_CORE}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${GEN_CORE}Config.cmake.in
    ${${GEN_CORE}_CONFIG_BUILD_DIR}/${GEN_CORE}Config.cmake
    INSTALL_DESTINATION 
    ${${GEN_CORE}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${GEN_CORE}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${GEN_CORE}_CONFIG_BUILD_DIR}/${GEN_CORE}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Installation

# Install *Config and *ConfigVersion files
install(FILES ${${GEN_CORE}_CONFIG_BUILD_DIR}/${GEN_CORE}Config.cmake 
    ${${GEN_CORE}_CONFIG_BUILD_DIR}/${GEN_CORE}ConfigVersion.cmake
    DESTINATION ${${GEN_CORE}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${GEN_CORE}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${GEN_CORE}-lib
    DESTINATION ${${GEN_CORE}_LIB_INSTALL_DIR}
    EXPORT ${GEN_CORE}Targets )

# Export library into *Targets file
install(EXPORT ${GEN_CORE}Targets
    DESTINATION ${${GEN_CORE}_CONFIG_INSTALL_DIR} )
