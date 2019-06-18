# MAP Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MAP}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    )

# List header dependencies
set(MAP_INCLUDES
    #${LINALG}
    #${COM}
    #${SOL}/${DLA}
    #${SOL}/${NLA}
    #${SOL}/${TSA}
    #${SOL}/${SOL_CORE}
    #${MTK}
    #${HMR}
    #${MTK}/${MAP}
    #${FEM}/${INT}
    #${FEM}/${MDL}
    #${FEM}/${MSI}
    )
    
# Include depedency directories
foreach(MAP_INCLUDE ${MAP_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MAP_INCLUDE}/src)
endforeach()
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAP Library

# List library source files
set( LIB_SOURCES
    cl_MTK_Mapper.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${CON}-lib
    ${MTK}-lib
    ${NLA}-lib
    ${INT}-lib
    ${TSA}-lib
    ${MDL}-lib
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE LIB_DEPENDENCIES)
list(REMOVE_DUPLICATES LIB_DEPENDENCIES)
list(REVERSE LIB_DEPENDENCIES)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MAP}-lib STATIC ${LIB_SOURCES})
target_include_directories(${MAP}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries(${MAP}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MAP}-lib PROPERTIES OUTPUT_NAME ${MAP})

foreach(TPL ${MAP_LIB_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MAP}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${MAP}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${MAP}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# map Config Files

# Build *Config and *ConfigVersion files
#include(CMakePackageConfigHelpers)

#configure_package_config_file(
#    ${MAP}Config.cmake.in
#    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}Config.cmake
#    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MAP}
#    PATH_VARS LIB_INSTALL_DIR )

#write_basic_package_version_file(
#    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}ConfigVersion.cmake
#    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Map Installation

# Install *Config and *ConfigVersion files
#install(FILES ${${GEN}_CONFIG_BUILD_DIR}/${GEN}Config.cmake 
#    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}ConfigVersion.cmake
#    DESTINATION ${CONFIG_INSTALL_DIR}/${GEN} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MTK}/${MAP} )

# Install library and executable
install(TARGETS ${MAP}-lib
    EXPORT ${MAP}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${MTK}/${MAP})

# Export library in to *Targets file
install(EXPORT ${MAP}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${MTK}/${MAP} )
