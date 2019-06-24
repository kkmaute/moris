# Tools Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${TOL}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Debug.hpp
    cl_Enums.hpp
    cl_Geometry.hpp )

# List header dependencies
set(TOL_INCLUDES
    #${LINALG}
    #${ALG}
    )

# Include depedency directories
foreach(TOL_INCLUDE ${TOL_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TOL_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Library

# List source files
set(LIB_SOURCES
    cl_Debug.cpp
    cl_Geometry.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${ALG}-lib
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
add_library(${TOL}-lib STATIC ${LIB_SOURCES})
target_include_directories(${TOL}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries(${TOL}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${TOL}-lib PROPERTIES OUTPUT_NAME ${TOL})

foreach(TPL ${TOL_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${TOL}-lib
    	PUBLIC ${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${TOL}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${TOL}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${TOL}Config.cmake.in
    ${${TOL}_CONFIG_BUILD_DIR}/${TOL}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${TOL}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${TOL}_CONFIG_BUILD_DIR}/${TOL}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Installation

# Install *Config and *ConfigVersion files
install(FILES ${${TOL}_CONFIG_BUILD_DIR}/${TOL}Config.cmake 
    ${${TOL}_CONFIG_BUILD_DIR}/${TOL}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${TOL} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${TOL} )

# Install library
install(TARGETS ${TOL}-lib
    DESTINATION ${LIB_INSTALL_DIR}/${TOL}
    EXPORT ${TOL}Targets )

# Export library in to *Targets file
install(EXPORT ${TOL}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${TOL} )
