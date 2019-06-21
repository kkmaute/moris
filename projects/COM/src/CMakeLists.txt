# Communications Source ---------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${COM}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Communication_Enums.hpp
    cl_Communication_Manager.hpp
    cl_Communication_Tools.hpp )

# List header dependencies
set(COM_INCLUDES
    #${LINALG}
    "" )

foreach(COM_INCLUDE ${COM_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${COM_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Library

# List library source files
set(LIB_SOURCES
    cl_Communication_Manager.cpp
    cl_Communication_Tools.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE LIB_DEPENDENCIES)
list(REMOVE_DUPLICATES LIB_DEPENDENCIES)
list(REVERSE LIB_DEPENDENCIES)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${COM}-lib STATIC ${LIB_SOURCES})
target_include_directories(${COM}-lib 
	PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(${COM}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${COM}-lib PROPERTIES OUTPUT_NAME ${COM})

# Link third party libraries
foreach(TPL ${COM_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${COM}-lib
    	PUBLIC ${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${COM}-lib 
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${COM}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${COM}Config.cmake.in
    ${${COM}_CONFIG_BUILD_DIR}/${COM}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${COM}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${COM}_CONFIG_BUILD_DIR}/${COM}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Communications Installation

# Install *Config and *ConfigVersion files
install(FILES ${${COM}_CONFIG_BUILD_DIR}/${COM}Config.cmake 
    ${${COM}_CONFIG_BUILD_DIR}/${COM}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${COM} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${COM} )

# Install library and executable
install(TARGETS ${COM}-lib
    EXPORT ${COM}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${COM}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${COM}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${COM} )

# Export library in to *Targets file
install(EXPORT ${COM}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${COM} )
