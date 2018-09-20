# Algorithms Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${ALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    algorithms.hpp
    fn_equal_to.hpp
    fn_equal_to.tpp
    fn_hash_value.hpp
    fn_num_digits.hpp )

# List header dependencies
set(ALG_INCLUDES
    ${LINALG})

# Include depedency directories
foreach(ALG_INCLUDE ${ALG_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${ALG_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Library

# List source files
set(LIB_SOURCES
    fn_equal_to.cpp )

# List library dependencies
set(LIB_DEPENDENCIES "")

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${ALG}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${ALG}-lib ${LIB_DEPENDENCIES})
set_target_properties(${ALG}-lib PROPERTIES OUTPUT_NAME ${ALG})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${ALG}Config.cmake.in
    ${${ALG}_CONFIG_BUILD_DIR}/${ALG}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${ALG}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${ALG}_CONFIG_BUILD_DIR}/${ALG}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Installation

# Install *Config and *ConfigVersion files
install(FILES ${${ALG}_CONFIG_BUILD_DIR}/${ALG}Config.cmake 
    ${${ALG}_CONFIG_BUILD_DIR}/${ALG}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${ALG} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${ALG} )

# Install library
install(TARGETS ${ALG}-lib
    DESTINATION ${LIB_INSTALL_DIR}/${ALG}
    EXPORT ${ALG}Targets )

# Export library in to *Targets file
install(EXPORT ${ALG}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${ALG} )
