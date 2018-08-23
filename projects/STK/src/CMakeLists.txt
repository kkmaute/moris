# STK Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${STK}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    mesh.hpp
    cl_STK_Implementation.hpp )

# List header dependencies
set(STK_INCLUDES
    ${LNA}
    ${ALG}
    ${MTK}
    ${COM}
    ${TOL}
    ${XTK})

foreach(STK_INCLUDE ${STK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${STK_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# STK Library

# List library source files
set(LIB_SOURCES
    cl_STK_Implementation.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_BOOST_LIBS} )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${STK}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${STK}-lib ${LIB_DEPENDENCIES})
set_target_properties(${STK}-lib PROPERTIES OUTPUT_NAME ${STK})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# STK Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# STK Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${STK}Config.cmake.in
    ${${STK}_CONFIG_BUILD_DIR}/${STK}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${STK}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${STK}_CONFIG_BUILD_DIR}/${STK}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# STK Installation

# Install *Config and *ConfigVersion files
install(FILES ${${STK}_CONFIG_BUILD_DIR}/${STK}Config.cmake 
    ${${STK}_CONFIG_BUILD_DIR}/${STK}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${STK} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${STK} )

# Install library and executable
install(TARGETS ${STK}-lib
    EXPORT ${STK}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${STK}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${STK}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${STK} )

# Export library in to *Targets file
install(EXPORT ${STK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${STK} )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Add Hierarchical Library

add_subdirectory(${HIA})
