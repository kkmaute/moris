# Tools Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${TOL}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    tools.hpp
    cl_Debug.hpp
    cl_Enums.hpp
    cl_FunctionFactory.hpp
    cl_Function.hpp
    cl_Function_Levelset.hpp
    cl_Geometry.hpp
    cl_Interpolation.hpp
    cl_Pairing.hpp )

# List header dependencies
set(TOL_INCLUDES
    ${LNA}
    ${IOS}
    ${ALG})

# Include depedency directories
foreach(TOL_INCLUDE ${TOL_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TOL_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Library

# List source files
set(LIB_SOURCES
	cl_Interpolation.cpp
	cl_Function_Levelset.cpp
    cl_Function.cpp
    cl_FunctionFactory.cpp
    cl_Pairing.cpp
    cl_Debug.cpp
    cl_Geometry.cpp )

# List library dependencies
set(LIB_DEPENDENCIES "")

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${TOL}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${TOL}-lib ${LIB_DEPENDENCIES})
set_target_properties(${TOL}-lib PROPERTIES OUTPUT_NAME ${TOL})

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
