# Optimization Source -----------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${OPT}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Opt_Alg.hpp
    cl_Opt_Alg_API.hpp
    cl_Opt_Alg_GCMMA.hpp
    cl_Opt_Alg_LBFGS.hpp
    cl_Opt_Alg_SQP.hpp
    cl_Opt_Alg_Sweep.hpp
    cl_Opt_Manager.hpp
    cl_Opt_Prob.hpp
    Opt_Input.hpp )

# List header dependencies
set(OPT_INCLUDES
    ${LNA}
    ${IOS})

# Include depedency directories
foreach(OPT_INCLUDE ${OPT_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${OPT_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Library

# List source files
set(LIB_SOURCES
    cl_Opt_Alg.cpp
    cl_Opt_Alg_API.cpp
    cl_Opt_Alg_GCMMA.cpp
    cl_Opt_Alg_LBFGS.cpp
    cl_Opt_Alg_SQP.cpp
    cl_Opt_Alg_Sweep.cpp
    cl_Opt_Manager.cpp
    cl_Opt_Prob.cpp
    Opt_Input.cpp )

# List library dependencies
set(LIB_DEPENDENCIES "")

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${OPT}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${OPT}-lib ${LIB_DEPENDENCIES})
set_target_properties(${OPT}-lib PROPERTIES OUTPUT_NAME ${OPT})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Optimization Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${OPT}Config.cmake.in
#     ${${OPT}_CONFIG_BUILD_DIR}/${OPT}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${OPT}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${OPT}_CONFIG_BUILD_DIR}/${OPT}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# # Optimization Installation
# 
# # Install *Config and *ConfigVersion files
# install(FILES ${${OPT}_CONFIG_BUILD_DIR}/${OPT}Config.cmake 
#     ${${OPT}_CONFIG_BUILD_DIR}/${OPT}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${OPT} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${OPT} )

# Install library
install(TARGETS ${OPT}-lib
    DESTINATION ${LIB_INSTALL_DIR}/${OPT}
    EXPORT ${OPT}Targets )

# Export library in to *Targets file
install(EXPORT ${OPT}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${OPT} )
