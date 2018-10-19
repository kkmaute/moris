#  Model Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${SDF}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_SDF_Triangle_Vertex.hpp
    cl_SDF_Triangle.hpp
    cl_SDF_Object.hpp
    )

# List header dependencies
set(SDF_INCLUDES
    ${LINALG}
    ${TOL}
    ${IOS}
    ${MTK}
    ${COM})

# Include depedency directories
foreach(SDF_INCLUDE ${SDF_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${SDF_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Library

# List library source files
set(LIB_SOURCES
    cl_SDF_Arguments.cpp
    cl_SDF_Triangle_Vertex.cpp
    cl_SDF_Triangle.cpp
    cl_SDF_Object.cpp
    cl_SDF_Cell.cpp
    cl_SDF_Vertex.cpp
    cl_SDF_Core.cpp
    cl_SDF_Data.cpp
    cl_SDF_Mesh.cpp
    cl_SDF_Generator.cpp
    cl_SDF_Parameters.cpp
    cl_SDF_STK.cpp
    cl_SDF_Field.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_BASE_LIBS})

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${SDF}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${SDF}-lib ${LIB_DEPENDENCIES})
set_target_properties(${SDF}-lib PROPERTIES OUTPUT_NAME ${SDF})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${SDF}Config.cmake.in
#     ${${SDF}_CONFIG_BUILD_DIR}/${SDF}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${SDF}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${SDF}_CONFIG_BUILD_DIR}/${SDF}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Installation

# # Install *Config and *ConfigVersion files
# install(FILES ${${SDF}_CONFIG_BUILD_DIR}/${SDF}Config.cmake 
#     ${${SDF}_CONFIG_BUILD_DIR}/${SDF}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${SDF} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${SDF} )

# Install library and executable
install(TARGETS ${SDF}-lib
    EXPORT ${SDF}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${SDF}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${SDF}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${SDF} )

# Export library in to *Targets file
install(EXPORT ${SDF}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${SDF} )
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# project Executable

# List executable source files
set(EXE_SOURCES
    main.cpp )

# List executable dependencies
set(EXE_DEPENDENCIES
#    ${GEN}-lib
    ${SDF}-lib
    ${MTK}-lib
#    ${TOL}-lib why?
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )

# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${SDF}-exe ${EXE_SOURCES})
target_link_libraries(${SDF}-exe ${EXE_DEPENDENCIES})
set_target_properties(${SDF}-exe PROPERTIES OUTPUT_NAME sdf )
