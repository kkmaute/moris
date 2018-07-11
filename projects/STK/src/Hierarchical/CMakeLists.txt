# STK Hierarchical Source -------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${HIA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Base_Mesh_Edge.hpp
    cl_Base_Mesh_Element.hpp
    cl_Base_Mesh_Face.hpp
    cl_Hierarchical_Mesh_Basis.hpp
    cl_Hierarchical_Mesh_Edge.hpp
    cl_Hierarchical_Mesh_Element.hpp
    cl_Hierarchical_Mesh_Face.hpp
    cl_Hierarchical_Mesh_Filter.hpp
    cl_Hierarchical_Mesh_Holes.hpp
    cl_Hierarchical_Mesh.hpp
    cl_Hierarchical_Mesh_Input.hpp
    cl_Hierarchical_Mesh_Main.hpp
    cl_Hierarchical_Mesh_MPI.hpp
    cl_Hierarchical_Mesh_MTK.hpp
    cl_Hierarchical_Mesh_Output.hpp
    cl_Hierarchical_Mesh_Refinement.hpp
    cl_Hierarchical_Mesh_Sets.hpp
    cl_Hierarchical_Mesh_Solver_Input.hpp
    cl_Hierarchical_Mesh_TMatrix.hpp
    cl_HMR_Implementation.hpp
    cl_Lagrange_Basis.hpp
    cl_Lagrange_Edge.hpp
    cl_Lagrange_Element.hpp
    cl_Lagrange_Face.hpp
    cl_Lagrange_Filter.hpp
    cl_Mesh_Element_Bucket.hpp )

# List header dependencies
set(HIA_INCLUDES
    ${MOD}
    ${DLA}
    ${STK})

foreach(HIA_INCLUDE ${HIA_INCLUDES})
    include_directories(../../../${HIA_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HIA Library

# List library source files
set(LIB_SOURCES
    cl_Base_Mesh_Edge.cpp
    cl_Base_Mesh_Element.cpp
    cl_Base_Mesh_Face.cpp
    cl_Hierarchical_Mesh_Basis.cpp
    cl_Hierarchical_Mesh.cpp
    cl_Hierarchical_Mesh_Edge.cpp
    cl_Hierarchical_Mesh_Element.cpp
    cl_Hierarchical_Mesh_Face.cpp
    cl_Hierarchical_Mesh_Filter.cpp
    cl_Hierarchical_Mesh_Holes.cpp
    cl_Hierarchical_Mesh_Input.cpp
    cl_Hierarchical_Mesh_Main.cpp
    cl_Hierarchical_Mesh_MPI.cpp
    cl_Hierarchical_Mesh_MTK.cpp
    cl_Hierarchical_Mesh_Output.cpp
    cl_Hierarchical_Mesh_Refinement.cpp
    cl_Hierarchical_Mesh_Sets.cpp
    cl_Hierarchical_Mesh_Solver_Input.cpp
    cl_Hierarchical_Mesh_TMatrix.cpp
    cl_HMR_Implementation.cpp
    cl_Lagrange_Basis.cpp
    cl_Lagrange_Edge.cpp
    cl_Lagrange_Element.cpp
    cl_Lagrange_Face.cpp
    cl_Lagrange_Filter.cpp
    cl_Mesh_Element_Bucket.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    "" )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${HIA}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${HIA}-lib ${LIB_DEPENDENCIES})
set_target_properties(${HIA}-lib PROPERTIES OUTPUT_NAME ${HIA})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HIA Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HIA Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${HIA}Config.cmake.in
#     ${${HIA}_CONFIG_BUILD_DIR}/${HIA}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${STK}/${HIA}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${HIA}_CONFIG_BUILD_DIR}/${HIA}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HIA Installation

# # Install *Config and *ConfigVersion files
# install(FILES ${${HIA}_CONFIG_BUILD_DIR}/${HIA}Config.cmake 
#     ${${HIA}_CONFIG_BUILD_DIR}/${HIA}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${STK}/${HIA} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${STK}/${HIA} )

# Install library and executable
install(TARGETS ${HIA}-lib
    EXPORT ${HIA}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${STK}/${HIA}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${STK}/${HIA}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${STK}/${HIA} )

# Export library in to *Targets file
install(EXPORT ${HIA}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${STK}/${HIA} )
