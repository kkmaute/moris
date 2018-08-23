# Integration and Interpolation Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${INT}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_FEM_Element.hpp
    cl_FEM_Enums.hpp
    cl_FEM_Geometry_Interpolator.hpp
    cl_FEM_Integration_Coeffs_Bar_1.hpp
    cl_FEM_Integration_Coeffs_Bar_2.hpp
    cl_FEM_Integration_Coeffs_Bar_3.hpp
    cl_FEM_Integration_Coeffs_Bar_4.hpp
    cl_FEM_Integration_Coeffs_Bar_5.hpp
    cl_FEM_Integration_Coeffs_Base.hpp
    cl_FEM_Integration_Coeffs_Hex_2x2x2.hpp
    cl_FEM_Integration_Coeffs_Hex_3x3x3.hpp
    cl_FEM_Integration_Coeffs_Hex_4x4x4.hpp
    cl_FEM_Integration_Coeffs_Hex_5x5x5.hpp
    cl_FEM_Integration_Coeffs.hpp
    cl_FEM_Integration_Coeffs_Quad_2x2.hpp
    cl_FEM_Integration_Coeffs_Quad_3x3.hpp
    cl_FEM_Integration_Coeffs_Quad_4x4.hpp
    cl_FEM_Integration_Coeffs_Quad_5x5.hpp
    cl_FEM_Integration_Coeffs_Tri_1.hpp
    cl_FEM_Integration_Coeffs_Tri_3.hpp
    cl_FEM_Integration_Coeffs_Tri_6.hpp
    cl_FEM_Integration_Rule.hpp
    cl_FEM_Integrator.hpp
    cl_FEM_Integrator_Test_Polynomial.hpp
    cl_FEM_Interpolation_Function_Base.hpp
    cl_FEM_Interpolation_Function_Factory.hpp
    cl_FEM_Interpolation_Function.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar2.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex20.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex64.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad16.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad4.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad8.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad9.hpp
    cl_FEM_Interpolation_Matrix.hpp
    cl_FEM_Interpolation_Rule.hpp
    cl_FEM_Interpolator.hpp
    cl_FEM_Node.hpp
    )

# List header dependencies
set(INT_INCLUDES
    ${LNA}
    ${IOS} #LNA
    ${MTK}
    ${COM}
    ${FEM}/${MSI} )

# Include depedency directories
foreach(INT_INCLUDE ${INT_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${INT_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Library

# List library source files
set(LIB_SOURCES
    cl_FEM_Element.cpp
    cl_FEM_Geometry_Interpolator.cpp
    cl_FEM_Integration_Rule.cpp
    cl_FEM_Integrator.cpp
    cl_FEM_Integrator_Test_Polynomial.cpp
    cl_FEM_Interpolation_Matrix.cpp
    cl_FEM_Interpolation_Function_Factory.cpp
    cl_FEM_Interpolation_Rule.cpp
    cl_FEM_Interpolator.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_LDLIBS})

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${INT}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${INT}-lib ${LIB_DEPENDENCIES})
set_target_properties(${INT}-lib PROPERTIES OUTPUT_NAME ${INT})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${INT}Config.cmake.in
#     ${${INT}_CONFIG_BUILD_DIR}/${INT}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${INT}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${INT}_CONFIG_BUILD_DIR}/${INT}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Installation

# # Install *Config and *ConfigVersion files
# install(FILES ${${INT}_CONFIG_BUILD_DIR}/${INT}Config.cmake 
#     ${${INT}_CONFIG_BUILD_DIR}/${INT}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${INT} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${INT} )

# Install library and executable
install(TARGETS ${INT}-lib
    EXPORT ${INT}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${INT}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${INT}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${INT} )

# Export library in to *Targets file
install(EXPORT ${INT}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${INT} )
