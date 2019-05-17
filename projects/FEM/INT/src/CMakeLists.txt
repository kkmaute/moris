# Integration and Interpolation Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${INT}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_FEM_Element_Factory.hpp
    cl_FEM_ElementProxy.hpp
    cl_FEM_NodeProxy.hpp
    cl_FEM_Cluster.hpp
    cl_FEM_Element.hpp
    cl_FEM_Element_Bulk.hpp
    cl_FEM_Set.hpp
    cl_FEM_Element_Sideset.hpp
    cl_FEM_Element_Time_Sideset.hpp
    cl_FEM_IWG_Factory.hpp
    cl_FEM_IWG.hpp
    cl_FEM_IWG_L2.hpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Sideset.hpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Neumann.hpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost_Left.hpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost_Right.hpp
    cl_FEM_IWG_Helmholtz_Bulk.hpp
    cl_FEM_IWG_Helmholtz_Bulk2.hpp
    cl_FEM_IWG_Helmholtz_Interface.hpp
    cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp
    cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp
    cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp
    cl_FEM_IWG_LSNormal_Bulk.hpp
    cl_FEM_IWG_Olsson_CLS_Bulk.hpp
    cl_FEM_IWG_Olsson_CLS_Interface.hpp
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
    cl_FEM_Integration_Coeffs_Tet_1.hpp
    cl_FEM_Integration_Coeffs_Tet_4.hpp
    cl_FEM_Integration_Coeffs_Tet_5.hpp
    cl_FEM_Integration_Coeffs_Tet_11.hpp
    cl_FEM_Integration_Coeffs_Tet_15.hpp
    cl_FEM_Integration_Rule.hpp
    cl_FEM_Integrator.hpp
    cl_FEM_Integrator_Test_Polynomial.hpp
    cl_FEM_Interpolation_Function_Base.hpp
    cl_FEM_Interpolation_Function_Factory.hpp
    cl_FEM_Interpolation_Function.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar2.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp
    cl_FEM_Interpolation_Function_Lagrange_Bar4.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex64.hpp
    cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad16.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad4.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad8.hpp
    cl_FEM_Interpolation_Function_Lagrange_Quad9.hpp
    cl_FEM_Interpolation_Function_Lagrange_Tri3.hpp
    cl_FEM_Interpolation_Function_Lagrange_Tri6.hpp
    cl_FEM_Interpolation_Function_Lagrange_Tri10.hpp
    cl_FEM_Interpolation_Function_Lagrange_Tet4.hpp
    cl_FEM_Interpolation_Function_Lagrange_Tet10.hpp
    cl_FEM_Interpolation_Function_Lagrange_Tet20.hpp
    cl_FEM_Interpolation_Function_Constant_Bar2.hpp
    cl_FEM_Interpolation_Rule.hpp
    cl_FEM_Node_Base.hpp
    cl_FEM_Node.hpp
    cl_FEM_Field_Interpolator.hpp
    )

# List header dependencies
set(INT_INCLUDES
    ${SOL}/${DLA}
    ${SOL}/${NLA}
    ${LINALG}
    ${MTK}
    ${COM}
    ${FEM}/${MSI}
    ${FEM}/${MDL} )

# Include depedency directories
foreach(INT_INCLUDE ${INT_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${INT_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Library

# List library source files
set(LIB_SOURCES
    cl_FEM_Cluster.cpp
    cl_FEM_Element_Factory.cpp
    cl_FEM_Element_Bulk.cpp
    cl_FEM_Set.cpp
    cl_FEM_Element_Sideset.cpp
    cl_FEM_Element_Time_Sideset.cpp
    cl_FEM_IWG_Factory.cpp
    cl_FEM_IWG_L2.cpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.cpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Sideset.cpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.cpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Neumann.cpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost_Left.cpp
    cl_FEM_IWG_Isotropic_Spatial_Diffusion_Ghost_Right.cpp
    cl_FEM_IWG_Helmholtz_Bulk.cpp
    cl_FEM_IWG_Helmholtz_Bulk2.cpp
    cl_FEM_IWG_Helmholtz_Interface.cpp
    cl_FEM_IWG_Hamilton_Jacobi_Bulk.cpp
    cl_FEM_IWG_Hamilton_Jacobi_Bulk2.cpp
    cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.cpp
    cl_FEM_IWG_LSNormal_Bulk.cpp
    cl_FEM_IWG_Olsson_CLS_Bulk.cpp
    cl_FEM_IWG_Olsson_CLS_Interface.cpp
    cl_FEM_IWG_Test.cpp
    cl_FEM_Interpolation_Function_Factory.cpp
    cl_FEM_Interpolation_Rule.cpp
    cl_FEM_Geometry_Interpolator.cpp
    cl_FEM_Field_Interpolator.cpp
    cl_FEM_Integration_Rule.cpp
    cl_FEM_Integrator.cpp
    cl_FEM_Integrator_Test_Polynomial.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MSI}-lib
    ${MDL}-lib)

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
