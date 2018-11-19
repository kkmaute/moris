# HMR Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${HMR}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_HMR_Arguments.hpp
    cl_HMR_Background_Edge.hpp
    cl_HMR_Background_Element_Base.hpp
    cl_HMR_Background_Element.hpp
    cl_HMR_Background_Facet.hpp
    cl_HMR_Background_Mesh_1D.hpp
    cl_HMR_Background_Mesh_2D.hpp
    cl_HMR_Background_Mesh_3D.hpp
    cl_HMR_Background_Mesh_Base.hpp
    cl_HMR_Background_Mesh.hpp
    cl_HMR_Basis.hpp
    cl_HMR_BSpline_Element_Hex27.hpp
    cl_HMR_BSpline_Element_Hex64.hpp
    cl_HMR_BSpline_Element_Hex8.hpp
    cl_HMR_BSpline_Element.hpp
    cl_HMR_BSpline_Element_Quad16.hpp
    cl_HMR_BSpline_Element_Quad4.hpp
    cl_HMR_BSpline_Element_Quad9.hpp
    cl_HMR_BSpline.hpp
    cl_HMR_BSpline_Mesh_Base.hpp
    cl_HMR_BSpline_Mesh.hpp
    cl_HMR_Database.hpp
    cl_HMR_Domain.hpp
    cl_HMR_Edge.hpp
    cl_HMR_Element.hpp
    cl_HMR_Facet.hpp
    cl_HMR_Factory.hpp
    cl_HMR_Field.hpp
    cl_HMR_Fields.hpp
    cl_HMR_File.hpp
    cl_HMR.hpp
    cl_HMR_Knot.hpp
    cl_HMR_Lagrange_Edge2.hpp
    cl_HMR_Lagrange_Edge3.hpp
    cl_HMR_Lagrange_Edge4.hpp
    cl_HMR_Lagrange_Edge.hpp
    cl_HMR_Lagrange_Element_Hex27.hpp
    cl_HMR_Lagrange_Element_Hex64.hpp
    cl_HMR_Lagrange_Element_Hex8.hpp
    cl_HMR_Lagrange_Element.hpp
    cl_HMR_Lagrange_Element_Quad16.hpp
    cl_HMR_Lagrange_Element_Quad4.hpp
    cl_HMR_Lagrange_Element_Quad9.hpp
    cl_HMR_Lagrange_Facet.hpp
    cl_HMR_Lagrange_Facet_Line2.hpp
    cl_HMR_Lagrange_Facet_Line3.hpp
    cl_HMR_Lagrange_Facet_Line4.hpp
    cl_HMR_Lagrange_Facet_Quad16.hpp
    cl_HMR_Lagrange_Facet_Quad4.hpp
    cl_HMR_Lagrange_Facet_Quad9.hpp
    cl_HMR_Lagrange_Mesh_Base.hpp
    cl_HMR_Lagrange_Mesh.hpp
    cl_HMR_Lagrange_Node.hpp
    cl_HMR_Lagrange_Node_Interpolation.hpp
    cl_HMR_Mesh_Base.hpp
    cl_HMR_Mesh.hpp
    cl_HMR_Parameters.hpp
    cl_HMR_State.hpp
    cl_HMR_STK.hpp
    cl_HMR_T_Matrix.hpp
    fn_HMR_Background_Element_Edges_3D.hpp
    fn_HMR_Background_Element_Neighbors_2D.hpp
    fn_HMR_Background_Element_Neighbors_3D.hpp
    fn_HMR_get_basis_neighbors_2d.hpp
    fn_HMR_get_basis_neighbors_3d.hpp
    HMR_Globals.hpp
    HMR_Tools.hpp
)

# List header dependencies
set(HMR_INCLUDES
    ${COM}
    ${LINALG}
    ${MTK}
    ${MTK}/${MAP}
#    ${ALG} #MTK
#    ${TOL} #MTK
    ${XTK} #MTK
#    ${STK} #MTK
    ${FEM}/${INT} # for L2
    ${FEM}/${MSI} # for L2
    ${FEM}/${MDL} # for L2
    )

# Include depedency directories
foreach(HMR_INCLUDE ${HMR_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${HMR_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Library

# List library source files
set(LIB_SOURCES
    cl_HMR_Arguments.cpp
    cl_HMR_Background_Edge.cpp
    cl_HMR_Background_Facet.cpp
    cl_HMR_Background_Mesh_Base.cpp
    cl_HMR_BSpline_Mesh_Base.cpp
    cl_HMR_Database.cpp
    cl_HMR.cpp
    cl_HMR_Edge.cpp
    cl_HMR_Element.cpp
    cl_HMR_Facet.cpp
    cl_HMR_Factory.cpp
    cl_HMR_Field.cpp
    cl_HMR_File.cpp
    cl_HMR_Lagrange_Mesh_Base.cpp
    cl_HMR_Lagrange_Node_Interpolation.cpp
    cl_HMR_Mesh_Base.cpp
    cl_HMR_Mesh.cpp
    cl_HMR_Parameters.cpp
    cl_HMR_STK.cpp
    cl_HMR_T_Matrix.cpp
   )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${HMR}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${HMR}-lib ${LIB_DEPENDENCIES})
set_target_properties(${HMR}-lib PROPERTIES OUTPUT_NAME ${HMR})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# project Executable

# List executable source files
set(EXE_SOURCES
    main.cpp )

# List executable dependencies
set(EXE_DEPENDENCIES
    ${HMR}-lib
    ${MAP}-lib
    ${MTK}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
#    ${ALG}-lib #MTK
#    ${TOL}-lib #MTK
    ${SDF}-lib # for SDF
#    ${XTK}-lib #MTK
#    ${STK}-lib #MTKz
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
add_executable(${HMR}-exe ${EXE_SOURCES})
target_link_libraries(${HMR}-exe ${EXE_DEPENDENCIES})
set_target_properties(${HMR}-exe PROPERTIES OUTPUT_NAME hmr )




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${HMR}Config.cmake.in
    ${${HMR}_CONFIG_BUILD_DIR}/${HMR}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${HMR}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${HMR}_CONFIG_BUILD_DIR}/${HMR}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Installation

# Install *Config and *ConfigVersion files
install(FILES ${${HMR}_CONFIG_BUILD_DIR}/${HMR}Config.cmake 
    ${${HMR}_CONFIG_BUILD_DIR}/${HMR}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${HMR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${HMR} )

# Install library and executable
install(TARGETS ${HMR}-lib
    EXPORT ${HMR}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${HMR}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${HMR}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${HMR} )

# Export library in to *Targets file
install(EXPORT ${HMR}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${HMR} )
