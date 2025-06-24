/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Surface_Mesh.cpp
 *
 */

#include "catch.hpp"
#include <cmath>
#include "paths.hpp"
#include <Kokkos_Core.hpp>

#include "parameters.hpp"
#include "cl_GEN_Surface_Mesh_Geometry.hpp"

#include "fn_check_equal.hpp"

#if MORIS_HAVE_ARBORX
namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Surface Mesh Regularization - Isotropic Laplacian", "[gen], [pdv], [regularization], [isotropic laplacian]," )
    {
        if ( par_size() == 1 )
        {
            // get root from environment
            std::string tMorisRoot = moris::get_base_moris_dir();

            // Create parameter list for surface mesh geometry
            Submodule_Parameter_Lists tFieldParameterLists( "GEOMETRIES" );
            tFieldParameterLists.add_parameter_list( prm::create_surface_mesh_geometry_parameter_list() );
            tFieldParameterLists.set( "file_path", tMorisRoot + "projects/GEN/test/data/triangle_sensitivity_oblique.obj" );
            tFieldParameterLists.set( "intersection_tolerance", 1e-9 );
            tFieldParameterLists.set( "regularization_type", gen::Regularization_Type::ISOTROPIC_LAPLACIAN );
            tFieldParameterLists.set( "regularization_factors", 0.1 );

            // Create surface mesh parameters struct from parameter list
            Surface_Mesh_Parameters tParameters( tFieldParameterLists( 0 ) );

            // Dummy mesh
            mtk::Mesh* tMesh = nullptr;

            // Other parameters needed to create surface mesh, all dummy
            Vector< ADV >                 tADVs;
            ADV_Manager                   tADVManager;
            Node_Manager                  tNodeManager( tMesh );
            std::shared_ptr< Library_IO > tLibrary = nullptr;

            // Create surface mesh geometry
            Surface_Mesh_Geometry tGeom( tParameters, tNodeManager, tADVs, tADVManager, tLibrary );

            Matrix< DDRMat > tVertexCoordinatesExpected = { { 2.1, 0.3875, 0.9750, 1.2875 }, { 1.1750, 0.5, -.1750, 0.5 } };

            // Regularize the surface mesh
            tGeom.regularize();

            // Get the regularized vertex coordinates
            Matrix< DDRMat > tNewCoords = tGeom.get_all_vertex_coordinates();

            CHECK_EQUAL( tNewCoords, tVertexCoordinatesExpected, );
        }
    }
    // TEST_CASE( "Surface Mesh Regularization - Isotropic Laplacian Sensitivity", "[gen], [pdv], [regularization], [isotropic laplacian sensitivity]," )
    // {
    //     if ( par_size() == 1 )
    //     {
    //         // get root from environment
    //         std::string tMorisRoot = moris::get_base_moris_dir();

    //         // Create parameter list for surface mesh geometry
    //         Submodule_Parameter_Lists tFieldParameterLists( "GEOMETRIES" );
    //         tFieldParameterLists.add_parameter_list( prm::create_surface_mesh_geometry_parameter_list() );
    //         tFieldParameterLists.set( "file_path", tMorisRoot + "projects/GEN/test/data/triangle_sensitivity_oblique.obj" );
    //         tFieldParameterLists.set( "intersection_tolerance", 1e-9 );
    //         tFieldParameterLists.set( "regularization_type", gen::Regularization_Type::ISOTROPIC_LAPLACIAN );

    //         // Create surface mesh parameters struct from parameter list
    //         Surface_Mesh_Parameters tParameters( tFieldParameterLists( 0 ) );

    //         // Dummy mesh
    //         mtk::Mesh* tMesh = nullptr;

    //         // Other parameters needed to create surface mesh, all dummy
    //         Vector< ADV >                 tADVs;
    //         ADV_Manager                   tADVManager;
    //         Node_Manager                  tNodeManager( tMesh );
    //         std::shared_ptr< Library_IO > tLibrary = nullptr;

    //         // Create surface mesh geometry
    //         Surface_Mesh_Geometry tGeom( tParameters, tNodeManager, tADVs, tADVManager, tLibrary );

    //         Matrix< DDRMat > tVertexCoordinatesExpected = { { 0.75, 1.625, 0.75, 1.625 }, { 0.5, 0.5, 0.5, 0.5 } };

    //         // Regularize the surface mesh
    //         tGeom.regularize();

    //         // Get the regularized vertex coordinates
    //         Matrix< DDRMat > tNewCoords = tGeom.get_all_vertex_coordinates();

    //         CHECK_EQUAL( tNewCoords, tVertexCoordinatesExpected, );
    //     }
    // }
    
}    // namespace moris::gen

#endif
