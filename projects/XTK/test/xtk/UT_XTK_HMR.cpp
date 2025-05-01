/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_HMR.cpp
 *
 */

#include "catch.hpp"
#include <memory>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Vertex_Enrichment.hpp"
#include "cl_XTK_Diagnostics.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR.hpp"

#include "fn_PRM_HMR_Parameters.hpp"
#include "cl_GEN_User_Defined_Field.hpp"

namespace moris
{

    moris::real
    LevelSetSphereCylinder( const moris::Matrix< DDRMat >& aPoint )
    {
        moris::Matrix< DDRMat > aCenter = { { 0.01 }, { 0.01 }, { 0.01 } };
        moris::Matrix< DDRMat > aAxis   = { { 0.0 }, { 1.0 }, { 0.0 } };
        moris::real             aRad    = 0.77;
        moris::real             aLength = 5;

        MORIS_ASSERT( aCenter.numel() == 3, "Centers need to have length 3" );
        MORIS_ASSERT( aAxis.numel() == 3, "axis need to have length 3" );

        Vector< moris::real > relativePosition = { ( aPoint( 0 ) - aCenter( 0 ) ), ( aPoint( 1 ) - aCenter( 1 ) ), ( aPoint( 2 ) - aCenter( 2 ) ) };
        moris::real           lsFromLeft       = ( relativePosition( 0 ) * ( -aAxis( 0 ) ) + relativePosition( 1 ) * ( -aAxis( 1 ) ) + relativePosition( 2 ) * ( -aAxis( 2 ) ) ) - aLength / 2.0;
        moris::real           lsFromRight      = ( relativePosition( 0 ) * ( aAxis( 0 ) ) + relativePosition( 1 ) * ( aAxis( 1 ) ) + relativePosition( 2 ) * ( aAxis( 2 ) ) ) - aLength / 2.0;

        moris::real           axialCrd  = ( relativePosition( 0 ) * ( aAxis( 0 ) ) + relativePosition( 1 ) * ( aAxis( 1 ) ) + relativePosition( 2 ) * ( aAxis( 2 ) ) );
        Vector< moris::real > radDir    = { ( relativePosition( 0 ) - aAxis( 0 ) * axialCrd ), ( relativePosition( 1 ) - aAxis( 1 ) * axialCrd ), ( relativePosition( 2 ) - aAxis( 2 ) * axialCrd ) };
        moris::real           radDist   = std::pow( radDir( 0 ) * radDir( 0 ) + radDir( 1 ) * radDir( 1 ) + radDir( 2 ) * radDir( 2 ), 0.5 );
        moris::real           lsFromRad = radDist - aRad;

        return -std::max( std::max( lsFromLeft, lsFromRight ), lsFromRad );
    }

    real LevelSetSphereCylinderGeometry( const Matrix< DDRMat >& aCoordinates, const Vector< real >& aParameters )
    {
        return LevelSetSphereCylinder( aCoordinates );
    }

    moris::real
    LevelSetPlaneFunction( const moris::Matrix< DDRMat >& aPoint )
    {

        real mXn = 0.0;
        real mYn = 0.0;
        real mZn = 1.0;
        real mXc = 0.0;
        real mYc = 0.0;
        real mZc = 1.8;
        return mXn * ( aPoint( 0 ) - mXc ) + mYn * ( aPoint( 1 ) - mYc ) + mZn * ( aPoint( 2 ) - mZc );
    }

    TEST_CASE( "XTK HMR Test", "[XTK_HMR]" )
    {

        // iterate through linear quadratic cubic
        for ( moris::uint iOrder = 1; iOrder < 2; iOrder++ )
        {

            Module_Parameter_Lists tParameterlist( Module_Type::HMR );
            tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
            tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", 6, 6, 6 );
            tParameterlist( 0 )( 0 ).set( "domain_dimensions", 2.0, 2.0, 2.0 );
            tParameterlist( 0 )( 0 ).set( "domain_offset", -1.0, 1.0, -2.0 );
            tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
            tParameterlist( 0 )( 0 ).set( "lagrange_orders", std::to_string( iOrder ) );
            tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
            tParameterlist( 0 )( 0 ).set( "bspline_orders", std::to_string( iOrder ) );
            tParameterlist( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );
            tParameterlist( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
            tParameterlist( 0 )( 0 ).set( "refinement_buffer", 1 );
            tParameterlist( 0 )( 0 ).set( "staircase_buffer", 1 );

            std::shared_ptr< hmr::HMR > tHMR = std::make_shared< hmr::HMR >( tParameterlist );

            // initialize a mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();

            tMeshManager->set_performer( tHMR );

            tHMR->set_performer( tMeshManager );

            tHMR->perform_initial_refinement();
            tHMR->perform();

            auto                                       tField          = std::make_shared< moris::gen::User_Defined_Field >( &( LevelSetSphereCylinderGeometry ), Vector< gen::ADV >() );
            Vector< std::shared_ptr< gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tField ) };

            size_t                                 tModelDimension = 3;
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tMeshManager->get_interpolation_mesh( 0 ), tGeometryEngineParameters );
            xtk::Model                  tXTKModel( tModelDimension, tMeshManager->get_interpolation_mesh( 0 ), &tGeometryEngine );
            tXTKModel.mVerbose = true;

            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };

            // Do the cutting
            tXTKModel.decompose( tDecompositionMethods );

            CHECK( xtk::interpolated_coordinate_check( tXTKModel.get_cut_integration_mesh() ) );
        }
    }
}    // namespace moris
