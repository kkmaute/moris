/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_Snapping.cpp
 *
 */

#include "catch.hpp"
#include <string>
#include "cl_Vector.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Geometry_Engine.hpp"

#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_Module_Parameter_Lists.hpp"

namespace moris::xtk
{

    TEST_CASE( "Plane Aligned with Background 2d", "[Snapping_2d]" )
    {
        if ( par_size() == 1 )
        {
            // Intersecting plane test parameters
            real tPlaneAngle                  = 90;      // CCW, degrees
            real tIntersectionLocalCoordinate = -0.5;    // (-1, 1), where the top middle vertical edge is intersected
            // bool tSnapFromRight = true;
            // bool tSnapBoth = false;
            bool tBilinear = false;

            // XTK parameter list
            Module_Parameter_Lists tXTKParams( Module_Type::XTK );
            tXTKParams( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
            tXTKParams( 0 )( 0 ).set( "decompose", true );
            tXTKParams( 0 )( 0 ).set( "decomposition_type", "conformal" );
            tXTKParams( 0 )( 0 ).set( "enrich", true );
            tXTKParams( 0 )( 0 ).set( "basis_rank", "bspline" );
            tXTKParams( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
            tXTKParams( 0 )( 0 ).set( "ghost_stab", true );
            tXTKParams( 0 )( 0 ).set( "multigrid", false );
            tXTKParams( 0 )( 0 ).set( "verbose", true );
            tXTKParams( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
            tXTKParams( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
            tXTKParams( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
            tXTKParams( 0 )( 0 ).set( "probe_bg_cells", "1,2" );

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // GEN Parameters
            Module_Parameter_Lists tGENParams( Module_Type::GEN );

            // Calculations
            real tOffset  = 0.5 * ( 1 + tIntersectionLocalCoordinate );
            real tXNormal = sin( -tPlaneAngle * M_PI / 180 );
            real tYNormal = cos( -tPlaneAngle * M_PI / 180 );

            tGENParams( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
            tGENParams( 0 )( 0 ).set( "output_mesh_file", "GEN_Snapping.exo" );

            // Geometry parameter lists
            moris::uint tGeoCounter        = 0;
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_threshold", 1e-16 );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "intersection_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "center_x", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "center_y", tOffset );
            tGENParams( 1 )( tGeoCounter ).set( "normal_x", tXNormal );
            tGENParams( 1 )( tGeoCounter ).set( "normal_y", tYNormal );
            tGENParams( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tBilinear );
            tGeoCounter++;

            // Geometry parameter lists
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_threshold", 1e-16 );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "intersection_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "center_x", -0.5 );
            tGENParams( 1 )( tGeoCounter ).set( "center_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tBilinear );
            tGeoCounter++;

            // Geometry parameter lists
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_threshold", 1e-16 );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "intersection_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "center_x", 1.0 );
            tGENParams( 1 )( tGeoCounter ).set( "center_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tBilinear );
            tGeoCounter++;

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // HMR parameters
            std::string tInterpolationOrder = "1";

            int                                      tRefineBuffer = 1;
            Module_Parameter_Lists tHMRParams( Module_Type::HMR );
            tHMRParams( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
            tHMRParams( 0 )( 0 ).set( "number_of_elements_per_dimension", 2, 1 );
            tHMRParams( 0 )( 0 ).set( "domain_dimensions", 2.0, 2.0 );
            tHMRParams( 0 )( 0 ).set( "domain_offset", -1.0, -1.0 );
            tHMRParams( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
            tHMRParams( 0 )( 0 ).set( "lagrange_orders", tInterpolationOrder );
            tHMRParams( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
            tHMRParams( 0 )( 0 ).set( "bspline_orders", tInterpolationOrder );
            tHMRParams( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );
            tHMRParams( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
            tHMRParams( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
            tHMRParams( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
            tHMRParams( 0 )( 0 ).set( "initial_refinement", "0" );
            tHMRParams( 0 )( 0 ).set( "initial_refinement_pattern", "0" );
            tHMRParams( 0 )( 0 ).set( "use_number_aura", 1 );
            tHMRParams( 0 )( 0 ).set( "use_multigrid", 0 );
            tHMRParams( 0 )( 0 ).set( "severity_level", 0 );
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // PSEUDO Workflow
            // HMR initialize
            std::shared_ptr< hmr::HMR > pHMR = std::make_shared< hmr::HMR >( tHMRParams( 0 )( 0 ) );

            // Geometry engine initialize
            std::shared_ptr< gen::Geometry_Engine > pGEN = std::make_shared< gen::Geometry_Engine >( tGENParams, nullptr );

            // Initialize  Background Mesh Mesh
            std::shared_ptr< mtk::Mesh_Manager > pBGMTK = std::make_shared< mtk::Mesh_Manager >();

            // Initialize output mesh from XTK
            std::shared_ptr< mtk::Mesh_Manager > pOutputMTK = std::make_shared< mtk::Mesh_Manager >();

            // XTK initialize
            std::shared_ptr< xtk::Model > pXTK = std::make_shared< xtk::Model >( tXTKParams( 0 )( 0 ) );

            // Set performer to HMR
            pHMR->set_performer( pBGMTK );

            // Set XTK Cooperations
            pXTK->set_geometry_engine( pGEN.get() );
            pXTK->set_input_performer( pBGMTK );
            pXTK->set_output_performer( pOutputMTK );

            // uniform initial refinement
            pHMR->perform_initial_refinement();

            // HMR finalize
            pHMR->perform();
            //
            pGEN->distribute_advs( pBGMTK->get_mesh_pair( 0 ) );

            // Output GEN fields, if requested
            pGEN->output_fields( pBGMTK->get_interpolation_mesh( 0 ) );

            // XTK perform - decompose - enrich - ghost - multigrid
            pXTK->perform_decomposition();
            pXTK->perform_enrichment();
        }
    }

}    // namespace moris::xtk
