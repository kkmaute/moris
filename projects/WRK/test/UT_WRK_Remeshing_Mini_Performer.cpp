/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Remeshing_Mini_Performer.cpp
 *
 */

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

#include "cl_Library_IO.hpp"

#include "moris_typedefs.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Field_Analytic.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

#include "cl_Matrix.hpp"    //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"    // ALG/src

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Factory.hpp"               //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Database.hpp"              //HMR/src

#include "cl_WRK_perform_refinement.hpp"
#include "cl_WRK_perform_remeshing.hpp"

#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"

using namespace moris;

moris::real
tCircleForHMR( const Matrix< DDRMat > &aCoords )
{
    return norm( aCoords ) - 1.5237;
}

moris::real
tCircle(
        const Matrix< DDRMat > &aCoords,
        const Matrix< DDRMat > &aParameters )
{
    return norm( aCoords ) - aParameters( 0 );
}

void DummyDerivativeFunction(
        const moris::Matrix< DDRMat > &aCoordinates,
        const moris::Matrix< DDRMat > &aParameters,
        moris::Matrix< DDRMat >       &aReturnValue )
{
}

TEST_CASE( "WRK L2 test", "[WRK_L2_test]" )
{
    if ( par_size() <= 1 )
    {
        std::string tFieldName = "Geometry";

        //----- HRM parameter list --------

        Parameter_List tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", 4, 4 );
        tParameters.set( "domain_dimensions", 4.0, 4.0 );
        tParameters.set( "domain_offset", -2.0, -2.0 );
        tParameters.set( "lagrange_output_meshes", "0" );

        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_to_bspline", "0" );

        tParameters.set( "refinement_buffer", 0 );
        tParameters.set( "staircase_buffer", 1 );
        tParameters.set( "initial_refinement", "1" );

        Parameter_List tRefinementParameters( "Refinement" );
        prm::create_refinement_parameterlist( tRefinementParameters );
        tRefinementParameters.set( "field_names", "Circle" );
        tRefinementParameters.set( "levels_of_refinement", "1" );
        tRefinementParameters.set( "refinement_pattern", "0" );

        Parameter_List tRemeshingParameters( "Remeshing" );
        prm::create_remeshing_parameterlist( tRemeshingParameters );
        tRemeshingParameters.set( "mode", "ab_initio" );
        tRemeshingParameters.set( "remeshing_field_names", "Circle" );
        tRemeshingParameters.set( "remeshing_levels_of_refinement", "1" );
        tRemeshingParameters.set( "remeshing_refinement_pattern", "0" );

        //---------------------------------------------------------------------------------------
        //                               Stage 1: HMR refinement
        //---------------------------------------------------------------------------------------

        std::shared_ptr< hmr::HMR > tHMRPerformer = std::make_shared< hmr::HMR >( tParameters );

        // uniform initial refinement
        tHMRPerformer->perform_initial_refinement();

        std::shared_ptr< hmr::Database > tHMRDatabase = tHMRPerformer->get_database();

        hmr::Interpolation_Mesh_HMR *tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                tHMRDatabase,
                1,
                0,
                1,
                0 );

        mtk::Mesh_Pair tMeshPair( tInterpolationMesh, nullptr, false );

        // Define two analytic MTK fields
        Vector< std::shared_ptr< mtk::Field > > tFields( 1, nullptr );
        Vector< std::shared_ptr< mtk::Field > > tFieldsOut;

        moris::Matrix< DDRMat > tMat = { { 1.5237 } };

        tFields( 0 ) = std::make_shared< mtk::Field_Analytic >( tCircle, DummyDerivativeFunction, tMat, tMeshPair, 1 );

        tFields( 0 )->set_label( "Circle" );

        tFields( 0 )->save_field_to_exodus( "Remeshing_Field1.exo" );

        wrk::Refinement_Mini_Performer tRefinementPerformer( tRefinementParameters );
        tRefinementPerformer.perform_refinement( tFields, tHMRPerformer );

        // delete( tInterpolationMesh );

        std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer_HMR = std::make_shared< mtk::Mesh_Manager >();
        tHMRPerformer->set_performer( tMTKPerformer_HMR );

        hmr::Interpolation_Mesh_HMR *tInterpolationMeshNew = new hmr::Interpolation_Mesh_HMR(
                tHMRDatabase,
                1,
                0,
                1,
                0 );

        mtk::Mesh_Pair tMeshPairNew( tInterpolationMeshNew, nullptr, false );

        tFields( 0 )->unlock_field();
        tFields( 0 )->set_mesh_pair( tMeshPairNew );
        tFields( 0 )->save_field_to_exodus( "Remeshing_Field2.exo" );

        tFields( 0 )->unlock_field();
        tFields( 0 )->set_coefficients( { { 0.8713 } } );
        tFields( 0 )->save_field_to_exodus( "Remeshing_Field3.exo" );

        // print( tFields( 0 )->get_values(), "Val1");

        tHMRDatabase->get_background_mesh()->update_database();
        tHMRDatabase->update_bspline_meshes();
        tHMRDatabase->update_lagrange_meshes();

        // HMR finalize
        tHMRPerformer->perform();

        Vector< std::shared_ptr< hmr::HMR > >          tHMRPerformers( 1, tHMRPerformer );
        Vector< std::shared_ptr< mtk::Mesh_Manager > > tMTKPerformers( 1 );

        wrk::Remeshing_Mini_Performer tRemeshingPerformer( tRemeshingParameters );
        tRemeshingPerformer.perform_remeshing( tFields, tHMRPerformers, tMTKPerformers, tFieldsOut );

        std::shared_ptr< hmr::Database > tHMRDatabaseNew = tHMRPerformers( 0 )->get_database();

        hmr::Interpolation_Mesh_HMR *tInterpolationMeshNewMesh = new hmr::Interpolation_Mesh_HMR(
                tHMRDatabaseNew,
                1,
                0,
                1,
                0 );

        mtk::Mesh_Pair tMeshPairNewMesh( tInterpolationMeshNewMesh, nullptr, false );
        tFields( 0 )->unlock_field();
        tFields( 0 )->set_mesh_pair( tMeshPairNewMesh );
        tFields( 0 )->save_field_to_exodus( "Remeshing_Field4.exo" );

        // explicitly delete meshes
        delete tInterpolationMesh;
        tInterpolationMesh = nullptr;

        delete tInterpolationMeshNewMesh;
        tInterpolationMeshNewMesh = nullptr;

        delete tInterpolationMeshNew;
        tInterpolationMeshNew = nullptr;
    }
}
