/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Field.cpp
 *
 */

#include "catch.hpp"

#include "paths.hpp"

// implementations to test
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src

#include "fn_PRM_HMR_Parameters.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

#define protected public
#define private public
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Analytic.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"
#undef protected
#undef private

namespace moris
{
    //---------------------------------------------------------------------------
    moris::real
    LevelSetFunctionForHMR( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        return norm( aPoint ) - 1.2;
    }

    moris::real
    LevelSetFunction(
            const moris::Matrix< moris::DDRMat >& aPoint,
            const moris::Matrix< moris::DDRMat >& aParameters )
    {
        return norm( aPoint ) - aParameters( 0 );
    }

    moris::real
    LevelSetPlaneFunction(
            const moris::Matrix< moris::DDRMat >& aPoint,
            const moris::Matrix< moris::DDRMat >& aParameters )
    {
        return aPoint( 0 ) - aParameters( 0 );
    }

    void
    DummyDerivativeFunction(
            const moris::Matrix< DDRMat >& aCoordinates,
            const moris::Matrix< DDRMat >& aParameters,
            moris::Matrix< DDRMat >&       aReturnValue )
    {
    }

    //---------------------------------------------------------------------------

    namespace mtk
    {
        TEST_CASE( "MTK Field", "[MTK],[MTK_Field]" )
        {
            if ( par_size() == 1 )
            {
                uint        tLagrangeMeshIndex     = 0;
                uint        tLagrangeMeshIndex_Out = 1;
                std::string tFieldName             = "Cylinder";

                Parameter_List tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4" ) );
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes", std::string( "0" ) );

                tParameters.set( "lagrange_orders", std::string( "1,1" ) );
                tParameters.set( "lagrange_pattern", std::string( "0,1" ) );
                tParameters.set( "bspline_orders", std::string( "1,1" ) );
                tParameters.set( "bspline_pattern", std::string( "0,1" ) );

                tParameters.set( "lagrange_to_bspline", "0;1" );

                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 0 );
                tParameters.set( "staircase_buffer", 1 );
                tParameters.set( "initial_refinement", "0,0" );
                tParameters.set( "initial_refinement_pattern", "0,1" );

                tParameters.set( "use_number_aura", 0 );

                tParameters.set( "use_multigrid", 0 );
                tParameters.set( "severity_level", 2 );

                hmr::HMR tHMR( tParameters );

                // initial refinement
                tHMR.perform_initial_refinement();

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                //// create field
                std::shared_ptr< moris::hmr::Field > tFieldHMR = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );

                for ( uint k = 0; k < 2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );
                }

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for ( uint k = 0; k < 2; ++k )
                {
                    for ( uint Ik = 0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
                    {
                        tHMR.get_database()->get_background_mesh()->get_element( Ik )->remove_from_refinement_queue();
                    }
                    tHMR.get_database()->get_background_mesh()->get_element( 5 )->put_on_refinement_queue();

                    // refine mesh
                    tHMR.get_database()->get_background_mesh()->perform_refinement( 1 );
                }

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 0 );

                tHMR.get_database()->update_bspline_meshes();
                tHMR.get_database()->update_lagrange_meshes();

                tHMR.finalize();

                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

                // Create integration mesh
                mtk::Integration_Mesh* tIntegrationMesh =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh );

                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_Out );

                // Create integration mesh
                mtk::Integration_Mesh* tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex_Out );

                // Create mesh manager
                mtk::Mesh_Pair tMeshPair_In( tInterpolationMesh, tIntegrationMesh, true );

                mtk::Mesh_Pair tMeshPair_Out( tInterpolationMesh_Out, tIntegrationMesh_Out, true );

                // Define two analtyic MTK fields
                Matrix< DDRMat > tParam = { { 0.2 } };

                mtk::Field_Analytic* tField_In = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_In );

                mtk::Field_Analytic* tField_Out = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_Out );

                // Get vector of nodal values and check for correct size
                // std::cout<<"Field_In size:  "<<tField_In->get_values().numel()<<std::endl;
                // std::cout<<"Field_Out size: "<<tField_Out->get_values().numel()<<std::endl;

                CHECK( equal_to( tField_In->get_values().numel(), 125 ) );
                CHECK( equal_to( tField_Out->get_values().numel(), 46 ) );

                //                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );
                //                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                delete tField_In;
                delete tField_Out;
            }
        }

        TEST_CASE( "MTK Map Field Linear", "[MTK],[MTK_Map_Field_Linear]" )
        {
            if ( par_size() == 1 )
            {
                uint        tLagrangeMeshIndex   = 0;
                uint        tLagrangeMeshIndex_2 = 1;
                std::string tFieldName           = "Cylinder";

                Parameter_List tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4" ) );
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes", std::string( "0" ) );

                tParameters.set( "lagrange_orders", std::string( "1,1" ) );
                tParameters.set( "lagrange_pattern", std::string( "0,1" ) );
                tParameters.set( "bspline_orders", std::string( "1,1" ) );
                tParameters.set( "bspline_pattern", std::string( "0,1" ) );

                tParameters.set( "lagrange_to_bspline", "0;1" );

                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 0 );
                tParameters.set( "staircase_buffer", 1 );
                tParameters.set( "initial_refinement", "0,0" );
                tParameters.set( "initial_refinement_pattern", "0,1" );

                tParameters.set( "use_number_aura", 0 );

                tParameters.set( "use_multigrid", 0 );
                tParameters.set( "severity_level", 2 );

                hmr::HMR tHMR( tParameters );

                // initial refinement
                tHMR.perform_initial_refinement();

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                //// create field
                std::shared_ptr< moris::hmr::Field > tFieldHMR = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );

                for ( uint k = 0; k < 2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );
                }

                //-------- next mesh ---------------

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for ( uint k = 0; k < 2; ++k )
                {
                    for ( uint Ik = 0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
                    {
                        tHMR.get_database()->get_background_mesh()->get_element( Ik )->remove_from_refinement_queue();
                    }
                    tHMR.get_database()->get_background_mesh()->get_element( 10 )->put_on_refinement_queue();

                    // refine mesh
                    tHMR.get_database()->get_background_mesh()->perform_refinement( 1 );
                }

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 0 );

                tHMR.get_database()->update_bspline_meshes();
                tHMR.get_database()->update_lagrange_meshes();

                tHMR.finalize();

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                mtk::Integration_Mesh*              tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex );

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_In = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_2 );
                mtk::Integration_Mesh*              tIntegrationMesh_In =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_In, tLagrangeMeshIndex_2 );

                // Create mesh manager
                mtk::Mesh_Pair tMeshPair_In( tInterpolationMesh_In, tIntegrationMesh_In, true );

                mtk::Mesh_Pair tMeshPair_Out( tInterpolationMesh_Out, tIntegrationMesh_Out, true );

                // Define analytic MTK field as input field
                Matrix< DDRMat > tParam = { { 0.2 } };

                mtk::Field_Analytic* tField_In = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_In );

                // Define discretized MTK field as output field
                mtk::Field_Discrete* tField_Out = new mtk::Field_Discrete(
                        tMeshPair_Out );

                // check that input field can return vector nodal values and has correct size
                CHECK( equal_to( tField_In->get_values().numel(), 46 ) );

                // Use mapper
                mtk::Mapper tMapper;
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tMapper.perform_mapping(
                        tField_Out,
                        EntityRank::BSPLINE,
                        EntityRank::NODE );

                tFieldHMR->get_node_values() = tField_Out->get_values();

                // tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                // tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                // create analytic MTK field on output mesh
                mtk::Field* tField_Ref = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_Out );

                CHECK( equal_to( tField_Out->get_values().numel(), 125 ) );
                CHECK( equal_to( tField_Ref->get_values().numel(), 125 ) );

                for ( uint Ik = 0; Ik < tField_Out->get_values().numel(); Ik++ )
                {
                    CHECK( equal_to( tField_Out->get_values()( Ik ), tField_Ref->get_values()( Ik ) ) );
                }

                delete tField_In;
                delete tField_Out;
                delete tField_Ref;
            }
        }

        TEST_CASE( "MTK Map Field Quad to Linear", "[MTK],[MTK_Map_Field_Quad_to_Linear]" )
        {
            if ( par_size() == 1 )
            {
                uint        tLagrangeMeshIndex   = 0;
                uint        tLagrangeMeshIndex_2 = 1;
                std::string tFieldName           = "Cylinder";

                Parameter_List tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4" ) );
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes", std::string( "0" ) );

                tParameters.set( "lagrange_orders", std::string( "1,2" ) );
                tParameters.set( "lagrange_pattern", std::string( "0,1" ) );
                tParameters.set( "bspline_orders", std::string( "1,2" ) );
                tParameters.set( "bspline_pattern", std::string( "0,1" ) );

                tParameters.set( "lagrange_to_bspline", "0;1" );

                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 0 );
                tParameters.set( "staircase_buffer", 1 );
                tParameters.set( "initial_refinement", "0,0" );
                tParameters.set( "initial_refinement_pattern", "0,1" );

                tParameters.set( "use_number_aura", 0 );

                tParameters.set( "use_multigrid", 0 );
                tParameters.set( "severity_level", 2 );

                hmr::HMR tHMR( tParameters );

                // initial refinement
                tHMR.perform_initial_refinement();

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                //// create field
                std::shared_ptr< moris::hmr::Field > tFieldHMR = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );

                for ( uint k = 0; k < 2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );
                }

                //-------- next mesh ---------------

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for ( uint k = 0; k < 2; ++k )
                {
                    for ( uint Ik = 0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
                    {
                        tHMR.get_database()->get_background_mesh()->get_element( Ik )->remove_from_refinement_queue();
                    }
                    tHMR.get_database()->get_background_mesh()->get_element( 10 )->put_on_refinement_queue();

                    // refine mesh
                    tHMR.get_database()->get_background_mesh()->perform_refinement( 1 );
                }

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 0 );

                tHMR.get_database()->update_bspline_meshes();
                tHMR.get_database()->update_lagrange_meshes();

                tHMR.finalize();

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                mtk::Integration_Mesh*              tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex );

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_In = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_2 );
                mtk::Integration_Mesh*              tIntegrationMesh_In =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_In, tLagrangeMeshIndex_2 );

                mtk::Mesh_Pair tMeshPair_In( tInterpolationMesh_In, tIntegrationMesh_In, true );

                mtk::Mesh_Pair tMeshPair_Out( tInterpolationMesh_Out, tIntegrationMesh_Out, true );

                // Define analytic MTK field as input field
                Matrix< DDRMat > tParam = { { 0.2 } };

                mtk::Field_Analytic* tField_In = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_In );

                // Define discretized MTK field as output field
                mtk::Field_Discrete* tField_Out = new mtk::Field_Discrete(
                        tMeshPair_Out );

                CHECK( equal_to( tField_In->get_values().numel(), 217 ) );

                // Use mapper
                mtk::Mapper tMapper;
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tMapper.perform_mapping(
                        tField_Out,
                        EntityRank::BSPLINE,
                        EntityRank::NODE );

                // create analytic MTK field on output mesh
                mtk::Field* tField_Ref = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_Out );

                // create analytic MTK field on output mesh for loading exodus data
                mtk::Field* tField_Exodus = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_Out );

                tField_In->save_field_to_exodus( "./mtk_field_test_in.e" );
                tField_Out->save_field_to_exodus( "./mtk_field_test_out.e" );
                tField_Ref->save_field_to_exodus( "./mtk_field_test_ref.e" );

                tField_Exodus->load_field_from_exodus( "./mtk_field_test_ref.e" );

                CHECK( equal_to( tField_Out->get_values().numel(), 133 ) );
                CHECK( equal_to( tField_Ref->get_values().numel(), 133 ) );

                // compare mapped field against reference field (direct and loaded from exodus)
                for ( uint Ik = 0; Ik < tField_Out->get_values().numel(); Ik++ )
                {
                    CHECK( equal_to( tField_Out->get_values()( Ik ), tField_Ref->get_values()( Ik ) ) );
                    CHECK( equal_to( tField_Out->get_values()( Ik ), tField_Exodus->get_values()( Ik ) ) );
                }

                delete tField_In;
                delete tField_Out;
                delete tField_Ref;
                delete tField_Exodus;
            }
        }

        TEST_CASE( "MTK Map Field Linear to Quad", "[MTK],[MTK_Map_Field_Linear_to_Quad]" )
        {
            if ( par_size() == 1 )
            {
                uint        tLagrangeMeshIndex   = 0;
                uint        tLagrangeMeshIndex_2 = 1;
                std::string tFieldName           = "Cylinder";

                Parameter_List tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4" ) );
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes", std::string( "0" ) );

                tParameters.set( "lagrange_orders", std::string( "2,1" ) );
                tParameters.set( "lagrange_pattern", std::string( "0,1" ) );
                tParameters.set( "bspline_orders", std::string( "2,1" ) );
                tParameters.set( "bspline_pattern", std::string( "0,1" ) );

                tParameters.set( "lagrange_to_bspline", "0;1" );

                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 0 );
                tParameters.set( "staircase_buffer", 1 );
                tParameters.set( "initial_refinement", "0,0" );
                tParameters.set( "initial_refinement_pattern", "0,1" );

                tParameters.set( "use_number_aura", 0 );

                tParameters.set( "use_multigrid", 0 );
                tParameters.set( "severity_level", 2 );

                hmr::HMR tHMR( tParameters );

                // initial refinement
                tHMR.perform_initial_refinement();

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                //// create field
                std::shared_ptr< moris::hmr::Field > tFieldHMR = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );

                for ( uint k = 0; k < 2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunctionForHMR );
                }

                //-------- next mesh ---------------

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for ( uint k = 0; k < 2; ++k )
                {
                    for ( uint Ik = 0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
                    {
                        tHMR.get_database()->get_background_mesh()->get_element( Ik )->remove_from_refinement_queue();
                    }
                    tHMR.get_database()->get_background_mesh()->get_element( 10 )->put_on_refinement_queue();

                    // refine mesh
                    tHMR.get_database()->get_background_mesh()->perform_refinement( 1 );
                }

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 0 );

                tHMR.get_database()->update_bspline_meshes();
                tHMR.get_database()->update_lagrange_meshes();

                tHMR.finalize();

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                mtk::Integration_Mesh*              tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex );

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh_In = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_2 );
                mtk::Integration_Mesh*              tIntegrationMesh_In =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_In, tLagrangeMeshIndex_2 );

                mtk::Mesh_Pair tMeshPair_In( tInterpolationMesh_In, tIntegrationMesh_In, true );

                mtk::Mesh_Pair tMeshPair_Out( tInterpolationMesh_Out, tIntegrationMesh_Out, true );

                // Define analytic MTK field as input field
                Matrix< DDRMat > tParam = { { 0.2 } };

                mtk::Field_Analytic* tField_In = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_In );

                // Define discretized MTK field as output field
                mtk::Field_Discrete* tField_Out = new mtk::Field_Discrete(
                        tMeshPair_Out );

                CHECK( equal_to( tField_In->get_values().numel(), 63 ) );

                // Use mapper
                mtk::Mapper tMapper;
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tMapper.perform_mapping(
                        tField_Out,
                        EntityRank::BSPLINE,
                        EntityRank::NODE );

                tFieldHMR->get_node_values() = tField_Out->get_values();

                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                // create analytic MTK field on output mesh
                mtk::Field* tField_Ref = new mtk::Field_Analytic(
                        LevelSetPlaneFunction,
                        DummyDerivativeFunction,
                        tParam,
                        tMeshPair_Out );

                CHECK( equal_to( tField_Out->get_values().numel(), 465 ) );

                for ( uint Ik = 0; Ik < tField_Out->get_values().numel(); Ik++ )
                {
                    CHECK( equal_to( tField_Out->get_values()( Ik ), tField_Ref->get_values()( Ik ) ) );
                }

                delete tField_In;
                delete tField_Out;
                delete tField_Ref;
            }
        }

    }    // namespace mtk
}    // namespace moris
