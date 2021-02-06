/*
 * UT_MTK_Field.cpp
 *
 *  Created on: Jan 21, 2021
 *      Author: schmidt
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
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "fn_PRM_HMR_Parameters.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

#include "cl_MTK_BSpline_Field.hpp"
#include "cl_MTK_Field_Proxy.hpp"
#include "cl_MTK_Mapper.hpp"

namespace moris
{
    moris::real
    LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
    {
        return norm( aPoint ) - 1.2;
    }

    moris::real
    LevelSetFunction1( const moris::Matrix< moris::DDRMat > & aPoint )
    {
        return norm( aPoint ) - 0.2;
    }

    moris::real
    LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
    {
        return  aPoint( 0 )  - 0.2;
    }

    namespace mtk
    {
        TEST_CASE("MTK Field","[MTK],[MTK_Field]")
                {
            if(par_size() ==1)
            {
                uint tLagrangeMeshIndex = 0;
                uint tLagrangeMeshIndex_Out = 1;
                std::string tFieldName = "Cylinder";

                ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4"));
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes",std::string( "0") );

                tParameters.set( "lagrange_orders", std::string("1,1" ));
                tParameters.set( "lagrange_pattern", std::string("0,1" ));
                tParameters.set( "bspline_orders", std::string("1,1" ));
                tParameters.set( "bspline_pattern", std::string("0,1" ));

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

                tFieldHMR->evaluate_scalar_function( LevelSetFunction );

                for( uint k=0; k<2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunction );
                }

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for( uint k=0; k<2; ++k )
                {
                    for( uint Ik=0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
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

                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

                // Create integration mesh
                mtk::Integration_Mesh* tIntegrationMesh =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh );

                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_Out );

                // Create integration mesh
                mtk::Integration_Mesh* tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex_Out );

                // Create mesh manager
                std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();

                // Register mesh pair
                uint tMeshIndex_In = tMeshManager->register_mesh_pair( tInterpolationMesh, tIntegrationMesh );
                uint tMeshIndex_Out = tMeshManager->register_mesh_pair( tInterpolationMesh_Out, tIntegrationMesh_Out );

                mtk::Field_Proxy tField_In( tMeshManager, tMeshIndex_In );
                mtk::Field_Proxy tField_Out( tMeshManager, tMeshIndex_Out );

                tField_In.evaluate_scalar_function( LevelSetPlaneFunction );
                tField_Out.evaluate_scalar_function( LevelSetPlaneFunction );

                std::cout<<"Field_In size: "<<tField_In.get_nodal_values(tInterpolationMesh).numel()<<std::endl;
                std::cout<<"Field_Out size: "<<tField_Out.get_nodal_values(tInterpolationMesh_Out).numel()<<std::endl;

                CHECK(equal_to( tField_In.get_nodal_values(tInterpolationMesh).numel(), 125));
                CHECK(equal_to( tField_Out.get_nodal_values(tInterpolationMesh_Out).numel(), 46));
            }
                }

        TEST_CASE("MTK Map Field Linear","[MTK],[MTK_Map_Field_Linear]")
        {
            if(par_size() ==1)
            {
                uint tLagrangeMeshIndex = 0;
                uint tLagrangeMeshIndex_2 = 1;
                std::string tFieldName = "Cylinder";

                ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4"));
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes",std::string( "0") );

                tParameters.set( "lagrange_orders", std::string("1,1" ));
                tParameters.set( "lagrange_pattern", std::string("0,1" ));
                tParameters.set( "bspline_orders", std::string("1,1" ));
                tParameters.set( "bspline_pattern", std::string("0,1" ));

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

                tFieldHMR->evaluate_scalar_function( LevelSetFunction );

                for( uint k=0; k<2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunction );
                }

                //-------- next mesh ---------------

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for( uint k=0; k<2; ++k )
                {
                    for( uint Ik=0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
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
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                mtk::Integration_Mesh* tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex );

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_In = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_2 );
                mtk::Integration_Mesh* tIntegrationMesh_In =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_In, tLagrangeMeshIndex_2 );

                // Create mesh manager
                std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();
                uint tMeshIndex_Out = tMeshManager->register_mesh_pair( tInterpolationMesh_Out, tIntegrationMesh_Out );
                uint tMeshIndex_In  = tMeshManager->register_mesh_pair( tInterpolationMesh_In, tIntegrationMesh_In );

                mtk::Field_Proxy tField_In( tMeshManager, tMeshIndex_In, 0 );
                mtk::BSpline_Field tField_Out( tMeshManager, tMeshIndex_Out, 0 );

                tField_In.evaluate_scalar_function( LevelSetPlaneFunction );

                CHECK(equal_to( tField_In.get_nodal_values(tInterpolationMesh_In).numel(), 46));;

                // Use mapper
                mtk::Mapper tMapper(tMeshManager, tMeshIndex_In);
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tMapper.perform_mapping(
                        tField_Out,
                        EntityRank::BSPLINE,
                        EntityRank::NODE);

                //tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                //tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                mtk::Field_Proxy tField_Ref( tMeshManager, tMeshIndex_Out, 0 );
                tField_Ref.evaluate_scalar_function( LevelSetPlaneFunction );

                CHECK(equal_to( tField_Out.get_nodal_values().numel(), 125));

                for( uint Ik = 0; Ik < tField_Out.get_nodal_values().numel(); Ik++ )
                {
                    CHECK(equal_to( tField_Out.get_nodal_values()( Ik ), tField_Ref.get_nodal_values(tInterpolationMesh_Out)( Ik )));
                }
            }
        }

        TEST_CASE("MTK Map Field Quad to Linear","[MTK],[MTK_Map_Field_Quad_to_Linear]")
        {
            if(par_size() ==1)
            {
                uint tLagrangeMeshIndex = 0;
                uint tLagrangeMeshIndex_2 = 1;
                std::string tFieldName = "Cylinder";

                ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4"));
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes",std::string( "0") );

                tParameters.set( "lagrange_orders", std::string("1,2" ));
                tParameters.set( "lagrange_pattern", std::string("0,1" ));
                tParameters.set( "bspline_orders", std::string("1,2" ));
                tParameters.set( "bspline_pattern", std::string("0,1" ));

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

                tFieldHMR->evaluate_scalar_function( LevelSetFunction );

                for( uint k=0; k<2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunction );
                }

                //-------- next mesh ---------------

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for( uint k=0; k<2; ++k )
                {
                    for( uint Ik=0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
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
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                mtk::Integration_Mesh* tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex );

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_In = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_2 );
                mtk::Integration_Mesh* tIntegrationMesh_In =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_In, tLagrangeMeshIndex_2 );

                // Create mesh manager
                std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();
                uint tMeshIndex_Out = tMeshManager->register_mesh_pair( tInterpolationMesh_Out, tIntegrationMesh_Out );
                uint tMeshIndex_In  = tMeshManager->register_mesh_pair( tInterpolationMesh_In, tIntegrationMesh_In );

                mtk::Field_Proxy tField_In( tMeshManager, tMeshIndex_In, 0 );
                mtk::BSpline_Field tField_Out( tMeshManager, tMeshIndex_Out, 0 );

                tField_In.evaluate_scalar_function( LevelSetPlaneFunction );

                CHECK(equal_to( tField_In.get_nodal_values(tInterpolationMesh_In).numel(), 217));;

                // Use mapper
                mtk::Mapper tMapper(tMeshManager, tMeshIndex_In);
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tMapper.perform_mapping(
                        tField_Out,
                        EntityRank::BSPLINE,
                        EntityRank::NODE);

                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                mtk::Field_Proxy tField_Ref( tMeshManager, tMeshIndex_Out, 0 );
                tField_Ref.evaluate_scalar_function( LevelSetPlaneFunction );

                CHECK(equal_to( tField_Out.get_nodal_values().numel(), 133));

                for( uint Ik = 0; Ik < tField_Out.get_nodal_values().numel(); Ik++ )
                {
                    CHECK(equal_to( tField_Out.get_nodal_values()( Ik ), tField_Ref.get_nodal_values(tInterpolationMesh_Out)( Ik )));
                }
            }
        }

        TEST_CASE("MTK Map Field Linear to Quad","[MTK],[MTK_Map_Field_Linear_to_Quad]")
        {
            if(par_size() ==1)
            {
                uint tLagrangeMeshIndex = 0;
                uint tLagrangeMeshIndex_2 = 1;
                std::string tFieldName = "Cylinder";

                ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::string( "4, 4"));
                tParameters.set( "domain_dimensions", "2, 2" );
                tParameters.set( "domain_offset", "-1.0, -1.0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes",std::string( "0") );

                tParameters.set( "lagrange_orders", std::string("2,1" ));
                tParameters.set( "lagrange_pattern", std::string("0,1" ));
                tParameters.set( "bspline_orders", std::string("2,1" ));
                tParameters.set( "bspline_pattern", std::string("0,1" ));

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

                tFieldHMR->evaluate_scalar_function( LevelSetFunction );

                for( uint k=0; k<2; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tFieldHMR );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tFieldHMR->evaluate_scalar_function( LevelSetFunction );
                }

                //-------- next mesh ---------------

                // manually select output pattern
                tHMR.get_database()->set_activation_pattern( 1 );

                for( uint k=0; k<2; ++k )
                {
                    for( uint Ik=0; Ik < tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc(); ++Ik )
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
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_Out = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                mtk::Integration_Mesh* tIntegrationMesh_Out =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_Out, tLagrangeMeshIndex );

                // Create integration first meshes
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh_In = tHMR.create_interpolation_mesh( tLagrangeMeshIndex_2 );
                mtk::Integration_Mesh* tIntegrationMesh_In =
                        create_integration_mesh_from_interpolation_mesh( MeshType::HMR, tInterpolationMesh_In, tLagrangeMeshIndex_2 );

                // Create mesh manager
                std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();
                uint tMeshIndex_Out = tMeshManager->register_mesh_pair( tInterpolationMesh_Out, tIntegrationMesh_Out );
                uint tMeshIndex_In  = tMeshManager->register_mesh_pair( tInterpolationMesh_In, tIntegrationMesh_In );

                mtk::Field_Proxy tField_In( tMeshManager, tMeshIndex_In, 0 );
                mtk::BSpline_Field tField_Out( tMeshManager, tMeshIndex_Out, 0 );

                tField_In.evaluate_scalar_function( LevelSetPlaneFunction );

                CHECK(equal_to( tField_In.get_nodal_values(tInterpolationMesh_In).numel(), 63));;

                // Use mapper
                mtk::Mapper tMapper(tMeshManager, tMeshIndex_In);
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tMapper.perform_mapping(
                        tField_Out,
                        EntityRank::BSPLINE,
                        EntityRank::NODE);

                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                mtk::Field_Proxy tField_Ref( tMeshManager, tMeshIndex_Out, 0 );
                tField_Ref.evaluate_scalar_function( LevelSetPlaneFunction );

                CHECK(equal_to( tField_Out.get_nodal_values().numel(), 465));

                for( uint Ik = 0; Ik < tField_Out.get_nodal_values().numel(); Ik++ )
                {
                    CHECK(equal_to( tField_Out.get_nodal_values()( Ik ), tField_Ref.get_nodal_values(tInterpolationMesh_Out)( Ik )));
                }
            }
        }


    }
}
