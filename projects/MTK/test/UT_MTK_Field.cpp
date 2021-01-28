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

#define protected public
#define private   public
#include "cl_MTK_Field_Proxy.hpp"
#include "cl_MTK_Mapper.hpp"
#undef protected
#undef private


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
    LevelSetPlainFunction( const moris::Matrix< moris::DDRMat > & aPoint )
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


                mtk::Field * tField_In = new mtk::Field_Proxy( tMeshManager, tMeshIndex_In );
                mtk::Field * tField_Out = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out );

                reinterpret_cast< mtk::Field_Proxy* >( tField_In )->evaluate_scalar_function( LevelSetPlainFunction );
                reinterpret_cast< mtk::Field_Proxy* >( tField_Out )->evaluate_scalar_function( LevelSetPlainFunction );

                std::cout<<"Field_In size: "<<tField_In->get_node_values().numel()<<std::endl;
                std::cout<<"Field_Out size: "<<tField_Out->get_node_values().numel()<<std::endl;

                CHECK(equal_to( tField_In->get_node_values().numel(), 125));
                CHECK(equal_to( tField_Out->get_node_values().numel(), 46));

                //                tFieldHMR->get_node_values() = tField_In->get_node_values();
                //                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );
                //                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                delete tField_In;
                delete tField_Out;
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

                mtk::Field * tField_In = new mtk::Field_Proxy( tMeshManager, tMeshIndex_In, 0 );
                mtk::Field * tField_Out = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out, 0 );

                reinterpret_cast< mtk::Field_Proxy* >( tField_In )->evaluate_scalar_function( LevelSetPlainFunction );

                CHECK(equal_to( tField_In->get_node_values().numel(), 46));;

                // Use mapper
                mapper::Mapper tMapper;
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tField_Out->evaluate_node_values();

                tFieldHMR->get_node_values() = tField_Out->get_node_values();

                //tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                //tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                mtk::Field * tField_Ref = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out, 0 );
                reinterpret_cast< mtk::Field_Proxy* >( tField_Ref )->evaluate_scalar_function( LevelSetPlainFunction );

                CHECK(equal_to( tField_Out->get_node_values().numel(), 125));

                for( uint Ik = 0; Ik < tField_Out->get_node_values().numel(); Ik++ )
                {
                    CHECK(equal_to( tField_Out->get_node_values()( Ik ), tField_Ref->get_node_values()( Ik )));
                }

                delete tField_In;
                delete tField_Out;
                delete tField_Ref;
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

                mtk::Field * tField_In = new mtk::Field_Proxy( tMeshManager, tMeshIndex_In, 0 );
                mtk::Field * tField_Out = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out, 0 );

                reinterpret_cast< mtk::Field_Proxy* >( tField_In )->evaluate_scalar_function( LevelSetPlainFunction );

                CHECK(equal_to( tField_In->get_node_values().numel(), 217));;

                // Use mapper
                mapper::Mapper tMapper;
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tField_Out->evaluate_node_values();

                tFieldHMR->get_node_values() = tField_Out->get_node_values();

                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                mtk::Field * tField_Ref = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out, 0 );
                reinterpret_cast< mtk::Field_Proxy* >( tField_Ref )->evaluate_scalar_function( LevelSetPlainFunction );

                CHECK(equal_to( tField_Out->get_node_values().numel(), 133));

                for( uint Ik = 0; Ik < tField_Out->get_node_values().numel(); Ik++ )
                {
                    CHECK(equal_to( tField_Out->get_node_values()( Ik ), tField_Ref->get_node_values()( Ik )));
                }

                delete tField_In;
                delete tField_Out;
                delete tField_Ref;
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

                mtk::Field * tField_In = new mtk::Field_Proxy( tMeshManager, tMeshIndex_In, 0 );
                mtk::Field * tField_Out = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out, 0 );

                reinterpret_cast< mtk::Field_Proxy* >( tField_In )->evaluate_scalar_function( LevelSetPlainFunction );

                CHECK(equal_to( tField_In->get_node_values().numel(), 63));;

                // Use mapper
                mapper::Mapper tMapper;
                tMapper.map_input_field_to_output_field( tField_In, tField_Out );

                tField_Out->evaluate_node_values();

                tFieldHMR->get_node_values() = tField_Out->get_node_values();

                tHMR.save_to_exodus( 0, "./mtk_field_test.e" );

                tHMR.save_to_exodus( 1, "./mtk_field_test_1.e" );

                mtk::Field * tField_Ref = new mtk::Field_Proxy( tMeshManager, tMeshIndex_Out, 0 );
                reinterpret_cast< mtk::Field_Proxy* >( tField_Ref )->evaluate_scalar_function( LevelSetPlainFunction );

                CHECK(equal_to( tField_Out->get_node_values().numel(), 465));

                for( uint Ik = 0; Ik < tField_Out->get_node_values().numel(); Ik++ )
                {
                    CHECK(equal_to( tField_Out->get_node_values()( Ik ), tField_Ref->get_node_values()( Ik )));
                }

                delete tField_In;
                delete tField_Out;
                delete tField_Ref;
            }
        }


    }
}
