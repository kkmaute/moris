/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Multigrid.cpp
 *
 */

#include "catch.hpp"
#include "moris_typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"

#include "cl_Communication_Tools.hpp"
#include "cl_Communication_Manager.hpp"

#define protected public
#define private   public
#include "cl_MSI_Multigrid.hpp"
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Pdof_Host.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Node_Proxy.hpp"
#undef protected
#undef private

#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"

#include "cl_FEM_Node_Base.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Element_Bulk.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "fn_r2.hpp"

#include "fn_PRM_MSI_Parameters.hpp"

moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

namespace moris::MSI
{
    TEST_CASE("MSI_Multigrid","[MSI],[multigrid]")
    {
        if( moris::par_size() == 1 )
        {
            moris::uint tLagrangeMeshIndex = 0;
            moris::uint tBSplineMeshIndex = 0;

            // create parameter object
            moris::hmr::Parameters tParameters;
            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );

            tParameters.set_severity_level( 0 );
            tParameters.set_multigrid( true );
            tParameters.set_bspline_truncation( true );

            tParameters.set_output_meshes( {{ {0} }} );

            tParameters.set_lagrange_orders  ( { {1} });
            tParameters.set_lagrange_patterns({ {0} });

            tParameters.set_bspline_orders   ( { {1} } );
            tParameters.set_bspline_patterns ( { {0} } );

            tParameters.set_union_pattern( 2 );
            tParameters.set_working_pattern( 3 );

            tParameters.set_refinement_buffer( 1 );
            tParameters.set_staircase_buffer( 1 );

            // create HMR object
            moris::hmr::HMR tHMR( tParameters );

            // flag first element for refinement
            tHMR.flag_element( 0 );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tHMR.flag_element( 0 );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tHMR.finalize();

            //tHMR.save_mesh_relations_to_hdf5_file( "Mesh_Dependencies_1.hdf5" );

             // grab pointer to output field
             //std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );
             hmr::Interpolation_Mesh_HMR * tInterpolationMesh =  tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
             hmr::Integration_Mesh_HMR *   tIntegrationMesh =  tHMR.create_integration_mesh( 1, 0, tInterpolationMesh );

//             std::cout << std::endl;
//             // List of B-Spline Meshes
//             for( uint k=0; k<tHMR.get_database()->get_number_of_bspline_meshes(); ++k )
//             {
//                 // get pointer to mesh
//                 moris::hmr::BSpline_Mesh_Base * tMesh = tHMR.get_database()->get_bspline_mesh_by_index( k );
//
//                 std::cout << "BSpline Mesh " << k <<
//                         ": active pattern " << tMesh->get_activation_pattern() <<
//                         " order " << tMesh->get_order() <<
//                         " active basis " << tMesh->get_number_of_active_basis_on_proc()
//                         << std::endl;
//
//             }
//             std::cout << std::endl;
//
//             // List of Lagrange Meshes
//             for( uint k=0; k<tHMR.get_database()->get_number_of_lagrange_meshes(); ++k )
//             {
//                 // get pointer to mesh
//                 moris::hmr::Lagrange_Mesh_Base * tMesh = tHMR.get_database()->get_lagrange_mesh_by_index( k );
//
//                 std::cout << "Lagrange Mesh " << k <<
//                         ": active pattern " << tMesh->get_activation_pattern() <<
//                         " order " << tMesh->get_order() <<
//                         " active basis " << tMesh->get_number_of_nodes_on_proc() << std::endl;
//
//             }
//             std::cout << std::endl;

//              tHMR.save_bsplines_to_vtk("BSplines_for_matlab.vtk");

              //tHMR.save_mesh_relations_to_hdf5_file( "Mesh_Dependencies_matlab.hdf5" );

             //-------------------------------------------------------------------------------------------
             // create a L2 IWG
             fem::IWG_Factory tIWGFactory;
             std::shared_ptr< fem::IWG > tIWGL2 = tIWGFactory.create_IWG( fem::IWG_Type::L2 );
             tIWGL2->set_residual_dof_type( { { MSI::Dof_Type::L2 } } );
             tIWGL2->set_dof_type_list( {{ MSI::Dof_Type::L2 }}, mtk::Leader_Follower::LEADER );

             // define set info
              Vector< fem::Set_User_Info > tSetInfo( 1 );
//              tSetInfo( 0 ).set_mesh_index( 0 );
              tSetInfo( 0 ).set_mesh_set_name( "HMR_dummy_c_p0" );
              tSetInfo( 0 ).set_IWGs( { tIWGL2 } );

             map< moris_id, moris_index >   tCoefficientsMap;
             Vector< fem::Node_Base* >        tNodes;
             Vector< MSI::Equation_Object* >  tElements;

             // get map from mesh
             tInterpolationMesh->get_adof_map( tBSplineMeshIndex, tCoefficientsMap );

             // ask mesh about number of nodes on proc
             luint tNumberOfNodes = tInterpolationMesh->get_num_nodes();

             // create node objects
             tNodes.resize( tNumberOfNodes, nullptr );

             for( luint k = 0; k < tNumberOfNodes; ++k )
             {
                 tNodes( k ) = new fem::Node( &tInterpolationMesh->get_mtk_vertex( k ) );
             }

             Vector< MSI::Equation_Set * >      tElementBlocks(1,nullptr);

             // init the fem set counter
             moris::uint tFemSetCounter = 0;

//             uint tNumBlocks = tIntegrationMesh->get_num_blocks();

             // loop over the used mesh block-set
             for( luint Ik = 0; Ik < 1; ++Ik )
             {
                 // create a list of cell clusters (this needs to stay in scope somehow)
                 moris::mtk::Set * tBlockSet = tIntegrationMesh->get_set_by_index( 0 );

                 // create new fem set
                 tElementBlocks( Ik ) = new fem::Set( nullptr,
                                                      tBlockSet,
                                                      tSetInfo( Ik ),
                                                      tNodes );

                 // collect equation objects associated with the block-set
                 tElements.append( tElementBlocks( tFemSetCounter )->get_equation_object_list() );

                 // update fem set counter
                 tFemSetCounter++;
             }
             moris::Parameter_List tMSIParameters = prm::create_msi_parameter_list();
             tMSIParameters.set( "L2", 0 );
             tMSIParameters.set( "multigrid", true );

             MSI::Model_Solver_Interface * tMSI = new moris::MSI::Model_Solver_Interface( tMSIParameters,
                                                                                          tElementBlocks,
                                                                                          tInterpolationMesh->get_communication_table(),
                                                                                          tCoefficientsMap,
                                                                                          tInterpolationMesh->get_max_num_coeffs_on_proc( tBSplineMeshIndex ),
                                                                                          tInterpolationMesh );

             tElementBlocks( 0 )->finalize( tMSI );

             tMSI->finalize();

             moris::Matrix< DDSMat > tExternalIndices( 9, 1 );
             tExternalIndices( 0, 0 ) = 17;
             tExternalIndices( 1, 0 ) = 18;
             tExternalIndices( 2, 0 ) = 21;
             tExternalIndices( 3, 0 ) = 22;
             tExternalIndices( 4, 0 ) = 23;
             tExternalIndices( 5, 0 ) = 25;
             tExternalIndices( 6, 0 ) = 27;
             tExternalIndices( 7, 0 ) = 26;
             tExternalIndices( 8, 0 ) = 28;

             moris::Matrix< DDSMat > tInternalIndices;

             tMSI->read_multigrid_maps( 2, tExternalIndices, 0, tInternalIndices );

             CHECK( equal_to( tInternalIndices( 0, 0 ), 0 ) );
             CHECK( equal_to( tInternalIndices( 1, 0 ), 1 ) );
             CHECK( equal_to( tInternalIndices( 2, 0 ), 2 ) );
             CHECK( equal_to( tInternalIndices( 3, 0 ), 3 ) );
             CHECK( equal_to( tInternalIndices( 4, 0 ), 4 ) );
             CHECK( equal_to( tInternalIndices( 5, 0 ), 5 ) );
             CHECK( equal_to( tInternalIndices( 6, 0 ), 6 ) );
             CHECK( equal_to( tInternalIndices( 7, 0 ), 7 ) );
             CHECK( equal_to( tInternalIndices( 8, 0 ), 8 ) );

             delete tMSI;
             //delete tIWGs( 0 );

             for( luint k = 0; k < tNumberOfNodes; ++k )
             {
                 delete tNodes( k );
             }
             delete tIntegrationMesh;
             delete tInterpolationMesh;

        }
    }

//        TEST_CASE("MSI_Multigrid1","[MSI],[multigrid1]")
//        {
//            if( moris::par_size() == 1 )
//            {
//                 moris::hmr::Parameters tParameters;
//
//                 tParameters.set_number_of_elements_per_dimension( { { 2} , { 2 } } );
//                 tParameters.set_domain_offset( 0, 0 );
//
//                 uint tOrder = 2;
//
//                 tParameters.set_multigrid( true );
//
//                 // create HMR object
//                 moris::hmr::HMR tHMR( tParameters );
//
//                 // flag first element
//                 tHMR.flag_element( 0 );
//                 tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
//
//                 // finish mesh
//                 tHMR.finalize();
//
//                 // grab pointer to output field
//                 std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );
//
//                 // create field
//                 std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tOrder );
//                 std::shared_ptr< moris::hmr::Field > tExact = tMesh->create_field( "Exact", tOrder );
//
//                 // evaluate node values
//                 tField->evaluate_scalar_function( LevelSetFunction );
//                 tExact->get_node_values() = tField->get_node_values();
//
//                 // create mapper
//                 moris::mtk::Mapper tMapper( tMesh );
//
//                 // call mapping function
//                 tMapper.perform_mapping(
//                         tField->get_label(),
//                         mtk::EntityRank::NODE,
//                         tField->get_label(),
//                         tField->get_bspline_rank() );
//
//                 tField->evaluate_nodal_values();
//
//                 // save field to hdf5
//                 tField->save_field_to_hdf5("Circle.hdf5");
//
//                 // determine coefficient of determination
//                 moris::real tR2 = moris::r2( tExact->get_node_values(),
//                         tField->get_node_values() );
//
//                 std::cout << "R2 " << tR2 << std::endl;
//            }
//    }
    }
