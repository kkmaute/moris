/*
 * cl_MSI_Multigrid.cpp
 *
 *  Created on: Jul 14, 2018
 *      Author: schmidt
 */

#include "catch.hpp"
#include "typedefs.hpp"
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
#include "cl_FEM_Element.hpp"

#include "cl_MTK_Mapper.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "fn_r2.hpp"

moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

namespace moris
{
    namespace MSI
    {
    TEST_CASE("MSI_Multigrid","[MSI],[multigrid]")
    {
        if( moris::par_size() == 1 )
        {
            // order for this example
            moris::uint tOrder = 1;

            // create parameter object
            moris::hmr::Parameters tParameters;
            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );
            tParameters.set_verbose( false );
            tParameters.set_multigrid( true );
            tParameters.set_bspline_truncation( true );
            tParameters.set_mesh_orders_simple( tOrder );

            // create HMR object
            moris::hmr::HMR tHMR( tParameters );

            // flag first element for refinement
            tHMR.flag_element( 0 );
            tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );

            tHMR.flag_element( 0 );
            tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );

            tHMR.finalize();

            //tHMR.save_mesh_relations_to_hdf5_file( "Mesh_Dependencies_1.hdf5" );

             // grab pointer to output field
             std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );

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

             // tHMR.save_bsplines_to_vtk("BSplines.vtk");

             moris::map< moris::moris_id, moris::moris_index > tMap;
             tMesh->get_adof_map( tOrder, tMap );
             //tMap.print("Adof Map");

             //-------------------------------------------------------------------------------------------

             // create IWG object
             Cell< fem::IWG* > tIWGs ( 1, nullptr );
             tIWGs( 0 ) = new moris::fem::IWG_L2( );

             map< moris_id, moris_index >   tCoefficientsMap;
             Cell< fem::Node_Base* >        tNodes;
             Cell< MSI::Equation_Object* >  tElements;

             // get map from mesh
             tMesh->get_adof_map( tOrder, tCoefficientsMap );

             // ask mesh about number of nodes on proc
             luint tNumberOfNodes = tMesh->get_num_nodes();

             // create node objects
             tNodes.resize( tNumberOfNodes, nullptr );

             for( luint k = 0; k < tNumberOfNodes; ++k )
             {
                 tNodes( k ) = new fem::Node( &tMesh->get_mtk_vertex( k ) );
             }

             // ask mesh about number of elements on proc
             luint tNumberOfElements = tMesh->get_num_elems();

             // create equation objects
             tElements.resize( tNumberOfElements, nullptr );

             for( luint k=0; k<tNumberOfElements; ++k )
             {
                 // create the element
                 tElements( k ) = new fem::Element( & tMesh->get_mtk_cell( k ),
                                                    tIWGs,
                                                    tNodes );
             }

             MSI::Model_Solver_Interface * tMSI = new moris::MSI::Model_Solver_Interface( tElements,
                                                                                          tMesh->get_communication_table(),
                                                                                          tCoefficientsMap,
                                                                                          tMesh->get_num_coeffs( tOrder ),
                                                                                          tMesh.get() );

             tMSI->set_param("L2")= 1;

             tMSI->finalize( true );

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
             delete tIWGs( 0 );

             for( luint k=0; k<tNumberOfElements; ++k )
             {
                 // create the element
                 delete tElements( k );
             }

             for( luint k = 0; k < tNumberOfNodes; ++k )
             {
                 delete tNodes( k );
             }

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
//                 tParameters.set_mesh_orders_simple( tOrder );
//                 tParameters.set_verbose( true );
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
//                 moris::mapper::Mapper tMapper( tMesh );
//
//                 // call mapping function
//                 tMapper.perform_mapping(
//                         tField->get_label(),
//                         EntityRank::NODE,
//                         tField->get_label(),
//                         tField->get_bspline_rank() );
//
//                 tField->evaluate_node_values();
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
}


