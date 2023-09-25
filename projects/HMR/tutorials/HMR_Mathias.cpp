/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Mathias.cpp
 *
 */

#include <memory>
#include <string>

// dynamik linker function
#include "dlfcn.h"

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp"

//------------------------------------------------------------------------------
// from LINALG
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_r2.hpp"
#include "fn_norm.hpp"
#include "HDF5_Tools.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_FEM_IWG_L2.hpp"

//#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

#include "../../FEM/INT/src/cl_FEM_Element_Bulk.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"

//------------------------------------------------------------------------------//
//HMR
#include "cl_HMR_Database.hpp"
#include "cl_MSI_Adof.hpp"
#include "cl_MSI_Dof_Manager.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Multigrid.hpp"

// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
//------------------------------------------------------------------------------

real
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}
//------------------------------------------------------------------------------
int
main(
        int    argc,
        char * argv[] )
{
//    // initialize MORIS global communication manager
//    gMorisComm = moris::Comm_Manager( &argc, &argv );
//
//    // Severity level 0 - all outputs
//    gLogger.initialize( 0 );
//
////------------------------------------------------------------------------------
//// this example creates a mesh and performs a manual refinement
////------------------------------------------------------------------------------
//    // order for this example
//    moris::uint tOrder = 1;
//
//    // create parameter object
//    moris::hmr::Parameters tParameters;
//    tParameters.set_number_of_elements_per_dimension( { { 1 }, { 1 } } );
//    tParameters.set_multigrid( true );
//    tParameters.set_bspline_truncation( true );
//
//    // create HMR object
//    moris::hmr::HMR tHMR( tParameters );
//
//    // flag first element for refinement
////    tHMR.flag_element( 0 );
////    tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
////
////    tHMR.flag_element( 0 );
////    tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
//
//    tHMR.finalize();
//
//     // grab pointer to output field
//     std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tOrder );
//
////             std::cout << std::endl;
////             // List of B-Spline Meshes
////             for( uint k=0; k<tHMR.get_database()->get_number_of_bspline_meshes(); ++k )
////             {
////                 // get pointer to mesh
////                 moris::hmr::BSpline_Mesh_Base * tMesh = tHMR.get_database()->get_bspline_mesh_by_index( k );
////
////                 std::cout << "BSpline Mesh " << k <<
////                         ": active pattern " << tMesh->get_activation_pattern() <<
////                         " order " << tMesh->get_order() <<
////                         " active basis " << tMesh->get_number_of_active_basis_on_proc()
////                         << std::endl;
////
////             }
////             std::cout << std::endl;
////
////             // List of Lagrange Meshes
////             for( uint k=0; k<tHMR.get_database()->get_number_of_lagrange_meshes(); ++k )
////             {
////                 // get pointer to mesh
////                 moris::hmr::Lagrange_Mesh_Base * tMesh = tHMR.get_database()->get_lagrange_mesh_by_index( k );
////
////                 std::cout << "Lagrange Mesh " << k <<
////                         ": active pattern " << tMesh->get_activation_pattern() <<
////                         " order " << tMesh->get_order() <<
////                         " active basis " << tMesh->get_number_of_nodes_on_proc() << std::endl;
////
////             }
////             std::cout << std::endl;
//
//     tHMR.save_bsplines_to_vtk("BSplines.vtk");
//
//     moris::map< moris::moris_id, moris::moris_index > tMap;
//     tMesh->get_adof_map( tOrder, tMap );
//     //tMap.print("Adof Map");
//
//     //-------------------------------------------------------------------------------------------
//
////      create IWG object
////     fem::IWG_L2 * tIWG = new moris::fem::IWG_L2( );
//
//     map< moris_id, moris_index >   tCoefficientsMap;
//     Cell< fem::Node_Base* >        tNodes;
//     Cell< MSI::Equation_Object* >  tElements;
//
//     // get map from mesh
//     tMesh->get_adof_map( tOrder, tCoefficientsMap );
//
//     // ask mesh about number of nodes on proc
//     luint tNumberOfNodes = tMesh->get_num_nodes();
//
//     // create node objects
//     tNodes.resize( tNumberOfNodes, nullptr );
//
//     for( luint k = 0; k < tNumberOfNodes; ++k )
//     {
//         tNodes( k ) = new fem::Node( &tMesh->get_mtk_vertex( k ) );
//     }
//
//     // ask mesh about number of elements on proc
//     luint tNumberOfElements = tMesh->get_num_elems();
//
//     // create equation objects
//     tElements.resize( tNumberOfElements, nullptr );
//
//
//     //fixme:This needs to be updated for integration/interpolation mesh
////     for( luint k=0; k<tNumberOfElements; ++k )
////     {
////         // create the element
////         tElements( k ) = new fem::Element( & tMesh->get_mtk_cell( k ),
////                                            tIWG,
////                                            tNodes );
////     }
////
////     MSI::Model_Solver_Interface * tMSI = new moris::MSI::Model_Solver_Interface( tElements,
////                                                                                  tMesh->get_communication_table(),
////                                                                                  tCoefficientsMap,
////                                                                                  tMesh->get_num_coeffs( tOrder ),
////                                                                                  tMesh.get() );
//
////     tMSI->set_param("L2")= (sint)tOrder;
//
////     tMSI->finalize();
//
////     std::cout << ( (int) ( tMSI != NULL ) ) << std::endl;
////------------------------------------------------------------------------------
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

