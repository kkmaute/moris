/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Facets_1.cpp
 *
 */

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Field.hpp"
//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// HMR
#include "cl_HMR_Parameters.hpp"
#define private public
#define protected public
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"
#undef private
#undef protected

//------------------------------------------------------------------------------

// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
//------------------------------------------------------------------------------

/*!
 * \section Facets Tutorial
 * This example creates a simple 2x2 mesh and tests MTK functionality.
 * This example was written for Keenan and will be turned into a test soon.
 */
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
//
//    /*!
//     * This setup creates a minimal mesh and assumes an element edge length
//     * of 1.
//     *
//     * HMR also supports parameter lists ( see tutorial )
//     * For internal tests, using the parameter object is more convenient
//     */
//    Parameters tParameters;
//
//    // create a 2x2 mesh
//    tParameters.set_number_of_elements_per_dimension( Matrix< DDLUMat >{ {2}, {2}, {2} } );
//
//    // create mesh order 1 ( XTK does only support 1st order so far
//    tParameters.set_mesh_order( 1 );
//
//    // we do not want tructation for this test
//    tParameters.set_bspline_truncation( false );
//
//    // we do not need a buffer if there is no truncation
//    tParameters->set_refinement_buffer( 0 );
//    tParameters->set_staircase_buffer( 0 );
//
////------------------------------------------------------------------------------
//
//    // create an HMR object and perform an arbitrary refinement
//    HMR tHMR( tParameters );
//    tHMR.flag_element( 0 );
//
//    // the optional flag resets the pattern
//    // ( if it is set, we HMR always starts with a tensor mesh )
//
//    tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
//    tHMR.flag_element( 1 );
//    tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
//
////------------------------------------------------------------------------------
//
//    // dump the Lagrange mesh
//    tHMR.save_mesh_to_vtk( "LagrangeMesh.vtk");
//
//    // dump the faces
//    tHMR.save_faces_to_vtk( "Faces.vtk" );
//
//    // dump the edges
//    tHMR.save_edges_to_vtk( "Edges.vtk" );
//
////------------------------------------------------------------------------------
//    // create mesh interface
//
//    auto tMesh = tHMR.create_mesh();
//
//    if( par_size() == 1 )
//    {
//        Matrix< IndexMat > tIndices;
//
//        if( tParameters.get_number_of_dimensions() == 2)
//        {
//            // ----------  test element neighbors
//            // must return 3, 4, 5, 7, 8 in 2D
//
//            tIndices = tMesh->get_elements_connected_to_element_loc_inds( 6 );
//
//            print( tIndices, "Elements connected to Element 6" );
//
//            // ----------  test elements connected to node
//
//            // must return 0, 3, 5, 6 in 2D
//            tIndices = tMesh->get_elements_connected_to_node_loc_inds( 2 );
//            print( tIndices, "Elements connected to Node 2" );
//
//            // must return 6, 7, 8, 9 in 2D
//            tIndices = tMesh->get_elements_connected_to_node_loc_inds( 13 );
//            print( tIndices, "Elements connected to Node 13" );
//
//            // ----------  test faces connected to node
//
//            // must return 28, 29, 30, 1 in 2D
//            tIndices = tMesh->get_faces_connected_to_element_loc_inds( 7 );
//            print( tIndices, "Faces connected to Element 7" );
//
//            // ----------  test nodes connected to element
//        }
//        else if( tParameters.get_number_of_dimensions() == 3 )
//        {
//            // ----------  test element neighbors
//            // must return 5, 6, 7, 8, 11, 14, 15, 18
//            tIndices = tMesh->get_elements_connected_to_element_loc_inds( 12 );
//            print( tIndices, "Elements connected to Element 12" );
//
//            // ----------  test elements connected to node
//            // must return 14, 15, 16, 17, 18, 19, 20, 21
//            tIndices = tMesh->get_elements_connected_to_node_loc_inds( 42 );
//            print( tIndices, "Elements connected to Node 42" );
//
//            // must return 12, 14
//            tIndices = tMesh->get_elements_connected_to_node_loc_inds( 45 );
//            print( tIndices, "Elements connected to Node 45" );
//
//            // ----------  test faces connected to node
//
//            // must return 68, 69, 70, 75, 77
//            tIndices = tMesh->get_faces_connected_to_node_loc_inds( 42 );
//            print( tIndices, "Faces connected to Node 42" );
//
//            // ----------  test nodes connected to element
//            // must return 6. 30, 36, 33, 39, 42, 45, 43
//            tIndices = tMesh->get_nodes_connected_to_element_loc_inds( 14 );
//            print( tIndices, "Nodes connected to element 14" );
//
//            //  ---------- test edges connected to element
//            // 31, 97, 98, 91, 101, 108, 117, 112, 111, 118, 119, 114
//            tIndices = tMesh->get_edges_connected_to_element_loc_inds( 14 );
//            print( tIndices, "Edges connected to element 14" );
//
//            // ------------ test faces connected to element
//            // 79, 75, 76, 71, 61, 77
//            tIndices = tMesh->get_faces_connected_to_element_loc_inds( 14 );
//            print( tIndices, "Faces connected to element 14" );
//        }
//    }
//
////------------------------------------------------------------------------------
//
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

