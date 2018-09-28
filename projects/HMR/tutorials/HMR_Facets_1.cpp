//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

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
//------------------------------------------------------------------------------

/*!
 * /section Tutorial 3: MTK Mesh demo
 * This example creates a simple 2x2 mesh and tests MTK functionality.
 * This example will be turned into a test sool.
 */
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------

    /*!
     * This setup creates a minimal mesh and assumes an element edge length
     * of 1.
     *
     * HMR also supports parameter lists ( see tutorial )
     * For internal tests, using the parameter object is more convenient
     */
    Parameters tParameters;

    // create a 2x2 mesh
    tParameters.set_number_of_elements_per_dimension( Matrix< DDLUMat >{ {2}, {2} } );

    // create mesh order 1 ( XTK does only support 1st order so far
    tParameters.set_mesh_order( 1 );

    // we do not want tructation for this test
    tParameters.set_bspline_truncation( false );

    // we do not need a buffer if there is no truncation
    tParameters.set_buffer_size( 0 );

    // do not print debug information during test
    tParameters.set_verbose( false );

//------------------------------------------------------------------------------

    // create an HMR object and perform an arbitrary refinement
    HMR tHMR( tParameters );
    tHMR.flag_element( 0 );

    // the optional flag resets the pattern
    // ( if it is set, we HMR always starts with a tensor mesh )

    tHMR.perform_refinement( true );
    tHMR.flag_element( 1 );
    tHMR.perform_refinement( false );

//------------------------------------------------------------------------------

    // These are private functions that I use until MTK can write STK again
    tHMR.mLagrangeMeshes( 1 )->save_faces_to_vtk("Faces.vtk");

    if( tParameters.get_number_of_dimensions() == 3)
    {
        tHMR.mLagrangeMeshes( 1 )->save_edges_to_vtk("Edges.vtk");
    }
    tHMR.mLagrangeMeshes( 1 )->save_to_vtk("Elements.vtk");

//------------------------------------------------------------------------------
    // create mesh interface

    mtk::Mesh * tMesh = tHMR.create_output_mesh();

    if( par_size() == 1 )
    {
        Matrix< IndexMat > tIndices;

        // ----------  Test element neighbors
        // must return 3, 4, 5, 7, 8 in 2D

        tIndices = tMesh->get_elements_connected_to_element_loc_inds( 6 );

        print( tIndices, "Elements connected to Element 6" );

        // ----------  Test elements connected to node

        // must return 0, 3, 5, 6 in 2D
        tIndices = tMesh->get_elements_connected_to_node_loc_inds( 2 );
        print( tIndices, "Elements connected to Node 2" );

        // must return 6, 7, 8, 9 in 2D
        tIndices = tMesh->get_elements_connected_to_node_loc_inds( 13 );
        print( tIndices, "Elements connected to Node 13" );

        // ----------  Test faces connected to node

        // must return 28, 29, 30, 1 in 2D
        tIndices = tMesh->get_faces_connected_to_element_loc_inds( 7 );
        print( tIndices, "Faces connected to Element 7" );

        // ----------  Test elements connected to element

        // must return 3, 4, 5, 7, 8
        tIndices = tMesh->get_elements_connected_to_element_loc_inds( 6 );
        print( tIndices, "Elements connected to Element 6" );

    }
    delete tMesh;

    // tHMR.mBackgroundMesh->save_to_vtk("Background.vtk");

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
