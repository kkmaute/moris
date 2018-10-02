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
#include "cl_HMR_Field.hpp"
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
//------------------------------------------------------------------------------

real
SphereFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

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

    // create a 2x2x2 mesh
    tParameters.set_number_of_elements_per_dimension( Matrix< DDLUMat >{ {1}, {1}, {1} } );

    // create mesh order 1 ( XTK does only support 1st order so far
    tParameters.set_mesh_order(2 );

    // we do not want tructation for this test
    // tParameters.set_bspline_truncation( false );

    // we do not need a buffer if there is no truncation
    // tParameters.set_buffer_size( 0 );

    // do not print debug information during test
    tParameters.set_verbose( true );

    tParameters.set_max_surface_level( 9 );

//------------------------------------------------------------------------------

    // create an HMR object and perform an arbitrary refinement
    HMR tHMR( tParameters );
    tHMR.save_mesh_to_vtk( "Primitive.vtk");
    // create mesh
    auto tMesh = tHMR.create_mesh();

    auto tField = tMesh->create_field("LevelSet");

    for( uint k=0; k<9; ++k )
    {
        // evaluate field
        tField->evaluate_scalar_function( SphereFunction );

        tHMR.flag_volume_and_surface_elements( tField );

        tHMR.perform_refinement();
        tHMR.update_refinement_pattern();
    }

    tHMR.finalize();

    // evaluate field one more time
    tField->evaluate_scalar_function( SphereFunction );

    tHMR.save_to_exodus("Mesh.exo");

//------------------------------------------------------------------------------

    // dump the Lagrange mesh
    tHMR.save_mesh_to_vtk( "LagrangeMesh.vtk");

    // dump the faces
    tHMR.save_faces_to_vtk( "Faces.vtk" );

    // dump the edges
    tHMR.save_edges_to_vtk( "Edges.vtk" );

    // dump the bsplines
    tHMR.save_edges_to_vtk( "BSplines.vtk" );

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
