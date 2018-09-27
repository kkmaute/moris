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

//------------------------------------------------------------------------------
// from MTK
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
 * the folowing function is used as level set for the tutorial
 *
 * \code{.cpp}
 * real
 * LevelSetFunction( const Matrix< DDRMat > & aPoint )
 * {
 *     return norm( aPoint ) - 1.2;
 * }
 * \endcode
 */
real
LevelSetFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

//------------------------------------------------------------------------------

/*!
 * /section Tutorial 1: Initialize a Simple 2D Mesh
 * This example creates a simple 2x3 mesh in serial and refines it according
 * to a given LevelSet.
 */
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------

    Parameters tParameters;

    tParameters.set_number_of_elements_per_dimension( Matrix< DDLUMat >{ {2}, {2}, {2} } );

    tParameters.set_mesh_order( 1 );

    tParameters.set_bspline_truncation( false );

    tParameters.set_buffer_size( 0 );

    // do not print debug information during test
    tParameters.set_verbose( false );

//------------------------------------------------------------------------------

    /*!
     * <b> Step 2: HMR object </b>
     */

    /*!
     * create the HMR main object from the parameters
     *
     * \code{.cpp}
     * HMR tHMR( tParameters );
     * \endcode
     */
    HMR tHMR( tParameters );
    tHMR.flag_element( 0 );
    tHMR.perform_refinement( true );
    //tHMR.flag_element( 1 );
    //tHMR.perform_refinement( false );

    tHMR.mLagrangeMeshes( 1 )->save_faces_to_vtk("Faces.vtk");
    tHMR.mLagrangeMeshes( 1 )->save_edges_to_vtk("Edges.vtk");
    tHMR.mLagrangeMeshes( 1 )->save_to_vtk("Elements.vtk");
    // tHMR.mBackgroundMesh->save_to_vtk("Background.vtk");

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
