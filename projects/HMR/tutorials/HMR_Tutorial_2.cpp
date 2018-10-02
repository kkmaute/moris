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
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

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
 * <h1>Tutorial 2 : Load 2D Mesh from File and refine again"</h1>
 * In the following example, a mesh and a nodal field are loaded
 * from existing files. The refinement is performed and
 * the new mesh is stored into an output.
 *
 * Note that the files created by Tutorial 1 must exist on the file system.
 */
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
