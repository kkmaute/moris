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
 * the following function is used as level set for the tutorial
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
 * \section Tutorial 1: Initialize a Simple 2D Mesh
 * In this example, we will
 *
 * - create a simple 2d mesh
 * - create a function
 * - save the mesh to a file
 * - save the function to a file
 *
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
      * <b> Step 1: create a parameter list </b>
      */

    /*!
     * The parameter list controls settings that are used by HMR, such
     * as the setup of the background mesh, polynomial degree et cetera.
     * The following function creates a default list
     * \code{.cpp}
     * ParameterList tParameters = create_parameter_list();
     * \endcode
     */
    ParameterList tParameters = create_parameter_list();

    /*!
       * the coarsest grid will have 3 elements in x-direction,
       * and 2 elements in y-direction
       *
       * \code{.cpp}
       * tParameters.set("number_of_elements_per_dimension", "3, 2" );
       * \endcode
       */
      tParameters.set("number_of_elements_per_dimension", "3, 2" );


      /*!
       * the domain will have a width of 3 and a height of 2 units
       *
       * \code{.cpp}
       * tParameters.set("domain_dimensions", "3, 2" );
       * \endcode
       */
      tParameters.set("domain_dimensions", "3, 2" );

      /*!
       * set the coordinates of the lower left node on the mesh
       *
       * \code{.cpp}
       * tParameters.set("domain_offset", "0, 0");
       * \endcode
       */
      tParameters.set("domain_offset", "0, 0");

//------------------------------------------------------------------------------

      /*!
       * <b> Step 2: HMR object </b>
       */

      /*!
       * All operations such as refining a mesh according to a field,
       * mapping a field onto a new mesh and providing the API of the new mesh
       * are handled by the HMR object.
       *
       * \code{.cpp}
       * HMR tHMR( tParameters );
       * \endcode
       */
      HMR tHMR( tParameters );


//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
