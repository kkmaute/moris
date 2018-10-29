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
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------
#include "fn_r2.hpp"
#include "fn_norm.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
// select namespaces
using namespace moris;
using namespace hmr;

real
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}


//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
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
     * ParameterList tParameters = create_hmr_parameter_list();
     * \endcode
     */
    ParameterList tParameters = create_hmr_parameter_list();

      tParameters.set( "number_of_elements_per_dimension", "4, 4, 4" );

      tParameters.set( "domain_offset", "-2, -2, -2" );
      tParameters.set( "domain_dimensions", "4, 4, 4" );
      uint tOrder = 2;
      tParameters.set( "interpolation_order", "2" );
      tParameters.set( "verbose", 1 );

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

      /*!
       * <b> Step 3: Creating a nodal field and refining according to it</b>
       */

      /*!
       * The following command creates a shared pointer to a field that is
       * called "LevelSet". The datatype is std::shared_ptr<moris::hmr::Field>
       *
       * \code{.cpp}
       * auto tField = tHMR.create_field( "Circle" );
       * \endcode
       */
       auto tField = tHMR.create_field( "Circle", tOrder, tOrder );

      /*!
       * This example uses an analytic level set, which is defined as follows
       *
       * \code{.cpp}
       * real
       * CircleFunction( const Matrix< DDRMat > & aPoint )
       * {
       *     return norm( aPoint ) - 1.2;
       * }
       * \endcode
       *
       * The pointer of this function is passed to the field.
       *
       * \code{.cpp}
       * tField->evaluate_scalar_function( CircleFunction );
       * \endcode
       */
      tField->evaluate_scalar_function( CircleFunction );

      /*!
       * In the next step, we use this field to identify elements that
       * are fully inside the level set, or intersected.
       *
       * This flagging can be repeated with an arbitrary number
       * of fields.
       *
       * \code{.cpp}
       * tHMR.flag_volume_and_surface_elements( tField );
       * \endcode
       */
       tHMR.flag_volume_and_surface_elements( tField );

      /*!
       * One all elements are flagged for the refinement, a procedure is
       * called which performs one refinement step and maps all fields
       * to the new mesh.
       *
       * \code{.cpp}
       * tHMR.perform_refinement_and_map_fields();
       * \endcode
       */
       tHMR.perform_refinement_and_map_fields();

//------------------------------------------------------------------------------


      auto tMesh = tHMR.create_mesh();

      auto tExact = tMesh->create_field( "Exact", tField->get_bspline_order() );

      tExact->evaluate_scalar_function( CircleFunction );


//      /print( tField->get_node_values(), "L2" );
      moris::real tR2 = moris::r2(
                              tExact->get_node_values(),
                              tField->get_node_values() );

                      std::cout << "R2 " << tR2 << std::endl;

                      tHMR.save_to_exodus( "Mesh1.exo" );
//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
