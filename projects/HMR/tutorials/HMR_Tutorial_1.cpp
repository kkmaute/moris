/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Tutorial_1.cpp
 *
 */

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "moris_typedefs.hpp"
#include "cl_Logger.hpp"
#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer.hpp" // MRS/IOS/src

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

//------------------------------------------------------------------------------

// geometry engine
//#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------

// HMR
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

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
CircleFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

//------------------------------------------------------------------------------

/*!
 * \section Tutorial_1: Initialize a Simple 2D Mesh
 * In this example, we will
 *
 * - create a simple 2d mesh
 * - create a nodal field
 * - refine the mesh with respect to that field
 * - save the mesh and field data to files
 *
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
//      * <b> Step 1: create a parameter list </b>
//      */
//
//    /*!
//     * The parameter list controls settings that are used by HMR, such
//     * as the setup of the background mesh, polynomial degree et cetera.
//     * The following function creates a default list
//     * \code{.cpp}
//     * ParameterList tParameters = create_hmr_parameter_list();
//     * \endcode
//     */
//    ParameterList tParameters = create_hmr_parameter_list();
//
//    /*!
//       * the coarsest grid will have 10 elements in x-direction,
//       * and 10 elements in y-direction
//       *
//       * \code{.cpp}
//       * tParameters.set("number_of_elements_per_dimension", "9, 9" );
//       * \endcode
//       */
//      tParameters.set( "number_of_elements_per_dimension", "9, 9" );
//
//
//      /*!
//       * the domain will have a width of 3 and a height of 3 units
//       *
//       * \code{.cpp}
//       * tParameters.set( "domain_dimensions", "3, 3" );
//       * \endcode
//       */
//      tParameters.set( "domain_dimensions", "3, 3" );
//
//      /*!
//       * set the coordinates of the lower left node on the mesh
//       *
//       * \code{.cpp}
//       * tParameters.set( "domain_offset", "-1.5, -1.5" );
//       * \endcode
//       */
//      tParameters.set( "domain_offset", "-1.5, -1.5" );
//
//
//      /*!
//       * By default, HMR is silent.
//       * If this output is set, HMR will print lots of debug information
//       * on screen.
//       *
//       * \code{.cpp}
//       * tParameters.set( "verbose", 1 );
//       * \endcode
//       */
//      tParameters.set( "truncate_bsplines", 1 );
//
//      tParameters.set( "bspline_orders", "2" );
//      tParameters.set( "lagrange_orders", "2" );
//
////------------------------------------------------------------------------------
//
//      /*!
//       * <b> Step 2: HMR object </b>
//       */
//
//      /*!
//       * All operations such as refining a mesh according to a field,
//       * mapping a field onto a new mesh and providing the API of the new mesh
//       * are handled by the HMR object.
//       *
//       * \code{.cpp}
//       * HMR tHMR( tParameters );
//       * \endcode
//       */
//      HMR tHMR( tParameters );
//
////------------------------------------------------------------------------------
//
//      /*!
//       * <b> Step 3: Creating a nodal field and refining according to it</b>
//       */
//
//      /*!
//       * The following command creates a shared pointer to a field that is
//       * called "LevelSet". The datatype is std::shared_ptr<moris::hmr::Field>
//       *
//       * \code{.cpp}
//       * auto tField = tHMR.create_field( "Circle" );
//       * \endcode
//       */
//      auto tField = tHMR.create_field( "Circle" );
//
//      /*!
//       * This example uses an analytic level set, which is defined as follows
//       *
//       * \code{.cpp}
//       * real
//       * CircleFunction( const Matrix< DDRMat > & aPoint )
//       * {
//       *     return norm( aPoint ) - 1.2;
//       * }
//       * \endcode
//       *
//       * The pointer of this function is passed to the field.
//       *
//       * \code{.cpp}
//       * tField->evaluate_scalar_function( CircleFunction );
//       * \endcode
//       */
//      tField->evaluate_scalar_function( CircleFunction );
//
//      /*!
//       * In the next step, we use this field to identify elements that
//       * are fully inside the level set, or intersected.
//       *
//       * This flagging can be repeated with an arbitrary number
//       * of fields.
//       *
//       * \code{.cpp}
//       * tHMR.flag_volume_and_surface_elements( tField );
//       * \endcode
//       */
//      tHMR.flag_volume_and_surface_elements_on_working_pattern( tField );
//
//      /*!
//       * One all elements are flagged for the refinement, a procedure is
//       * called which performs one refinement step and maps all fields
//       * to the new mesh.
//       *
//       * \code{.cpp}
//       * tHMR.perform_refinement_and_map_fields();
//       * \endcode
//       */
//       tHMR.perform_refinement_and_map_fields();
//
////------------------------------------------------------------------------------
//
//      /*!
//       * <b> Step4: Saving the data</b>
//       */
//
//      /*!
//       * This command creates an exodus file of the refined mesh.
//       *
//       * \code{.cpp}
//       * tHMR.save_to_exodus( "Mesh1.exo" );
//       * \endcode
//       */
//      tHMR.save_to_exodus( "Mesh1.exo" );
//
//      /*!
//       * One can also recall the state of the mesh before the last
//       * refinement, and save it into a file.
//       *
//       * \code{.cpp}
//       * tHMR.save_last_step_to_exodus( "LastStep.exo" );
//       * \endcode
//       */
//      tHMR.save_last_step_to_exodus( "LastStep.exo" );
//
//      /*!
//       * The database that defines the mesh in its refined state saved using this command.
//       * Note that field information is not stored.
//       *
//       * \code{.cpp}
//       * tHMR.save_to_hdf5( "Database.hdf5" );
//       * \endcode
//       */
//      tHMR.save_to_hdf5( "Database.hdf5" );
//
//      /*!
//       * For third party applications, the interpolation coefficients
//       * can be saved into an output file as well
//       * \code{.cpp}
//       * tHMR.save_coeffs_to_hdf5_file( "TMatrix.hdf5" );
//       * \endcode
//       */
//      tHMR.save_coeffs_to_hdf5_file( "TMatrix.hdf5" );
//
//      /*!
//       * Each individual field can be stored into an hdf5 file on its own.
//       * These files can be exchanged with other codes such as MATLAB
//       * \code{.cpp}
//       * tField->save_field_to_hdf5("Circle.hdf5");
//       * \endcode
//       */
//      tField->save_field_to_hdf5( "Circle.hdf5" );
//
////------------------------------------------------------------------------------
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

