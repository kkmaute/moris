/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Tutorial_2.cpp
 *
 */

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Logger.hpp"
#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer.hpp" // MRS/IOS/src

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

//------------------------------------------------------------------------------
// from MTK
//#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------

// HMR
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"

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
 * \section Tutorial_2: Load a mesh from a file
 *
 *  In this example we will
 *  - load a mesh from a file
 *  - load a field from a file
 *  - create two new fields
 *  - refine the mesh according to two fields
 *  - create an output mesh
 */
/*!
 * the following functions are used as level set for the tutorial
 *
 * \code{.cpp}
 * real
 * RosenbrockFunction( const Matrix< DDRMat > & aPoint )
 * {
 *     std::pow( 1.0 - aPoint( 0 ), 2 )
 *       + 100.0*std::pow( aPoint( 1 ) - aPoint( 0 )*aPoint( 0 ), 2 );
 * }
 * real
 * SimionescuFunction( const Matrix< DDRMat > & aPoint )
 * {
 *   return 0.1 * aPoint( 0 ) * aPoint ( 1 );
 * }
 * \endcode
 */
real
RosenbrockFunction( const Matrix< DDRMat > & aPoint )
{
    return std::pow( 1.0 - aPoint( 0 ), 2 )
     + 100.0*std::pow( aPoint( 1 ) - aPoint( 0 )*aPoint( 0 ), 2 );
}

real
SimionescuFunction( const Matrix< DDRMat > & aPoint )
{
    return 0.1 * aPoint( 0 ) * aPoint ( 1 );
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
//
//    /*!
//     * <b> Step 1: Load Mesh and Level Set Field from files </b>
//     */
//
//    /*!
//     * The Mesh and its parameters are loaded from a HDF5 Database.
//     *
//     * \code{.cpp}
//     * HMR tHMR( "Database.hdf5" );
//     * \endcode
//     */
//    HMR tHMR( "Database.hdf5" );
//
//    /*!
//     * The Circle Field is loaded from the HDF5 file
//     *
//     * \code{.cpp}
//     *  auto tCircle = tHMR.load_field_from_hdf5_file( "Circle", "Circle.hdf5" );
//     * \endcode
//     */
//    auto tCircle = tHMR.load_field_from_hdf5_file( "Circle", "Circle.hdf5" );
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 2: Create a New Field and evaluate a function on the loaded mesh</b>
//     *
//     * \code{.cpp}
//     *  auto tRosenbrock = tHMR.create_field( "Rosenbrock" );
//     *  tRosenbrock->evaluate_scalar_function( RosenbrockFunction );
//     * \endcode
//     */
//    auto tRosenbrock = tHMR.create_field( "Rosenbrock" );
//    tRosenbrock->evaluate_scalar_function( RosenbrockFunction );
//
//    auto tSimionescu = tHMR.create_field( "Simionescu" );
//    tSimionescu->evaluate_scalar_function( SimionescuFunction );
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 3: Refine the mesh again against the first field</b>
//     */
//
//    /*!
//     * Now we refine after the surfaces of both fields
//     * \code{.cpp}
//     * tHMR.flag_surface_elements( tSimionescu );
//     * tHMR.flag_surface_elements( tCircle );
//     * tHMR.perform_refinement_and_map_fields();
//     * \endcode
//     */
//
//     tHMR.flag_surface_elements( tSimionescu );
//     tHMR.flag_surface_elements( tCircle );
//
//     tHMR.perform_refinement_and_map_fields();
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 4: Save a second Exodus file</b>
//     *
//     * Investigate the refined mesh in paraview. Both fields are mapped
//     * onto the new mesh.
//     * \code{.cpp}
//     *  tHMR.save_to_exodus( "Mesh2.exo" );
//     * \endcode
//     */
//    tHMR.save_to_exodus( "Mesh2.exo" );
//
////------------------------------------------------------------------------------
//
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

