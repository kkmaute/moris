/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */

// MORIS header files.
#include <fn_r2.hpp>
#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_HMR.hpp" // HMR/src
#include "cl_FEM_Element.hpp" // FEM/INT/src
#include "cl_Model_Solver_Interface.hpp"
#include "cl_Equation_Object.hpp"
#include "cl_FEM_IWG_L2_Test.hpp"

moris::Comm_Manager gMorisComm;

using namespace moris;

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);

    // print welcome banner and system information
    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
    // In this part, an HMR object is created

    // create settings object
    /* ( this is actually my internal settings object. As far as I understood,
                 the user is supposed to control HMR using a parameter list.
                 This is not implemented yet ) */

    hmr::Parameters tParameters;

    // create a Mat for a 2D object
    Mat< luint > tNumberOfElements = { { 2 }, { 2 } };

    // mesh orders
    Mat< uint >  tMeshOrders = { { 1 }, { 1 }, { 1 } };

    // set interpolation degrees
    tParameters.set_mesh_orders( tMeshOrders );

    // pass number of elements to settings
    tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

    // make mesh output silent
    tParameters.set_verbose( true );

    // buffer size must be set at least to max polynomial if truncation is used
    //tParameters.set_buffer_size( tParameters.get_max_polynomial() );

    tParameters.set_buffer_size( 0 );
    tParameters.set_bspline_truncation( false );

    hmr::HMR tHMR( &tParameters );

    tHMR.set_active_pattern( 1 );
    tHMR.flag_element( 0 );
    tHMR.perform_refinement();

    tHMR.set_active_pattern( 2 );
    tHMR.flag_element( 3 );
    tHMR.perform_refinement();

    tHMR.unite_patterns( 1, 2, 0 );

    tHMR.update_meshes();
//------------------------------------------------------------------------------

    tHMR.set_active_pattern( 1 );


    //tHMR.flag_element( 3 );
    //tHMR.perform_refinement();

    // get number of vertices
    tHMR.save_to_exodus( 1, "Mesh.exo" );


//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
