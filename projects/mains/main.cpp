/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */

// MORIS header files.

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_HMR.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "fn_norm.hpp"
#include "fn_r2.hpp"

moris::Comm_Manager gMorisComm;

using namespace moris;

real
lvlset( const Mat< real > & aPoint )
{
    /*const real tAlpha = 5; // 8.5
    const real tBeta =  0.5; // 0.625
    const real tDelta = 0.5;    // 1.0
    const real tEpsilon = 0.05;
    const real tOmega = 8*3.141592653589793; */
    //real tX = aPoint( 0 );


    //real tX2 = 1.0;

    /* if ( tX < tX1 )
    {
        return 0.0;
    }
    else if ( tX < tX2 )
    {
        return std::cos( ( tX - tX0 ) * tOmega ) + 1;
    }
    else
    {
        return 0.0;
    } */


    //return tX * tX ;
    /*real tX0 = 0.625;
    real tX1 = 0.375;
    real tX2 = 0.875;

    //



    if( tX < tX1 )
    {
        return -1.0;
    }
    else if ( tX < tX2 )
    {
        return std::cos( ( tX - tX0 )*4*3.141592653589793 );
    }
    else
    {
        return -1.0;
    } */
    //return 1.0/( 1.0 + std::exp( - tAlpha* ( tX - tBeta ) ) ) - tDelta
   //         + tEpsilon * std::sin( tX*tOmega );


    return norm( aPoint ) - 0.5;
}

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // print welcome banner and system information
    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
//  HMR Parameters setup
//------------------------------------------------------------------------------

    // The parameter object controls the behavior of HMR.
    hmr::Parameters tParameters;

    //uint tElements = 0;
    //std::cin >> tElements;
    // We create a 2-Dimensional mesh with 2x2 elements ...
    tParameters.set_number_of_elements_per_dimension( 2, 2 );

    tParameters.set_max_volume_level( 3 );
    tParameters.set_max_surface_level( 4 );

    // we create a square of 2x2
    tParameters.set_domain_dimensions( 1.0, 1.0 );
    //tParameters.set_domain_offset( -3.0, -3.0 );

    // set default mesh order to two
    tParameters.set_mesh_order( 1 );
    //tParameters.set_mesh_orders_simple( 1 );
    // make mesh output silent
    tParameters.set_verbose( true );

    // B-Spline truncation is turned on by default. It is recommended to leave this
    // setting as is.
    tParameters.set_bspline_truncation( true );
    //tParameters.set_buffer_size( 2 );

//------------------------------------------------------------------------------
//  HMR Initialization
//------------------------------------------------------------------------------

    // We create the HMR object by passing the settings to the constructor
    hmr::HMR tHMR( tParameters );

    for( uint k=0; k<6; ++k )
    {
        auto tField = tHMR.create_field( "Circle", 1 );
        tField->evaluate_function( lvlset );

        tHMR.refine_against_nodal_field( tField->get_data() );
        delete tField;

    }
    auto tField = tHMR.create_field( "Circle", 1 );

    tField->evaluate_function( lvlset );
    tHMR.add_field( tField );
    // save mesh to file
    //tHMR.save_to_exodus( 0, "Mesh0.exo" );
    tHMR.save_to_exodus( 1, "Mesh1.exo" );
    //tHMR.save_to_exodus( 2, "Mesh2.exo" );
    //tHMR.save_to_exodus( 3, "Mesh3.exo" ); */
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
