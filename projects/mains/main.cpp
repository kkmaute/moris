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
#define private public
#define protected public
#include "cl_HMR.hpp"
#undef private
#undef public
#include "cl_MDL_Model.hpp"
#include "cl_FEM_IWG_L2.hpp"

#include "fn_norm.hpp"
#include "fn_r2.hpp"

moris::Comm_Manager gMorisComm;

using namespace moris;

real
lvlset( const Mat< real > & aPoint )
{
    const real tAlpha = 5; // 8.5
    const real tBeta =  0.5; // 0.625
    const real tDelta = 0.5;    // 1.0
    const real tEpsilon = 0.05;
    const real tOmega = 8*3.141592653589793;
    real tX = aPoint( 0 );


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
    return 1.0/( 1.0 + std::exp( - tAlpha* ( tX - tBeta ) ) ) - tDelta
            + tEpsilon * std::sin( tX*tOmega );



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

    uint tElements = 0;
    std::cin >> tElements;
    // We create a 2-Dimensional mesh with 2x2 elements ...
    tParameters.set_number_of_elements_per_dimension( tElements, 1 );


    // we create a square of 2x2
    tParameters.set_domain_dimensions( 2.0, 1.0 );
    //tParameters.set_domain_offset( -3.0, -3.0 );

    // set default mesh order to two
    tParameters.set_mesh_order( 1 );
    //tParameters.set_mesh_orders_simple( 1 );
    // make mesh output silent
    tParameters.set_verbose( false );

    // B-Spline truncation is turned on by default. It is recommended to leave this
    // setting as is.
    tParameters.set_bspline_truncation( true );
    //tParameters.set_buffer_size( 2 );

//------------------------------------------------------------------------------
//  HMR Initialization
//------------------------------------------------------------------------------

    // We create the HMR object by passing the settings to the constructor
    hmr::HMR tHMR( tParameters );

    // we create one refinement pattern
    tHMR.set_activation_pattern( 0 );

    // the following lines will be replaced by the refinement manager

    for( uint tLevel = 0; tLevel < 1; ++tLevel )
    {
        luint tNumberOfElements = tHMR.get_number_of_elements_on_proc();

        for( uint k=0; k<tNumberOfElements*0.5; ++k )
        {
            tHMR.flag_element( k );
        }

        tHMR.perform_refinement();
    }
    // refine for three levels
    //for( uint tLevel = 0; tLevel < 1; ++tLevel )
    //{
        // get number of elements on mesh

            //luint tNumberOfElements = tHMR.get_number_of_elements_on_proc();

      /*      for( uint k=0; k<8; ++k )
            {
                tHMR.flag_element( k );
            }

            // perform refinement
            tHMR.perform_refinement();

            //tNumberOfElements = tHMR.get_number_of_elements_on_proc();

            for( uint k=0; k<16; ++k )
            {
                tHMR.flag_element( k );
            }

            // perform refinement
            tHMR.perform_refinement(); */

        //for ( uint k=0; k<8; ++k )
        //{
        //    tHMR.flag_element( k );
//
  //      }

    //    tHMR.perform_refinement();
    //}

    // we create one refinement pattern
       tHMR.set_activation_pattern( 1  );

       // the following lines will be replaced by the refinement manager

       // refine for three levels
       for( uint tLevel = 0; tLevel < 1; ++tLevel )
       {


           luint tNumberOfElements = tHMR.get_number_of_elements_on_proc();
           for( uint k=tNumberOfElements*0.5; k<tNumberOfElements; ++k )
           {
               tHMR.flag_element( k );
               //tHMR.flag_element( k+8 );
           }
           // tHMR.flag_element( 0 );

           // perform refinement
           tHMR.perform_refinement();


       }

//------------------------------------------------------------------------------
//  Fields
//------------------------------------------------------------------------------


    // Create fields for all three patterns
    auto tField0 = tHMR.create_field( "LevelSet_Union", 0 );

    tField0->evaluate_function( lvlset );


    auto tField2 = tHMR.create_field( "LevelSet", 2 );



    tHMR.unite_patterns( 0, 1, 2 );
    tHMR.interpolate_field( tField0, tField2 );
    auto tField1 = tHMR.create_field( "LevelSet_Direct", 1 );
    tHMR.extract_field( tField2, tField1 );
    //tField1->evaluate_function( lvlset );
    // einmal rom
    real tError1 = 0;
    tField1->l2_project_coefficients(  tError1, lvlset  );

    // ond wieder nom
    tField1->evaluate_node_values();

   // auto tExact = tHMR.create_field( "Exact", 1 );
   // tExact->evaluate_function( lvlset );


  //  auto tExact2 = tHMR.create_field( "Exact", 2 );
  //  tExact2->evaluate_function( lvlset );


    //Mat< real > tData1 = tField1->get_data();
    //Mat< real > tRef = tExact->get_data() ;
    real tError3 = 0;
    auto tField3 = tHMR.map_field_to_output_mesh( tField0 , tError3, lvlset );

    //real tError1 = 0;
    //tField0->l2_project_coefficients( tError1, lvlset  );
    //tField0->evaluate_node_values();

    std::cout << tField1->get_label() << " " << tError1 << std::endl;
    std::cout << tField3->get_label() << " " << tError3 << std::endl;
    std::cout << tElements << ", " << tError1 <<  " " << tError3 << std::endl;
//------------------------------------------------------------------------------
//    Output
//------------------------------------------------------------------------------
    //auto tField3 = tHMR.create_field( "Project", 3 );
   //tField3->evaluate_node_values( tField1->get_coefficients() );

    // save mesh to file
    tHMR.save_to_exodus( 0, "Mesh0.exo" );
    tHMR.save_to_exodus( 1, "Mesh1.exo" );
    tHMR.save_to_exodus( 2, "Mesh2.exo" );
    //tHMR.save_to_exodus( 3, "Mesh3.exo" ); */
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
