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
distance( const Mat< real > & aPoint )
{
    return norm ( aPoint );
    //return std::cos( aPoint( 0 ) * 3.141592653589793 ) * std::pow(  aPoint( 1 ), 2 );
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

    // We create a 2-Dimensional mesh with 2x2 elements ...
    tParameters.set_number_of_elements_per_dimension( 4, 4 );

    // we create a square of 1x1
    // tParameters.set_domain_dimensions( 2, 2 );

    // set default mesh order to two
    tParameters.set_mesh_order( 2 );
    //tParameters.set_mesh_orders_simple( 1 );
    // make mesh output silent
    tParameters.set_verbose( false );

    // B-Spline truncation is turned on by default. It is recommended to leave this
    // setting as is.
    tParameters.set_bspline_truncation( false );

//------------------------------------------------------------------------------
//  HMR Initialization
//------------------------------------------------------------------------------

    // We create the HMR object by passing the settings to the constructor
    hmr::HMR tHMR( tParameters );

    // we create one refinement pattern
    tHMR.set_activation_pattern( 0 );

    // the following lines will be replaced by the refinement manager

    // refine for three levels
    for( uint tLevel = 0; tLevel < 4; ++tLevel )
    {
        // get number of elements on mesh
        //luint tNumberOfElements = tHMR.get_number_of_elements_on_proc();

        tHMR.flag_element( 0 );

        // perform refinement
        tHMR.perform_refinement();
    }


    tHMR.set_activation_pattern( 1 );

    // refine for three levels
    for( uint tLevel = 0; tLevel < 4; ++tLevel )
    {
        tHMR.flag_element(  tHMR.get_number_of_elements_on_proc()-1 );
        // perform refinement
        tHMR.perform_refinement();
    }

    //tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");
    //tHMR.mLagrangeMeshes( 1 )->save_to_vtk("Lagrange.vtk");
    // Unite both patterns and store them to pattern 2
    //tHMR.unite_patterns( 0, 1, 2 );
   // exit( 0 );
//------------------------------------------------------------------------------
//  Fields
//------------------------------------------------------------------------------


    // Create fields for all three patterns
       auto tField0 = tHMR.create_field( "Field", 0 );
       tField0->evaluate_function( distance );


    real tError = 0;

    auto tField1 = tHMR.map_field_to_output_mesh( tField0, tError, distance );

    std::cout << par_rank() << " error: " << tError << std::endl;

    //tField0->l2_project_coefficients(); //<-- this still causes an error

    auto tExact = tHMR.create_field( "Exact", 1 );
    tExact->evaluate_function( distance );

//------------------------------------------------------------------------------
//   Test error
//------------------------------------------------------------------------------


    real tR2 = r2( tExact->get_data(), tField1->get_data() );

    std::cout << "R2 : " << tR2 << std::endl;

//------------------------------------------------------------------------------
//    Output
//------------------------------------------------------------------------------

    // save mesh to file
    //tHMR.save_to_exodus( 0, "Mesh0.exo" );
    //tHMR.save_to_exodus( 1, "Mesh1.exo" );
    tHMR.add_field( tField1 );
    tHMR.add_field( tExact );

    tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");
    tHMR.save_to_exodus( "Mesh.exo" );

    delete tField0;
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
