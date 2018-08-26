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

    // We create a 2-Dimensional mesh with 10x10x10 elements ...
    tParameters.set_number_of_elements_per_dimension( 2, 2, 2);

    // we create a square of 1x1
    tParameters.set_domain_dimensions( 2, 2, 2 );

    // set default mesh order to two
    tParameters.set_mesh_order( 2 );
    //tParameters.set_mesh_orders_simple( 1 );
    // make mesh output silent
    tParameters.set_verbose( true );

    // B-Spline truncation is turned on by default. It is recommended to leave this
    // setting as is.
    tParameters.set_bspline_truncation( true );

//------------------------------------------------------------------------------
//  HMR Initialization
//------------------------------------------------------------------------------

    // We create the HMR object by passing the settings to the constructor
    hmr::HMR tHMR( tParameters );

    // we create one refinement pattern
    tHMR.set_active_pattern( 0 );

    // the following lines will be replaced by the refinement manager

    // refine for three levels
    for( uint tLevel = 0; tLevel < 5; ++tLevel )
    {
        // get number of elements on mesh
        luint tNumberOfElements = tHMR.get_number_of_elements_on_proc();

        uint k = 2;
        // some elements
        for( luint e=0; e<tNumberOfElements; e+=k )
        {
            tHMR.flag_element( e );
            k = k + 1;
            if( k == 13 )
            {
                k = 2;
            }
        }


        // perform refinement
        tHMR.perform_refinement();
    }

    //tHMR.update_meshes();

    //tHMR.mBackgroundMesh->save_to_vtk("BackgroundMesh.vtk");


    //tHMR.mBackgroundMesh->print_active_elements();



    tHMR.set_active_pattern( 1 );

    // refine for three levels
    for( uint tLevel = 0; tLevel < 5; ++tLevel )
    {
        // get number of elements on mesh
        luint tNumberOfElements = tHMR.get_number_of_elements_on_proc();
        uint k = 1;
        // some elements
        for( luint e=0; e<tNumberOfElements; e+=k )
        {
            tHMR.flag_element( e );
            k = k + 1;
            if( k == 7 )
            {
                k = 1;
            }
        }

        // perform refinement
        tHMR.perform_refinement();
    }

    // Unite both patterns and store them to pattern 2
    tHMR.unite_patterns( 0, 1, 2 );

   // exit( 0 );
//------------------------------------------------------------------------------
//  Fields
//------------------------------------------------------------------------------

    tHMR.set_active_pattern( 0 ); // 2
    tHMR.activate_all_t_matrices(); // < -- this is a problem and needs fixing

    //tHMR.update_meshes();

    // Create fields for all three patterns
    auto tField0 = tHMR.create_field( "Field", 0 );
    auto tField1 = tHMR.create_field( "Field", 1 );
    auto tField2 = tHMR.create_field( "Field", 2 );

    auto tExact = tHMR.create_field( "Exact", 1 );

    tField0->evaluate_function( distance );
    //tField0->l2_project_coefficients(); //<-- this still causes an error

    tExact->evaluate_function( distance );

    // the union

    tHMR.set_active_pattern( 2 ); // 2
    tHMR.activate_all_t_matrices(); // < -- this is a problem and needs fixing

    // map data from field 1 to field 2
    // interpolate field 1 onto field 2
    tHMR.interpolate_field( tField0, tField2 );
    //tField2->evaluate_function( distance );

    // create mesh interface
    auto tMesh = tHMR.create_interface( 2 );

    // create IWG object
    moris::fem::IWG_L2 tIWG;


    mdl::Model tModel( tMesh, tIWG, tField2->get_data(), tField1->get_coefficients() );



    tField1->evaluate_node_values();

//------------------------------------------------------------------------------
//   Test error
//------------------------------------------------------------------------------


    real tError = r2( tExact->get_data(), tField1->get_data() );

    std::cout << "R2 : " << tError << std::endl;

//------------------------------------------------------------------------------
//    Output
//------------------------------------------------------------------------------

    // save mesh to file
    tHMR.save_to_exodus( 0, "Mesh0.exo" );
    tHMR.save_to_exodus( 1, "Mesh1.exo" );
    tHMR.save_to_exodus( 2, "Mesh2.exo" );
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
