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

#include "cl_HMR.hpp"
#include "cl_FEM_Element.hpp" // FEM/INT/src
#include "cl_Model_Solver_Interface.hpp"
#include "cl_Equation_Object.hpp"
#include "cl_FEM_IWG_L2_Test.hpp"
#include "cl_MDL_Model.hpp"


moris::Comm_Manager gMorisComm;

using namespace moris;

real
distance( const Mat< real > & aPoint )
{
    real tValue = 0;
    uint tNumberOfDimensions = aPoint.length();

    for( uint k=0; k< tNumberOfDimensions; ++k )
    {
        tValue += std::pow( aPoint( k ), 2 );
    }
    return std::sqrt( tValue );
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
    // In this part, an HMR object is created

    // create settings object
    /* ( this is actually my internal settings object. As far as I understood,
                 the user is supposed to control HMR using a parameter list.
                 This is not implemented yet ) */

    hmr::Parameters tParameters;

    // create a Mat for a 2D object

    Mat< luint > tNumberOfElements = { { 2 }, { 2 } };

    // mesh orders
    Mat< uint >  tLagrangeOrders = { { 1 }, { 1 }, { 1 } };

    // set interpolation degrees
    tParameters.set_lagrange_orders( tLagrangeOrders );

    Mat< uint >  tBSplineOrders = { { 1 }, { 1 } };
    tParameters.set_bspline_orders( tBSplineOrders );

    // set patterns
    //Parameters.set_default_patterns();

    Mat< uint >  tLagrangePatterns = { { 0 }, { 1 }, { 2 } };
    tParameters.set_lagrange_patterns( tLagrangePatterns );

    // @fixme make sure that B-Spline pattern is always coarser
    Mat< uint >  tBSplinePatterns = { { 0 }, { 1 } };
    tParameters.set_bspline_patterns( tBSplinePatterns );

    Mat< uint >  tLinks = { { 0 },  { 1 },  { 1 } };
    tParameters.set_lagrange_to_bspline( tLinks );


    // pass number of elements to settings
    tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

    // make mesh output silent

    tParameters.set_verbose( true );

    // buffer size must be set at least to max polynomial if truncation is used
    //tParameters.set_buffer_size( tParameters.get_max_polynomial() );

    tParameters.set_buffer_size( 0 );
    tParameters.set_bspline_truncation( false );

    hmr::HMR tHMR( &tParameters );

    tHMR.set_active_pattern( 0 );
    tHMR.flag_element( 0 );
    tHMR.perform_refinement();

    tHMR.set_active_pattern( 1 );
    tHMR.flag_element( 3 );
    tHMR.perform_refinement();


    tHMR.unite_patterns( 0, 1, 2 );


//------------------------------------------------------------------------------

    tHMR.set_active_pattern( 2 );
    tHMR.update_meshes();
    //tHMR.activate_all_t_matrices();

    auto tField0 = tHMR.create_field( "Field0", 0 );

    tField0->evaluate_function( distance );

    auto tField1 = tHMR.create_field( "Field1", 1 );

    // the union
    auto tUnion = tHMR.create_field( "Union",  2 );

    // interpolate field 1 onto field 2
    tHMR.interpolate_field( tField0, tUnion );

    // create a mesh interface
    auto tMesh = tHMR.create_interface();

    // create IWG object
    moris::fem::IWG_L2_Test tIWG;

    // create Model
    mdl::Model tModel( tMesh, tIWG, tUnion->get_data(), tField1->get_coefficients() );

    // tHMR.save_to_exodus( 2, "Mesh.exo" );

    // auto tMesh = tHMR.get_lagrange_mesh_by_index( 2 );

    //tHMR.flag_element( 3 );
    //tHMR.perform_refinement();

    // get number of vertices
    //
    //tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSplines.vtk");

    //tHMR.save_to_hdf5("Mesh.h5");

    //hmr::HMR tHMR2("Mesh.h5");
    //tHMR2.save_to_exodus( 1, "Mesh.exo" );

//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
