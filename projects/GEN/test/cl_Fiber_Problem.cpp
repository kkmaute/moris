/*
 * cl_Fiber_Problem.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */


#include "catch.hpp"
#include "HDF5_Tools.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "HMR_Globals.hpp"

#include "fn_cylinder_with_end_caps.hpp"

using namespace moris;
namespace ge
{

TEST_CASE("fiber_problem_test_01","[GE],[fiber_test]")
{
///*
// * ------------------------------------------------------------------------------
// *                      HMR Mesh Using Parameter Options
// * ------------------------------------------------------------------------------
// */
//        uint tLagrangeMeshIndex = 0;
//        //  HMR Parameters setup
//        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//
//        tParameters.set( "number_of_elements_per_dimension", "100, 100, 2" );
//        tParameters.set( "domain_dimensions",                "160, 80, 8" );
//        tParameters.set( "domain_offset",                    "-80, -40, -4" );
////        tParameters.set( "number_of_elements_per_dimension", "100, 100, 2" );
////        tParameters.set( "domain_dimensions",                "80, 20, 1" );
////        tParameters.set( "domain_offset",                    "-40, 0, -0.5" );
//
//        tParameters.set( "truncate_bsplines", 1 );
//        tParameters.set( "lagrange_orders", "1" );
//        tParameters.set( "lagrange_pattern", "0" );
//        tParameters.set( "bspline_orders", "1" );
//        tParameters.set( "bspline_pattern", "0" );
//
//        tParameters.set( "lagrange_output_meshes", "0" );
//
//
//        tParameters.set( "lagrange_to_bspline", "0" );
//
//        tParameters.set( "use_multigrid", 0 );
//
//        tParameters.set( "refinement_buffer", 1 );
//        tParameters.set( "staircase_buffer", 1 );
//
//        //  HMR Initialization
//        moris::hmr::HMR tHMR( tParameters );
//
//        // std::shared_ptr< Database >
//        auto tDatabase = tHMR.get_database();
//
//        tHMR.perform_initial_refinement( 0 );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );   // HMR Lagrange mesh
//
////        tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_2_initial.vtk");
////------------------------------------------------------------------------------
//
////------------------------------------------------------------------------------
//moris::tic tTimer0;
//        std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "fibers", tLagrangeMeshIndex);
//moris::tic tTimer1;
//        tField->evaluate_scalar_function( moris::ge::cylinderWithEndCaps );
//real tElapsedTime1 = tTimer1.toc<moris::chronos::milliseconds>().wall;
//tElapsedTime1 /= 1000;
//std::cout<<"==============================================="<< std::endl;
//std::cout<<"Total to evaluate scalar function once: "<< tElapsedTime1 << std::endl << std::endl;
//std::cout<<"==============================================="<< std::endl;
//        for( uint k=0; k<2; ++k )
//        {
//            tHMR.flag_surface_elements_on_working_pattern( tField );
//            tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//            tField->evaluate_scalar_function( moris::ge::cylinderWithEndCaps );
//        }
//        tDatabase->get_background_mesh()->save_to_vtk("Backgroundmesh_initial.vtk");
////------------------------------------------------------------------------------
//
//    tDatabase->update_bspline_meshes();
//    tDatabase->update_lagrange_meshes();
//
//    // calculate T-Matrices etc
//    tDatabase->finalize();
////------------------------------------------------------------------------------
//    tField->evaluate_scalar_function( moris::ge::cylinderWithEndCaps );
//    tHMR.save_to_exodus( 0, "all_fibers.g" );
//
//real tElapsedTime0 = tTimer0.toc<moris::chronos::milliseconds>().wall;
//tElapsedTime0 /= 1000;
//std::cout<<"==============================================="<< std::endl;
//std::cout<<"Total time for evaluation: "<< tElapsedTime0 << std::endl << std::endl;
//std::cout<<"==============================================="<< std::endl;
//
}

}   // ge namespace
