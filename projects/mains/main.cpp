/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */
#include<iostream>
#include<vector>
#include<string>
#include<fn_print.hpp>

#include <cstdio>		// nicer than streams in some respects
// C system files
#include <unistd.h>
// C++ system files
#include <stdio.h>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
//// TPL header files
//#ifdef PARALLEL
//#include "mpi.h"
//#endif
// MORIS header files.
#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
// other header files
//#include <catch.hpp>
//#include "fn_equal_to.hpp" //ALG
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Stopwatch.hpp" //CHR/src
#include "op_move.hpp"

#include "cl_Matrix.hpp"

//#include "cl_Profiler.hpp" //profiler header
//#include <Eigen/Dense>
//#include <armadillo>
#include "cl_Logger.hpp" // MRS/IOS/src

#include "fn_Exec_load_user_library.hpp"

#include "cl_HMR_Paramfile.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

//#include "gperftools/profiler.h"

using namespace moris;
//---------------------------------------------------------------

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 2 );
//    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
//  The main executable is just for developing and testing.
//  Please do not push this file to git.
//------------------------------------------------------------------------------

    //dynamically linked file
    Library_IO tLibrary( argv[ 1 ] );

    moris::ParameterList tParameters = hmr::create_hmr_parameter_list();

    // load user defined function
    MORIS_PARAMETER_FUNCTION tHMRParameters = tLibrary.load_parameter_file( "HMR_Mesh_Paramters" );

    tHMRParameters( tParameters );

    hmr::HMR tHMR( tParameters );

   //initial refinement
    tHMR.perform_initial_refinement( 0 );

    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( 0 );

//    //  create field
//    std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, 0 );
//
//    tField->evaluate_scalar_function( LvlSetCircle_2D );
//
//    for( uint k=0; k<2; ++k )
//    {
//        tHMR.flag_surface_elements_on_working_pattern( tField );
//        tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//        tField->evaluate_scalar_function( LvlSetCircle_2D );
//    }

    tHMR.finalize();

    tHMR.save_to_exodus( 0, "./main_exo/main_test.e" );

    // finalize MORIS global communication manager
    gMorisComm.finalize();
    return 0;
}
