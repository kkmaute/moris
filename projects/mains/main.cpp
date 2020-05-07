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
#include "cl_Logger.hpp" // MRS/IOS/src

#include "fn_Exec_load_user_library.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"

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
    moris::print_banner( argc, argv );

    //------------------------------------------------------------------------------
    //  The main executable is just for developing and testing.
    //  Please do not push this file to git.
    //------------------------------------------------------------------------------

    if (argc < 2)
    {
        std::cout << "\n Error: input file required\n" << "\n";
        return -1;
    }

    std::string tInputArg = std::string(argv[ 1 ]);
    std::string tString = "Reading dynamically linked shared object " + tInputArg + ".";
    MORIS_LOG( tString.c_str() );

    //dynamically linked file
    std::shared_ptr< Library_IO >tLibrary = std::make_shared< Library_IO >( argv[ 1 ] );

    {
        wrk::Performer_Manager tPerformerManager( tLibrary );

        tPerformerManager.initialize_performers();

        tPerformerManager.set_performer_cooperations();
        {
            wrk::Workflow tWorkflow( &tPerformerManager );

            Matrix<DDRMat> tADVs(1, 1, 0.0);
            tWorkflow.perform(tADVs);
        }
    }

    // finalize MORIS global communication manager
    gMorisComm.finalize();
    return 0;
}
