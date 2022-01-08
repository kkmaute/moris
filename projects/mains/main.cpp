#include<iostream>
#include<ios>
#include<limits>

#include "../MRS/IOS/src/cl_Git_info.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Logger.hpp"                // MRS/IOS/src
#include "banner.hpp"                   // COR/src
#include <Kokkos_Core.hpp>

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;
//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

void moris_pause( int  & argc, char * argv[] )
{
    // go through user arguments and look for flags
    for( int k=0; k<argc; ++k )
    {
        // user requests delayed start
        if ( std::string( argv[ k ] ) == "--pause" || std::string( argv[ k ] ) == "-p" )
        {
            if ( par_rank() == 0 )
            {
                std::string dummy;

                std::cout << "Press enter to continue . . .\n" << std::flush;
                std::getline(std::cin, dummy);
            }
        }
    }

    // wait until all processors are ready to continue
    barrier();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // delay start of execution for attaching debugger to process
    moris_pause( argc, argv );

    // Kokkos::initialize(argc, argv);
    barrier();
    int one=1;

    int sumup=sum_all(one);
    if (par_rank() == 0)
    {
        std::cout << " in main 0 " << sumup << std::endl << std::flush ;
    }
    // set severity level 0 - all outputs
    gLogger.initialize( argc, argv );

    barrier();
    sumup=sum_all(one);
    if (par_rank() == 0)
    {
        std::cout << " in main 1 " << sumup << std::endl << std::flush ;
    }
    // print banner
    moris::print_banner( argc, argv );

    barrier();
    sumup=sum_all(one);
    if (par_rank() == 0)
    {
        std::cout << " in main 2 " << sumup << std::endl << std::flush ;
    }

    // print git branch and hash
    if ( par_rank() == 0 )
    {
        git_info tGitInfo;

        std::fprintf( stdout, "\n     GIT branch    : %s\n",  tGitInfo.get_git_branch().c_str() );
        std::fprintf( stdout,   "     GIT revision  : %s\n\n",tGitInfo.get_git_hash().c_str() );

        MORIS_LOG_SPEC("Par Rank",par_rank());
        MORIS_LOG_SPEC("Par Size",par_size());
    }

    barrier();
    sumup=sum_all(one);
    if (par_rank() == 0)
    {
        std::cout << " in main 3 " << sumup << std::endl << std::flush ;
    }


    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    //Kokkos::finalize_all();

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
