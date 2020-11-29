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

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // Kokkos::initialize(argc, argv);

    // set severity level 0 - all outputs
    gLogger.initialize( 2 );

    // print banner
    moris::print_banner( argc, argv );

    // print git branch and hash
    if ( par_rank() == 0 )
    {
        git_info tGitInfo;
        std::fprintf( stdout, "\n     GIT branch    : %s\n",  tGitInfo.get_git_branch().c_str() );
        std::fprintf( stdout,   "     GIT revision  : %s\n\n",tGitInfo.get_git_hash().c_str() );
    }
    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    //Kokkos::finalize_all();

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
