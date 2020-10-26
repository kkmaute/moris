//
// example specific interface to moris
//

#include <catch.hpp>

#include "cl_Logger.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE("SIMP",
        "[moris],[example],[optimization],[levelset_boxbeam]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );

    // Read Exodus file
    moris::mtk::Exodus_IO_Helper tExoIO("levelset_boxbeam.exo", 0, false, false);

    // Checks
}
