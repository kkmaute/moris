//
// example specific interface to moris
//

#include <catch.hpp>
#include "paths.hpp"

#include "cl_Logger.hpp"                // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"  // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

//// global variable for interpolation order
//uint gInterpolationOrder;
//
//// problem dimension: 2D or 3D
//uint gDim;
//
//// test case index
//uint gTestCaseIndex;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE("Field_example_write",
        "[moris],[example],[structure],[Plate_Temperature_Field_Write]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Beam_Temperature_Field.so";

    char * argv[2] = {tString1,tString2};

//    // set interpolation order
//    gInterpolationOrder = 1;
//
//    MORIS_LOG_INFO("");
//    MORIS_LOG_INFO("Executing Field_Example_Write - 2D: Interpolation order 1 - %i Processors.",par_size());
//    MORIS_LOG_INFO("");
//
//    // set dimension: 2D
//    gDim = 2;
//
//    // set test case index
//    gTestCaseIndex = 0;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

}
