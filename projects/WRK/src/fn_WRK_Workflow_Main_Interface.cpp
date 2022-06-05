#include<iostream>
#include<vector>
#include<string>
#include<fn_print.hpp>

#include <cstdio>// nicer than streams in some respects
// C system files
#include <unistd.h>
// C++ system files
#include <stdio.h>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
// other header files
//#include <catch.hpp>
//#include "fn_equal_to.hpp" //ALG
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Stopwatch.hpp" //CHR/src
#include "op_move.hpp"

#include "cl_Matrix.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src

#include "cl_WRK_Workflow_Factory.hpp"
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "cl_OPT_Manager.hpp"

#include "cl_Library_IO.hpp"

using namespace moris;

// Parameter function
typedef void ( *Parameter_Function ) ( moris::Cell< moris::Cell< moris::ParameterList > > & aParameterList );

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] )
{    
    if (argc < 2)
    {
        std::cout << "\n Error: input file required\n" << "\n";
        return -1;
    }

    // last input argument is assumed to be shared object file
    // FIXME: shared object file name should be identified with together with other command line options
    std::string tInputArg = std::string(argv[ argc - 1 ]);
    std::string tString = "Reading dynamically linked shared object " + tInputArg + ".";
    MORIS_LOG( tString.c_str() );

    //dynamically linked file
    std::shared_ptr< Library_IO >tLibrary = std::make_shared< Library_IO >( argv[ 1 ] );

    {
        // load the OPT parameter list
        std::string tOPTString = "OPTParameterList";
        Parameter_Function tOPTParameterListFunc = tLibrary->load_function<Parameter_Function>( tOPTString );
        moris::Cell< moris::Cell< ParameterList > > tOPTParameterList;
        tOPTParameterListFunc( tOPTParameterList );

        // Create performer manager
        wrk::Performer_Manager tPerformerManager( tLibrary );

        // FIXME: get this from parameter list"workflow"
        std::string tWRKFlowStr = tOPTParameterList( 0 )( 0 ).get< std::string >("workflow");

        moris::Cell<std::shared_ptr<moris::opt::Criteria_Interface>> tWorkflows = {
                wrk::create_workflow(tWRKFlowStr, &tPerformerManager) };

        if( tOPTParameterList( 0 )( 0 ).get< bool >("is_optimization_problem") )
        {
            moris::opt::Manager tManager( tOPTParameterList, tWorkflows );
            tManager.perform();
        }
        else
        {
            Matrix<DDRMat> tADVs(0, 0);
            Matrix<DDRMat> tDummyBounds;
            tWorkflows(0)->initialize(tADVs, tDummyBounds, tDummyBounds);
            Matrix<DDRMat> tIQIVal = tWorkflows(0)->get_criteria(tADVs);

            print(tIQIVal, "IQI values");
        }
    }

    return 0;
}
