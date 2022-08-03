

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_IO.hpp"
#include "cl_Communication_Tools.hpp"

using namespace moris;


TEST_CASE( "WRK_sensitivity_test ", "[moris],[WRK_sensitivity_test]" )
{
    if( par_size() == 1 )
    {
        std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_1.so";

        std::shared_ptr< Library_IO > tLibrary = std::make_shared< Library_IO >( tInputFilePath );

        wrk::Performer_Manager tPerformerManager( tLibrary );
        wrk::Workflow_HMR_XTK tWorkflow( &tPerformerManager );

        Matrix<DDRMat> tDummy(1, 1, 0.0);
        Matrix<IdMat> tDummy1(1, 1, 0.0);
        tWorkflow.initialize(tDummy, tDummy, tDummy,tDummy1);
        tWorkflow.perform(tDummy);

//        tWorkflow.get_criteria(tADVs);
    }
}
