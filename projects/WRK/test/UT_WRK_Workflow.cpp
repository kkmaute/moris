/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Workflow.cpp
 *
 */

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_Factory.hpp"
#include "cl_Communication_Tools.hpp"

using namespace moris;

TEST_CASE( "WRK_sensitivity_test", "[moris],[WRK_sensitivity_test]" )
{
    if ( par_size() == 1 )
    {
        std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_1.so";

        // Load library
        moris::Library_Factory        tLibraryFactory;
        std::shared_ptr< Library_IO > tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
        tLibrary->load_parameter_list( tInputFilePath, File_Type::SO_FILE );
        tLibrary->finalize();

        wrk::Performer_Manager tPerformerManager( tLibrary );
        wrk::Workflow_HMR_XTK  tWorkflow( &tPerformerManager );

        Matrix< DDRMat > tDummy( 1, 1, 0.0 );
        Matrix< IdMat >  tDummy1( 1, 1, 0.0 );
        tWorkflow.initialize( tDummy, tDummy, tDummy, tDummy1 );
        tWorkflow.perform( tDummy );

        // tWorkflow.get_criteria(tADVs);
    }
}
