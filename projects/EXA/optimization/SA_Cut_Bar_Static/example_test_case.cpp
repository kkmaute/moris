/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 */

#include <catch.hpp>

#include "cl_Logger.hpp" // MRS/IOS/src
#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE("SA_Cut_Bar_Static",
        "[moris],[example],[optimization],[sweep],[sweep_static]")
{
    // Tolerance for adjoint vs. FD sensitivities
    moris::real tToleranceSensties = 0.01;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "SA_Cut_Bar_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );

    // Sweep HDF5 file
    hid_t tFileID = open_hdf5_file( "SA_Cut_Bar_Static.hdf5" );
    herr_t tStatus = 0;

    // Declare sensitivity matrices for comparison
    Matrix<DDRMat> tObjectiveAnalytical;
    Matrix<DDRMat> tConstraintsAnalytical;
    Matrix<DDRMat> tObjectiveFD;
    Matrix<DDRMat> tConstraintsFD;

    // Read analytical sensitivities
    load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 analytical", tObjectiveAnalytical, tStatus);
    load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus);
    REQUIRE(tObjectiveAnalytical.length() == tConstraintsAnalytical.length()); // one objective and one constraint for this problem only

    // Read FD sensitivities and compare
    Vector<std::string> tFDTypes = {"fd_forward", "fd_backward", "fd_central"};
    for (uint tFDIndex = 0; tFDIndex < tFDTypes.size(); tFDIndex++)
    {
        load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 epsilon_1-1 " + tFDTypes(tFDIndex), tObjectiveFD, tStatus);
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 epsilon_1-1 " + tFDTypes(tFDIndex), tConstraintsFD, tStatus);

        REQUIRE(tObjectiveAnalytical.length() == tObjectiveFD.length());
        REQUIRE(tConstraintsAnalytical.length() == tConstraintsFD.length());

        for (uint tADVIndex = 0; tADVIndex < tObjectiveAnalytical.length(); tADVIndex++)
        {
            MORIS_LOG_INFO("Check derivative of objective  wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                    tADVIndex,
                    tObjectiveAnalytical(tADVIndex),
                    tFDTypes(tFDIndex).c_str(),
                    tObjectiveFD(tADVIndex),
                    100*std::abs((tObjectiveAnalytical(tADVIndex)-tObjectiveFD(tADVIndex))/tObjectiveFD(tADVIndex)));

            MORIS_LOG_INFO("Check derivative of constraint wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                    tADVIndex,
                    tConstraintsAnalytical(tADVIndex),
                    tFDTypes(tFDIndex).c_str(),
                    tConstraintsFD(tADVIndex),
                    100*std::abs((tConstraintsAnalytical(tADVIndex)-tConstraintsFD(tADVIndex))/tConstraintsFD(tADVIndex)));

            CHECK( std::abs( ( tObjectiveAnalytical(   tADVIndex ) - tObjectiveFD(   tADVIndex ) ) /
                    tObjectiveFD(   tADVIndex ) ) < tToleranceSensties );
            CHECK( std::abs( ( tConstraintsAnalytical( tADVIndex ) - tConstraintsFD( tADVIndex ) ) /
                    tConstraintsFD( tADVIndex ) ) < 0.001 );
        }
    }

    // close file
    close_hdf5_file( tFileID );
}

