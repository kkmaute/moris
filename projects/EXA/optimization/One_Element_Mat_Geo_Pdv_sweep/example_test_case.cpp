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

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

// global variable for interpolation order
uint gInterpolationOrder;

// problem dimension: 2D or 3D
uint gDim;

// test case index
uint gTestCaseIndex;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        std::string aExoFileName,
        std::string aHdf5FileName,
        uint        aTestCaseIndex )
{
}
TEST_CASE( "One_Element_Mat_Geo_Pdv_sweep",
        "[moris],[example],[optimization],[One_Element_Mat_Geo_Pdv_sweep]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "One_Element_Mat_Geo_Pdv_sweep.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( "shape_opt_test.hdf5" );
    herr_t tStatus = 0;

    if ( par_size() == 1 )
    {
        // perturbation of denominator when building relative error
        real tDeltaEps = 1.0e-14;

        // Declare sensitivity matrices for comparison
        Matrix< DDRMat > tObjectiveAnalytical;
        Matrix< DDRMat > tConstraintsAnalytical;
        Matrix< DDRMat > tObjectiveFD;
        Matrix< DDRMat > tConstraintsFD;

        // Read analytical sensitivities
        load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 analytical", tObjectiveAnalytical, tStatus );
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus );
        REQUIRE( tObjectiveAnalytical.length() == tConstraintsAnalytical.length() );    // one objective and one constraint for this problem only

        // Read FD sensitivities and compare
        Cell< std::string > tFDTypes = { "fd_forward", "fd_backward", "fd_central" };
        for ( uint tFDIndex = 0; tFDIndex < tFDTypes.size(); tFDIndex++ )
        {
            load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tObjectiveFD, tStatus );
            load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tConstraintsFD, tStatus );

            REQUIRE( tObjectiveAnalytical.length() == tObjectiveFD.length() );
            REQUIRE( tConstraintsAnalytical.length() == tConstraintsFD.length() );

            for ( uint tADVIndex = 0; tADVIndex < tObjectiveAnalytical.length(); tADVIndex++ )
            {
                real tRelObectiveDifference =
                        std::abs( tObjectiveAnalytical( tADVIndex ) - tObjectiveFD( tADVIndex ) ) /    //
                        ( std::abs( tObjectiveFD( tADVIndex ) ) + tDeltaEps );

                real tRelConstraintDifference =
                        std::abs( tConstraintsAnalytical( tADVIndex ) - tConstraintsFD( tADVIndex ) ) /    //
                        ( std::abs( tConstraintsFD( tADVIndex ) ) + tDeltaEps );

                MORIS_LOG_INFO( "Check derivative of objective  wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                        tADVIndex,
                        tObjectiveAnalytical( tADVIndex ),
                        tFDTypes( tFDIndex ).c_str(),
                        tObjectiveFD( tADVIndex ),
                        100 * tRelObectiveDifference );

                MORIS_LOG_INFO( "Check derivative of constraint wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                        tADVIndex,
                        tConstraintsAnalytical( tADVIndex ),
                        tFDTypes( tFDIndex ).c_str(),
                        tConstraintsFD( tADVIndex ),
                        100 * tRelConstraintDifference );

                if ( std::abs( tObjectiveFD( tADVIndex ) ) > 1e-6 )
                {
                    CHECK( tRelObectiveDifference < 1e-5 );
                }
                else
                {
                    CHECK( tObjectiveAnalytical( tADVIndex ) == Approx( tObjectiveFD( tADVIndex ) ).margin( 1e-6 ) );
                }

                if ( std::abs( tConstraintsFD( tADVIndex ) ) > 1e-6 )
                {
                    CHECK( tRelConstraintDifference < 1e-5 );
                }
                else
                {
                    CHECK( tConstraintsAnalytical( tADVIndex ) == Approx( tConstraintsFD( tADVIndex ) ).margin( 1e-6 ) );
                }
            }
        }

        // close file
        close_hdf5_file( tFileID );
    }
}
