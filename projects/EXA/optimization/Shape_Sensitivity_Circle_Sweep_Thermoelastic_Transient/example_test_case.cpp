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
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "paths.hpp"

#include "cl_Logger.hpp"    // MRS/IOS/src
#include "HDF5_Tools.hpp"
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src

using namespace moris;

extern "C" {
bool gThermoelasticUsePetsc = false;
bool gPrintReferenceValues  = false;
}

// test case index
uint gTestCaseIndex;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        const std::string& aExoFileName,
        uint               aTestCaseIndex )
{

    // PETSc uses a different vector-file representation. Its isolated case
    // checks the same physical objective/constraint derivatives against all
    // three FD schemes below; Epetra retains the stored full-adjoint reference.
    if ( !gThermoelasticUsePetsc )
    {
        Matrix< DDRMat > tAdjointRefValues;
        Matrix< DDRMat > tAdjointValues;

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tFieldRefPath = tPrefix + "/projects/EXA/optimization/Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient/Shape_Sensitivity_Thermoelastic_Transient_Ref.hdf5";
        std::string tLabel        = "LHS/Values";

        hid_t  tFileRef = open_hdf5_file( tFieldRefPath );
        herr_t tStatus  = 0;
        load_matrix_from_hdf5_file(
                tFileRef,
                tLabel,
                tAdjointRefValues,
                tStatus );

        tStatus = close_hdf5_file( tFileRef );

        std::string tFieldPath = "";

        if ( aTestCaseIndex == 0 )
        {
            tFieldPath = "./Shape_Sensitivity_Thermoelastic_Transient.hdf5";
        }
        else if ( aTestCaseIndex == 1 )
        {
            tFieldPath = "./Shape_Sensitivity_Thermoelastic_Transient_Staggered.hdf5";
        }

        hid_t tFile = open_hdf5_file( tFieldPath );
        tStatus     = 0;
        load_matrix_from_hdf5_file(
                tFile,
                tLabel,
                tAdjointValues,
                tStatus );

        tStatus = close_hdf5_file( tFile );

        MORIS_ERROR( tStatus == 0, "Field_Example: Status returned != 0, Error in reading values" );

        REQUIRE( tAdjointRefValues.n_rows() == tAdjointValues.n_rows() );
        REQUIRE( tAdjointRefValues.n_cols() == tAdjointValues.n_cols() );

        for ( uint Ik = 0; Ik < tAdjointValues.n_cols(); Ik++ )
        {
            for ( uint Ii = 0; Ii < tAdjointValues.n_rows(); Ii++ )
            {
                CHECK( std::abs( tAdjointValues( Ii, Ik ) - tAdjointRefValues( Ii, Ik ) ) < 1e-10 );
            }
        }
    }

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( "shape_opt_test.hdf5" );
    herr_t tStatus = 0;

    // Declare sensitivity matrices for comparison
    Matrix< DDRMat > tObjectiveAnalytical;
    Matrix< DDRMat > tConstraintsAnalytical;
    Matrix< DDRMat > tObjectiveFD;
    Matrix< DDRMat > tConstraintsFD;

    // Read analytical sensitivities
    load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 analytical", tObjectiveAnalytical, tStatus );
    load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus );
    REQUIRE( tObjectiveAnalytical.length() == tConstraintsAnalytical.length() );    // one objective and one constraint for this problem only

    // Shared physical and sensitivity references for both backends/solver
    // variants. Recorded from the one-rank, two-step Epetra staggered case
    // (2026-11-23), after its stored full-adjoint and FD checks passed. Geometry
    // partials use finite differences, so allow their solver-dependent roundoff
    // at 1e-5 relative accuracy, tighter than the 1e-4 FD comparisons below.
    const real tObjectiveGradientRef[]  = { 0.29689091680120339, 0.26082742822558808, -6.834070431977584 };
    const real tConstraintGradientRef[] = { 4.5183212796229277e-09, -0.014213082765570158, -8.8166036161426025 };
    REQUIRE( tObjectiveAnalytical.numel() == 3 );
    for ( uint i = 0; i < 3; ++i )
    {
        CHECK( tObjectiveAnalytical( i ) == Approx( tObjectiveGradientRef[ i ] ).epsilon( 1e-5 ).margin( 1e-8 ) );
        CHECK( tConstraintsAnalytical( i ) == Approx( tConstraintGradientRef[ i ] ).epsilon( 1e-5 ).margin( 1e-8 ) );
        if ( gPrintReferenceValues )
        {
            MORIS_LOG_INFO( "Gradient reference ADV(%u): objective %.17g, constraint %.17g",
                    i,
                    tObjectiveAnalytical( i ),
                    tConstraintsAnalytical( i ) );
        }
    }

    Matrix< DDRMat > tObjective, tConstraint;
    load_matrix_from_hdf5_file( tFileID, "objectives eval_1-1", tObjective, tStatus );
    REQUIRE( tStatus == 0 );
    load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraint, tStatus );
    REQUIRE( tStatus == 0 );
    REQUIRE( tObjective.numel() == 1 );
    REQUIRE( tConstraint.numel() == 1 );
    CHECK( tObjective( 0 ) == Approx( 17.813255834439584 ).epsilon( 1e-9 ) );
    CHECK( tConstraint( 0 ) == Approx( 17.905548482069083 ).epsilon( 1e-9 ) );
    if ( gPrintReferenceValues )
    {
        MORIS_LOG_INFO( "Physical references: objective %.17g, constraint %.17g", tObjective( 0 ), tConstraint( 0 ) );
    }

    // perturbation of denominator when building relative error
    real tDeltaEps = 1.0e-14;

    // Read FD sensitivities and compare
    Vector< std::string > tFDTypes = { "fd_forward", "fd_backward", "fd_central" };
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
                    100.0 * tRelObectiveDifference );

            MORIS_LOG_INFO( "Check derivative of constraint wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                    tADVIndex,
                    tConstraintsAnalytical( tADVIndex ),
                    tFDTypes( tFDIndex ).c_str(),
                    tConstraintsFD( tADVIndex ),
                    100.0 * tRelConstraintDifference );

            if ( std::abs( tObjectiveFD( tADVIndex ) ) > 1e-6 )
            {
                CHECK( tRelObectiveDifference < 1e-4 );
            }
            else
            {
                CHECK( tObjectiveAnalytical( tADVIndex ) == Approx( tObjectiveFD( tADVIndex ) ).margin( 1e-6 ) );
            }

            if ( std::abs( tConstraintsFD( tADVIndex ) ) > 1e-6 )
            {
                CHECK( tRelConstraintDifference < 1e-4 );
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

TEST_CASE( "Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient",
        "[moris],[example],[optimization],[sweep],[sweep_thermoelastic_transient]" )
{
    // Two-step thermoelastic circle sweep: compare the final adjoint with the
    // stored reference and objective/constraint gradients with three FD schemes.
    // The temporal carry must remain fixed through each monolithic adjoint solve.
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // set test case index
    gTestCaseIndex = 0;

    // perform check for Test Case 0
    check_results( "ShapeSensitivitiesThermoelasticTransient.exo", gTestCaseIndex );
}

TEST_CASE( "Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient_Staggered",
        "[moris],[example],[optimization],[sweep],[sweep_thermoelastic_transient_staggered]" )
{
    // Solve the same transient model by thermal/displacement subsystems. The
    // all-field temporal carry is updated only after the complete adjoint solve;
    // use the same adjoint reference and FD gradient tolerances as the monolithic case.
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient_Staggered.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // set test case index
    gTestCaseIndex = 1;

    // perform check for Test Case 0
    check_results( "ShapeSensitivitiesThermoelasticTransientStaggered.exo", gTestCaseIndex );
}

#ifdef MORIS_HAVE_PETSC
// The same two-step thermal/displacement model and FD tolerances exercise the
// all-state history carry entering both PETSc subsystem maps. Run separately:
// PETSc cannot be initialized again after a workflow finalizes it.
TEST_CASE( "PETSc staggered thermoelastic transient sensitivity",
        "[.][ThermoelasticPetscAdjoint]" )
{
    gThermoelasticUsePetsc = true;
    gTestCaseIndex         = 1;
    REQUIRE( moris::par_size() == 1 );

    // Remove every checked output so a failed workflow cannot pass using files
    // left by another backend or an earlier run in this dedicated directory.
    std::remove( "shape_opt_test.hdf5" );
    std::remove( "ShapeSensitivitiesThermoelasticTransientStaggered.exo" );
    std::string         tLibrary = MORIS_THERMOELASTIC_STAGGERED_LIBRARY;
    std::vector< char > tArgument( tLibrary.begin(), tLibrary.end() );
    tArgument.push_back( '\0' );
    char  tName[]      = "";
    char* tArguments[] = { tName, tArgument.data() };
    REQUIRE( fn_WRK_Workflow_Main_Interface( 2, tArguments ) == 0 );
    check_results( "ShapeSensitivitiesThermoelasticTransientStaggered.exo", gTestCaseIndex );
    gThermoelasticUsePetsc = false;
}
#endif
