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
#include "paths.hpp"

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

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        const std::string &aExoFileName,
        uint               aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 2, 2 };
    Vector< uint > tReferenceNumNodes = { 6067, 6067, 6067 };
    Vector< uint > tReferenceNumElems = { 1570, 1570, 1570 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( aTestCaseIndex ) );

    if ( aTestCaseIndex == 0 )
    {
        Matrix< DDRMat > tNodalRefFieldValues;
        Matrix< DDRMat > tNodalFieldValues;

        std::string tPrefix       = moris::get_base_moris_dir();
        std::string tFieldRefPath = tPrefix + "/projects/EXA/thermal/diffusion/Field_Example/Field_example_ref.hdf5";
        std::string tLabel        = "FieldNodalTEMP";

        hid_t  tFileRef = open_hdf5_file( tFieldRefPath );
        herr_t tStatus  = 0;
        load_matrix_from_hdf5_file(
                tFileRef,
                tLabel,
                tNodalRefFieldValues,
                tStatus );

        tStatus = close_hdf5_file( tFileRef );

        std::string tFieldPath = "./Field_example_write.hdf5";

        hid_t tFile = open_hdf5_file( tFieldPath );
        tStatus     = 0;
        load_matrix_from_hdf5_file(
                tFile,
                tLabel,
                tNodalFieldValues,
                tStatus );

        tStatus = close_hdf5_file( tFile );

        MORIS_ERROR( tStatus == 0, "Field_Example: Status returned != 0, Error in reading values" );

        CHECK( tNodalRefFieldValues.numel() == tNodalFieldValues.numel() );

        for ( uint Ik = 0; Ik < tNodalFieldValues.numel(); Ik++ )
        {
            CHECK( tNodalFieldValues( Ik ) - tNodalRefFieldValues( Ik ) < 1e-12 );
        }
    }

    if ( aTestCaseIndex == 1 )
    {
        uint tNodeId = 5992;
        // define reference coordinates for node aNodeId
        Matrix< DDRMat > tReferenceCoordinate = { { 0.175 }, { 0.775 } };

        // check nodal coordinates
        real tRelDiffNorm = moris::norm( tExoIO.get_nodal_coordinate( tNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );
        REQUIRE( tRelDiffNorm < 1.0e-12 );

        // check temperature at node aNodeId in first time step (diff is 3rd nodal field, first time step has index 0)
        real tReferenceDiff = 0.0;
        real tRelDifference = std::abs( ( tExoIO.get_nodal_field_value( tNodeId, 3, 0 ) - tReferenceDiff ) );
        REQUIRE( tRelDifference < 1.0e-12 );
    }

    if ( aTestCaseIndex == 2 )
    {
        uint tNodeId = 5992;
        // define reference coordinates for node aNodeId
        Matrix< DDRMat > tReferenceCoordinate = { { 0.175 }, { 0.775 } };

        // check nodal coordinates
        real tRelDiffNorm = moris::norm( tExoIO.get_nodal_coordinate( tNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );
        REQUIRE( tRelDiffNorm < 1.0e-12 );

        // check temperature at node aNodeId in first time step (diff is 3rd nodal field, first time step has index 0)
        real tReferenceDiff = 5.30216e-6;
        real tRelDifference = std::abs( ( tExoIO.get_nodal_field_value( tNodeId, 3, 0 ) - tReferenceDiff ) );
        REQUIRE( tRelDifference < 1.0e-10 );
    }
}

//---------------------------------------------------------------

TEST_CASE( "Field_example_write",
        "[moris],[example],[thermal],[Field_example_write]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Field_example_write.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Field_Example_Write - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 0;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 0
    check_results( "Field_example_write.exo", gTestCaseIndex );
}

//---------------------------------------------------------------

TEST_CASE( "Field_example_read",
        "[moris],[example],[thermal],[Field_example_read]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Field_example_read.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Field_Example_Read - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 1;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 0
    check_results( "Field_example_read.exo", gTestCaseIndex );
}

//---------------------------------------------------------------

TEST_CASE( "Field_example_compare",
        "[moris],[example],[thermal],[Field_example_compare]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Field_example_compare.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Field_Example_Compare - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 2;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 0
    check_results( "Field_example_compare.exo", gTestCaseIndex );
}
