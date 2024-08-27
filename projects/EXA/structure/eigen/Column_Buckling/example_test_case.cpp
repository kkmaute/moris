/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp.template
 *
 */

#include <catch.hpp>

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

//---------------------------------------------------------------

// interpolation order
std::string tOrder;

// spatial dimensions
uint tDim;

// stress used in geometric stiffness (PK2 or Linear)
std::string tStressType;

// output file name
std::string tOutputFileName;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        std::string aExoFileName,
        uint        aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( "Checking Results - Reading from %s", aExoFileName.c_str() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        // coordinates of reference point
        real tVolume = tExoIO.get_global_variable( 0, 0 );
        std::cout << "Volume:        " << tVolume << std::endl;
        std::cout << "Eigen value 1: " << tExoIO.get_global_variable( 1, 0 ) / tVolume << std::endl;
        std::cout << "Eigen value 2: " << tExoIO.get_global_variable( 2, 0 ) / tVolume << std::endl;
        std::cout << "Eigen value 3: " << tExoIO.get_global_variable( 3, 0 ) / tVolume << std::endl;
        std::cout << "Eigen value 4: " << tExoIO.get_global_variable( 4, 0 ) / tVolume << std::endl;
        std::cout << "Eigen value 5: " << tExoIO.get_global_variable( 5, 0 ) / tVolume << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    // clang-format off
    // test case                            0    1    2    3     4    5     6    7    8    9     10    11
    Vector< uint > tReferenceNumDims  = {   2,   2,   2,   2,    2,   2,    3,   3,   2,   2,     3,    3 };
    Vector< uint > tReferenceNumNodes = { 306, 104, 306, 104, 1111, 357, 1836, 624, 698, 195, 10184, 2803 };
    Vector< uint > tReferenceNumElems = { 250,  75, 250,  75,  250,  75, 1250, 375, 522, 133, 16168, 4120 };
    // clang-format on

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( aTestCaseIndex ) );

    // define reference volume
    // test case                            1   2   3   4   5   6   7   8   9
    real tReferenceVolume = 10.0;

    real tActualVolume = tExoIO.get_global_variable( 0, 0 );

    real tRelVolumeDifference = std::abs( tActualVolume - tReferenceVolume ) / tReferenceVolume;

    MORIS_LOG_INFO( "Check volume:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceVolume,
            tActualVolume,
            tRelVolumeDifference * 100.0 );

    REQUIRE( tRelVolumeDifference < 1.0e-5 );

    // check lowest five eigen values
    Vector< Matrix< DDRMat > > tReferenceEval;

    // clang-format off
    // test case                  0         1         2         3         4         5         6         7         8         9        10         11
    tReferenceEval.push_back( { { 0.208264, 0.208264, 0.209892, 0.209892, 0.204191, 0.204191, 0.209678, 0.209678, 0.208056, 0.208056, 0.209338, 0.209338 } } );
    tReferenceEval.push_back( { { 1.76392,  1.76392,  1.77737,  1.77737,  1.72616,  1.72616,  0.209678, 0.209678, 1.76131,  1.76131,  0.209338, 0.209338 } } );
    tReferenceEval.push_back( { { 4.39768,  4.39768,  4.42993,  4.42993,  4.28966,  4.28966,  1.77607,  1.77607,  4.38765,  4.38765,  1.77238,  1.77238  } } );
    tReferenceEval.push_back( { { 7.51133,  7.51133,  7.56415,  7.56415,  7.29788,  7.29788,  1.77607,  1.77607,  7.48695,  7.48695,  1.77238,  1.77238  } } );
    tReferenceEval.push_back( { { 10.6635,  10.6635,  10.7357,  10.7357,  10.3161,  10.3161,  4.42856,  4.42856,  10.618,   10.618,   4.41592,  4.41592  } } );
    // clang-format on

    Matrix< DDRMat >
            tActualEigenValues = { {
                    tExoIO.get_global_variable( 1, 0 ) / tActualVolume,
                    tExoIO.get_global_variable( 2, 0 ) / tActualVolume,
                    tExoIO.get_global_variable( 3, 0 ) / tActualVolume,
                    tExoIO.get_global_variable( 4, 0 ) / tActualVolume,
                    tExoIO.get_global_variable( 5, 0 ) / tActualVolume,
            } };

    for ( uint i = 0; i < tReferenceEval.size(); ++i )
    {
        real tRelEvalDifference = std::abs( tActualEigenValues( i ) - tReferenceEval( i )( aTestCaseIndex ) ) / tReferenceEval( i )( aTestCaseIndex );

        MORIS_LOG_INFO( "Check eigen value %i:  reference %12.5e, actual %12.5e, percent error %12.5e.",
                i + 1,
                tReferenceEval( i )( aTestCaseIndex ),
                tActualEigenValues( i ),
                tRelEvalDifference * 100.0 );

        REQUIRE( tRelEvalDifference < 1.0e-5 );
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------
#ifdef MORIS_HAVE_SLEPC

TEST_CASE( "Column_Buckling_2D",
        "[moris],[example],[structure],[eigen]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Column_Buckling.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // test case index
    uint tTestCase = par_size() == 1 ? 0 : 1;

    // interpolation order
    tOrder = "1";

    // spatial dimensions
    tDim = 2;

    // stress used in geometric stiffness (PK2 or Linear)
    tStressType = "Linear";

    // file name
    tOutputFileName = "Column_Buckling_" + tStressType + "_PolyOrd_" + tOrder + "_" + std::to_string( tDim ) + "D.exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Column Buckling in %iD with interpolation order %s using %s stress on %i Processors.",
            tDim,
            tOrder.c_str(),
            tStressType.c_str(),
            par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    if ( par_rank() == 0 )
    {
        check_results( tOutputFileName, tTestCase );
    }

    //----------------------------------------------------------------------------------------

    // test case index
    tTestCase = par_size() == 1 ? 2 : 3;

    // stress used in geometric stiffness (PK2 or Linear)
    tStressType = "PK2";

    // file name
    tOutputFileName = "Column_Buckling_" + tStressType + "_PolyOrd_" + tOrder + "_" + std::to_string( tDim ) + "D.exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Column Buckling in %iD with interpolation order %s using %s stress on %i Processors.",
            tDim,
            tOrder.c_str(),
            tStressType.c_str(),
            par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    if ( par_rank() == 0 )
    {
        check_results( tOutputFileName, tTestCase );
    }

    //----------------------------------------------------------------------------------------

    // test case index
    tTestCase = par_size() == 1 ? 4 : 5;

    // interpolation order
    tOrder = "2";

    // stress used in geometric stiffness (PK2 or Linear)
    tStressType = "Linear";

    // file name
    tOutputFileName = "Column_Buckling_" + tStressType + "_PolyOrd_" + tOrder + "_" + std::to_string( tDim ) + "D.exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Column Buckling in %iD with interpolation order %s using %s stress on %i Processors.",
            tDim,
            tOrder.c_str(),
            tStressType.c_str(),
            par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    if ( par_rank() == 0 )
    {
        check_results( tOutputFileName, tTestCase );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "Column_Buckling_3D",
        "[moris],[example],[structure],[eigen]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Column_Buckling.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // test case index
    uint tTestCase = par_size() == 1 ? 6 : 7;

    // interpolation order
    tOrder = "1";

    // spatial dimensions
    tDim = 3;

    // stress used in geometric stiffness (PK2 or Linear)
    tStressType = "Linear";

    // file name
    tOutputFileName = "Column_Buckling_" + tStressType + "_PolyOrd_" + tOrder + "_" + std::to_string( tDim ) + "D.exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Column Buckling in %iD with interpolation order %s using %s stress on %i Processors.",
            tDim,
            tOrder.c_str(),
            tStressType.c_str(),
            par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    if ( par_rank() == 0 )
    {
        check_results( tOutputFileName, tTestCase );
    }
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "Column_Buckling_Immersed_2D",
        "[moris],[example],[structure],[eigen]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Column_Buckling_Immersed.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // test case index
    uint tTestCase = par_size() == 1 ? 8 : 9;

    // interpolation order
    tOrder = "1";

    // spatial dimensions
    tDim = 2;

    // stress used in geometric stiffness (PK2 or Linear)
    tStressType = "Linear";

    // file name
    tOutputFileName = "Column_Buckling_Immersed_" + tStressType + "_PolyOrd_" + tOrder + "_" + std::to_string( tDim ) + "D.exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Column Buckling in %iD with interpolation order %s using %s stress on %i Processors.",
            tDim,
            tOrder.c_str(),
            tStressType.c_str(),
            par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    if ( par_rank() == 0 )
    {
        check_results( tOutputFileName, tTestCase );
    }
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "Column_Buckling_Immersed_3D",
        "[moris],[example],[structure],[eigen]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Column_Buckling_Immersed.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // test case index
    uint tTestCase = par_size() == 1 ? 10 : 11;

    // interpolation order
    tOrder = "1";

    // spatial dimensions
    tDim = 3;

    // stress used in geometric stiffness (PK2 or Linear)
    tStressType = "Linear";

    // file name
    tOutputFileName = "Column_Buckling_Immersed_" + tStressType + "_PolyOrd_" + tOrder + "_" + std::to_string( tDim ) + "D.exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Column Buckling in %iD with interpolation order %s using %s stress on %i Processors.",
            tDim,
            tOrder.c_str(),
            tStressType.c_str(),
            par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    if ( par_rank() == 0 )
    {
        check_results( tOutputFileName, tTestCase );
    }
}
#endif
