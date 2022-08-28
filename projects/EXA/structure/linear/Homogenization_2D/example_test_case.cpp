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
#include "cl_MTK_Exodus_IO_Helper.hpp"  // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

extern "C"
void check_linear_results(moris::mtk::Exodus_IO_Helper & aExoIO,uint aInnerNodeId, uint aOuterNodeId)
{

    std::cout<<"E_Inner "<<             aExoIO.get_nodal_field_value( aInnerNodeId, 5, 0 ) << std::endl;
    std::cout<<"E_Outer "<<             aExoIO.get_nodal_field_value( aOuterNodeId, 4, 0 ) << std::endl;

    if (gPrintReferenceValues)
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aInnerNodeId ), "Coordinates of reference point");
        moris::print( aExoIO.get_nodal_coordinate( aOuterNodeId ), "Coordinates of reference point");

        std::cout<<"E_Inner: "<<               aExoIO.get_nodal_field_value(aInnerNodeId, 4, 0 ) << std::endl;
        std::cout<<"E_Outer: "<<               aExoIO.get_nodal_field_value(aInnerNodeId, 5, 0 ) << std::endl;
        std::cout<<"C_1111_Inner "<<             aExoIO.get_nodal_field_value( aInnerNodeId, 6, 0 ) << std::endl;
        std::cout<<"C_1111_Outer: "<<         aExoIO.get_nodal_field_value( aInnerNodeId, 7, 0 ) << std::endl;
        return;
    }

    // define error threshold
    real tEpsilon = 1.0e-2;

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinateInner = { {0.5},{0.5} };
    Matrix< DDRMat > tReferenceCoordinateOuter = { {0.5},{0.1} };

    // check nodal coordinates
    real tRelDiffNormInner = moris::norm( aExoIO.get_nodal_coordinate( aInnerNodeId ) - tReferenceCoordinateInner )/ moris::norm(tReferenceCoordinateInner);
    real tRelDiffNormOuter = moris::norm( aExoIO.get_nodal_coordinate( aOuterNodeId ) - tReferenceCoordinateOuter )/ moris::norm(tReferenceCoordinateOuter);
    REQUIRE( tRelDiffNormInner <  tEpsilon );
    REQUIRE( tRelDiffNormOuter <  tEpsilon );

    // check reference displacements
    real tReferenceHomogenizedC1111 = 1.14;
    real tHomogenizedC1111 = aExoIO.get_global_variable(0, 0) + aExoIO.get_global_variable(1, 0);
    REQUIRE( std::abs( tHomogenizedC1111 - tReferenceHomogenizedC1111 ) <  tEpsilon );

}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Homogenization_2D.exo",0,false,false);

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

            std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

    if (gPrintReferenceValues)
    {
        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;
    }
    else
    {
        REQUIRE( tNumDims  ==  2   );
        REQUIRE( tNumNodes ==  337 );
        REQUIRE( tNumElems ==  256 );
    }

    // check results at node with most deformation
    uint tInnerNodeId = 10;
    uint tOuterNodeId = 140;

    check_linear_results(tExoIO,tInnerNodeId,tOuterNodeId);
}

//---------------------------------------------------------------

TEST_CASE("Axisymmetric_Problem_Linear",
        "[moris],[example],[structure],[homogenization]")
{
    // temporary
    gLogger.initialize( "Log.log" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Homogenization_2D.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            check_linear_results_serial();
            break;
        }
        default:
        {
            MORIS_ERROR(false,"This 2D Example can only be run in serial.",par_size());
        }
    }
}

