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

#include "cl_Logger.hpp"                // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"  // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

extern "C"
void check_linear_results(moris::mtk::Exodus_IO_Helper & aExoIO,uint aNodeId)
{
    if (gPrintReferenceValues)
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point");

        std::cout<<"Ux: "<<          aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;
        std::cout<<"Uy: "<<          aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) << std::endl;
        std::cout<<"Temperature: "<< aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) << std::endl;
        return;
    }

    // define error threshold.
    // larger deflection epsilon used since the
    real tEpsilon            = 1.0e-6;
    real tEpsilon_Deflection = 1.0e-3;

    // check nodal coordinates
    Matrix< DDRMat > tReferenceCoordinate = { {10.0},{0.5} };
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);
    REQUIRE( tRelDiffNorm <  tEpsilon );

    // check beam deflection from the temperature field
    // see matlab script or https://nptel.ac.in/content/storage2/courses/105101085/downloads/lec-26.pdf
    real tReferenceUy = -2.0;
    real tUy =  aExoIO.get_nodal_field_value( aNodeId, 3, 0 );
    REQUIRE( ( tUy - tReferenceUy ) <  tEpsilon_Deflection );

    // check to see if the nodal value at the coordinate location matches
    real tUy_coords = aExoIO.get_nodal_field_value_by_coords(tReferenceCoordinate, 3, 0);
    REQUIRE( ( tUy_coords - tReferenceUy ) <  tEpsilon_Deflection );

    // check field temperature at this node
    real tTemperatureField = 1.0;
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) - tTemperatureField < tEpsilon);

    // check to see if the nodal value at the coordinate location matches
    REQUIRE(  aExoIO.get_nodal_field_value_by_coords( tReferenceCoordinate, 4, 0 ) - tTemperatureField < tEpsilon);

}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Beam_Temperature_Field_Problem.exo",0,false,false);

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if (gPrintReferenceValues)
    {
        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;
    }
    else
    {
        REQUIRE( tNumDims  ==  2 );
        REQUIRE( tNumNodes ==  3975 );
        REQUIRE( tNumElems ==  924 );
    }

    // check results at node with most deformation
    uint tNodeId = 3973;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

TEST_CASE("Field_example_write",
        "[moris],[example],[structure],[Plate_Temperature_Field_Write]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Beam_Temperature_Field.so";

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

