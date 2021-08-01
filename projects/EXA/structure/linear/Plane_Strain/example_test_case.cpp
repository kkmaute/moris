//
// example specific interface to moris
//

#include <catch.hpp>

#include "cl_Logger.hpp" // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"  // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

//---------------------------------------------------------------

// global variable for interpolation order
uint gInterpolationOrder;

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

        std::cout<<"U1: "<<                aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;
        std::cout<<"U2: "<<                aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) << std::endl;
        std::cout<<"S_VM: "<<              aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) << std::endl;
        std::cout<<"S_11 (x): "<<          aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) << std::endl;
        std::cout<<"S_22 (y): "<<          aExoIO.get_nodal_field_value( aNodeId, 6, 0 ) << std::endl;
        std::cout<<"S_33 (transverse): "<< aExoIO.get_nodal_field_value( aNodeId, 7, 0 ) << std::endl;
        return;
    }

    // define error threshold
    real tEpsilon = 1.0e-5;

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { {1.0},{1.0} };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);
    REQUIRE( tRelDiffNorm <  tEpsilon );

    // check reference displacements
    Matrix< DDRMat > tReferenceU = { {0.78},{-1.82} };
    Matrix< DDRMat > tU =  {{ aExoIO.get_nodal_field_value(aNodeId, 2, 0)},{aExoIO.get_nodal_field_value(aNodeId, 3, 0) }};
    REQUIRE( moris::norm( tU - tReferenceU ) <  tEpsilon );

    // check stresses.
    real tReference_S_x =  0.0;
    real tReference_S_y = -1.0;
    real tReference_S_z = -0.3;
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) - tReference_S_x < tEpsilon);
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 6, 0 ) - tReference_S_y < tEpsilon);
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 7, 0 ) - tReference_S_z < tEpsilon);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Plane_Strain_Problem.exo",0,false,false);

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
        REQUIRE( tNumNodes ==  9 );
        REQUIRE( tNumElems ==  4 );
    }

    // check results at node with most deformation
    uint tNodeId = 7;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

TEST_CASE("Plane_Strain_Problem_Linear",
        "[moris],[example],[structure],[plane_strain]")
{
    // temporary
    gLogger.initialize( "Log.log" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Plane_Strain_Problem.so";

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

