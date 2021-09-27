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

        std::cout<<"U1: "<<               aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;
        std::cout<<"U2: "<<               aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) << std::endl;
        std::cout<<"S_VM: "<<             aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) << std::endl;
        std::cout<<"S_11 (x): "<<         aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) << std::endl;
        std::cout<<"S_22 (radial): "<<    aExoIO.get_nodal_field_value( aNodeId, 6, 0 ) << std::endl;
        std::cout<<"S_33 (azimuthal): "<< aExoIO.get_nodal_field_value( aNodeId, 7, 0 ) << std::endl;
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
    Matrix< DDRMat > tReferenceU = { {1.2},{-1.4} };
    Matrix< DDRMat > tU =  {{ aExoIO.get_nodal_field_value(aNodeId, 2, 0)},{aExoIO.get_nodal_field_value(aNodeId, 3, 0) }};
    REQUIRE( moris::norm( tU - tReferenceU ) <  tEpsilon );

    // check stresses.  Since this is a solid cylinder, the radial and azimuthal stress should be equal to the pressure (neumann traction)
    // pressure = 1.0 so radial and theta stress should be -1.0
    // https://en.wikipedia.org/wiki/Cylinder_stress
    real tReference_S_x     = 0.0;
    real tReference_S_rad   = -1.0;
    real tReference_S_theta = -1.0;
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) - tReference_S_x     < tEpsilon);
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 6, 0 ) - tReference_S_rad   < tEpsilon);
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 7, 0 ) - tReference_S_theta < tEpsilon);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_stress(moris::mtk::Exodus_IO_Helper & aExoIO,uint aNodeId)
{
    if (gPrintReferenceValues)
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point");

        std::cout<<"U1: "<<               aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;
        std::cout<<"U2: "<<               aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) << std::endl;
        std::cout<<"S_VM: "<<             aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) << std::endl;
        std::cout<<"S_11 (x): "<<         aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) << std::endl;
        std::cout<<"S_22 (radial): "<<    aExoIO.get_nodal_field_value( aNodeId, 6, 0 ) << std::endl;
        std::cout<<"S_33 (azimuthal): "<< aExoIO.get_nodal_field_value( aNodeId, 7, 0 ) << std::endl;
        return;
    }

    // define error threshold
    real tEpsilon = 1.0e-5;

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { {1.0},{1.0} };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);
    REQUIRE( tRelDiffNorm <  tEpsilon );

    // check stresses.  Since this is a solid cylinder, the radial and azimuthal stress should be equal to the pressure (neumann traction)
    // pressure = 1.0 so radial and theta stress should be -1.0
    // https://en.wikipedia.org/wiki/Cylinder_stress
    real tReference_S_x     = 0.0;
    real tReference_S_rad   = -1.0;
    real tReference_S_theta = -1.0;
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) - tReference_S_x     < tEpsilon);
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 6, 0 ) - tReference_S_rad   < tEpsilon);
    REQUIRE(  aExoIO.get_nodal_field_value( aNodeId, 7, 0 ) - tReference_S_theta < tEpsilon);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Axisymmetric_Problem.exo",0,false,false);
    moris::mtk::Exodus_IO_Helper tExoIO_thermal("Axisymmetric_Problem_ConstantThermalField.exo",0,false,false);

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
        REQUIRE( tNumDims  ==  2   );
        REQUIRE( tNumNodes ==  9 );
        REQUIRE( tNumElems ==  4 );
    }

    // check results at node with most deformation
    uint tNodeId = 7;

    check_linear_results(tExoIO,tNodeId);
    check_linear_results_stress(tExoIO_thermal,tNodeId);
}

//---------------------------------------------------------------

TEST_CASE("Axisymmetric_Problem_Linear",
        "[moris],[example],[structure],[axisymmetric]")
{
    // temporary
    gLogger.initialize( "Log.log" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Axisymmetric_Problem.so";

    char * argv[2] = {tString1,tString2};

    char tString1_thermal[] = "";
    char tString2_thermal[] = "./Axisymmetric_Problem_ConstantThermalField.so";

    char * argv_thermal[2] = {tString1_thermal,tString2_thermal};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );
    int tRet_thermal = fn_WRK_Workflow_Main_Interface( argc, argv_thermal );

    // catch test statements should follow
    REQUIRE( tRet         ==  0 );
    REQUIRE( tRet_thermal ==  0 );

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

