//
// example specific interface to moris
//

#include <catch.hpp>

#include "cl_Logger.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "HDF5_Tools.hpp"

using namespace moris;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

extern "C"
void check_results_serial(
        std::string aExoFileName,
        uint        aTestCaseIndex)
{

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Test Case %d on %i processor.",aTestCaseIndex,par_size());
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO(aExoFileName.c_str(),0,false,false);

    if (false)
    {
        std::cout << "Test case index: " << aTestCaseIndex << std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell<uint> tReferenceNumDims  = { 2,2};
    Cell<uint> tReferenceNumNodes = {10452,12401};
    Cell<uint> tReferenceNumElems = {8173,6388};

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims(aTestCaseIndex),tNumDims,std::abs((tNumDims-tReferenceNumDims(aTestCaseIndex))/tReferenceNumDims(aTestCaseIndex)*100.0));
    MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes(aTestCaseIndex),tNumNodes,std::abs((tNumNodes-tReferenceNumNodes(aTestCaseIndex))/tReferenceNumNodes(aTestCaseIndex)*100.0));
    MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElems(aTestCaseIndex),tNumElems,std::abs((tNumElems-tReferenceNumElems(aTestCaseIndex))/tReferenceNumElems(aTestCaseIndex)*100.0));

    REQUIRE( tNumDims  ==  tReferenceNumDims(aTestCaseIndex)  );
    REQUIRE( tNumNodes ==  tReferenceNumNodes(aTestCaseIndex) );
    REQUIRE( tNumElems ==  tReferenceNumElems(aTestCaseIndex) );

}

extern "C"
void check_results_parallel(
        std::string aExoFileName,
        uint        aTestCaseIndex)
{

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Test Case %d on %i processor.",aTestCaseIndex,par_size());
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO(aExoFileName.c_str(),0,false,false);

    if (false)
    {
        std::cout << "Test case index: " << aTestCaseIndex << " on proc: "<< par_rank()<< std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell<uint> tReferenceNumDims= { 2,2};
    Cell<uint> tReferenceNumNodes;
    Cell<uint> tReferenceNumElems;

    if( par_rank() == 0 )
    {
        tReferenceNumNodes = {2745,3185};
        tReferenceNumElems = {2035,1590};
    }
    if( par_rank() == 1 )
    {
        tReferenceNumNodes = {2729,3206};
        tReferenceNumElems = {2038,1604};
    }
    if( par_rank() == 2 )
    {
        tReferenceNumNodes = {2557,3059};
        tReferenceNumElems = {2058,1598};
    }
    if( par_rank() == 3 )
    {
        tReferenceNumNodes = {2586,3108};
        tReferenceNumElems = {2046,1596};
    }

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    //std::cout<<"nodes: "<<tNumNodes<<" ele: "<<tNumElems<<" proc: "<<par_rank()<<std::endl;

    MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims(aTestCaseIndex),tNumDims,std::abs((tNumDims-tReferenceNumDims(aTestCaseIndex))/tReferenceNumDims(aTestCaseIndex)*100.0));
    MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes(aTestCaseIndex),tNumNodes,std::abs((tNumNodes-tReferenceNumNodes(aTestCaseIndex))/tReferenceNumNodes(aTestCaseIndex)*100.0));
    MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElems(aTestCaseIndex),tNumElems,std::abs((tNumElems-tReferenceNumElems(aTestCaseIndex))/tReferenceNumElems(aTestCaseIndex)*100.0));

    REQUIRE( tNumDims  ==  tReferenceNumDims(aTestCaseIndex)  );
    REQUIRE( tNumNodes ==  tReferenceNumNodes(aTestCaseIndex) );
    REQUIRE( tNumElems ==  tReferenceNumElems(aTestCaseIndex) );

}

//---------------------------------------------------------------

TEST_CASE("Leveset Boxbeam Adaptive Refinement",
        "[moris],[example],[optimization],[levelset_boxbeam_adaptive_refinement]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam_Adaptive_Refinement.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );

    // set test case index
    uint tTestCaseIndex = 0;

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            // perform check for Test Case 0
            check_results_serial("Levelset_Boxbeam_Adaptive_Refinement.exo.e-s.0015",tTestCaseIndex);
            break;
        }
        case 4:
        {
            // perform check
            check_results_parallel( "Levelset_Boxbeam_Adaptive_Refinement.exo.e-s.0015",tTestCaseIndex );

            break;
        }
        default:
        {
            MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
        }
    }
}

//---------------------------------------------------------------

TEST_CASE("Leveset Boxbeam Adaptive Refinement Quadratic",
        "[moris],[example],[optimization],[levelset_boxbeam_adaptive_refinement_quadratic]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam_Adaptive_Refinement_Quadratic.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );

    // set test case index
    uint tTestCaseIndex = 1;

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            // perform check for Test Case 0
            check_results_serial("Levelset_Boxbeam_Adaptive_Refinement_Quadratic.exo.e-s.0015",tTestCaseIndex);
            break;
        }
        case 4:
        {

            // perform check
            check_results_parallel( "Levelset_Boxbeam_Adaptive_Refinement_Quadratic.exo.e-s.0015",tTestCaseIndex );
            break;
        }
        default:
        {
            MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
        }
    }
}


