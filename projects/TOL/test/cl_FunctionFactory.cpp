/*
 * cl_FunctionFactory.cpp
 *
 *  Created on: Nov 1, 2016
 *      Author: doble
 */
#include <catch.hpp>
#include "algorithms.hpp"

#include "ios.hpp"
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "cl_FunctionFactory.hpp" // TOL/src/
#include "cl_Enums.hpp" // TOL/src/
#include "cl_Function_Levelset.hpp" // TOL/src/

/*
 * This test class tests the function generation using the function factory and verifies the functions are correct
 * - Function type LEVELSET_SPHERE
 */

TEST_CASE(
        "moris::tools::FunctionFactory",
        "[moris],[tools],[FunctionFactory]")
{
    // Tests the factories ability to create a function correctly
    //initialize
    moris::tools::FunctionFactory tFcnFactory;
    moris::tools::Function* SphereLevelsetTest;

    //Specify parameters
    moris::real radius  =  1.0;
    moris::real xcenter =  1.0;
    moris::real ycenter =  2.0;
    moris::real zcenter = -1.0;

    //Create parameter list
    moris::Mat<moris::real> FunctionParameters(4,1);
    FunctionParameters(0,0) = radius;
    FunctionParameters(1,0) = xcenter;
    FunctionParameters(2,0) = ycenter;
    FunctionParameters(3,0) = zcenter;

    //Specify function type
    FunctionType TestFunctionType = FunctionType::LEVELSET_SPHERE;

    //Create Levelset Function using factory
    SphereLevelsetTest = tFcnFactory.create_explicit_function(TestFunctionType, FunctionParameters);

    //Test Functions Capabilities
    moris::Mat<moris::real>  tNodeCoords(1,3);
    tNodeCoords(0,0) = 1.5;  tNodeCoords(0,1) = 2; tNodeCoords(0,2) = -1;
    moris::real testLSV = -0.75;

    moris::real tempLSV;
    tempLSV = SphereLevelsetTest->evaluate_function(tNodeCoords);
    REQUIRE( moris::equal_to(tempLSV,-0.75));

    moris::Mat<moris::real> tDxDp;
    tNodeCoords(0,0) = 1.5;  tNodeCoords(0,1) = -2; tNodeCoords(0,2) = -1.5;
    tDxDp = SphereLevelsetTest->evaluate_dxdp(tNodeCoords);

//    REQUIRE(moris::equal_to(tDxDp(0,0),-.2560737598657920));
//    REQUIRE(moris::equal_to(tDxDp(0,1), .4264014327112208));
//    REQUIRE(moris::equal_to(tDxDp(0,2),-.2169304578186562));
//    REQUIRE(moris::equal_to(tDxDp(1,0),1));
//    REQUIRE(moris::equal_to(tDxDp(1,1),0));
//    REQUIRE(moris::equal_to(tDxDp(1,2),0));
//    REQUIRE(moris::equal_to(tDxDp(2,0),0));
//    REQUIRE(moris::equal_to(tDxDp(2,1),1));
//    REQUIRE(moris::equal_to(tDxDp(2,2),0));
//    REQUIRE(moris::equal_to(tDxDp(3,0),0));
//    REQUIRE(moris::equal_to(tDxDp(3,1),0));
//    REQUIRE(moris::equal_to(tDxDp(3,2),1));
    delete SphereLevelsetTest;


}


