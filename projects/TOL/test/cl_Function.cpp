/*
 * cl_Function.cpp
 *
 *  Created on: Nov 1, 2016
 *      Author: doble
 */




#include <catch.hpp>

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
        "moris::tools::FunctionLevelset",
        "[moris],[tools],[Function],[AnalyticalLevelset]")
{
    // Sphere test
}
