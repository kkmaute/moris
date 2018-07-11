/*
 * cl_Interpolation.cpp
 *
 *  Created on: Apr 4, 2017
 *      Author: doble
 */

#include <catch.hpp>
#include "algorithms.hpp"
#include "cl_Interpolation.hpp" // TOL/src/
#include "cl_Mat.hpp" // LNA/src

TEST_CASE("moris::tools::Interpolation",
          "[moris],[tools],[Interpolation]")
        {
            // Linear interpolation
            // To a specified location
            moris::Mat<moris::real> tInterpVars = {{0,1,5},{1,0,7}};  // can be interpreted as x,y,z coords
            moris::Mat<moris::real> tLocation   = {{0}};              // can be interpreted as a local coordinate along an edge

            moris::Mat<moris::real> tVars = moris::Interpolation::linear_interpolation_location(tInterpVars,tLocation);

            REQUIRE(moris::equal_to(tVars(0,0),0.5));
            REQUIRE(moris::equal_to(tVars(0,1),0.5));
            REQUIRE(moris::equal_to(tVars(0,2),6));

            // To a specified value
            moris::Mat<moris::real> tValue = {{0}};
            moris::Mat<moris::real> tLoc = moris::Interpolation::linear_interpolation_value(tInterpVars,tValue);

            REQUIRE(moris::equal_to(tLoc(0,0),-1.0));
            REQUIRE(moris::equal_to(tLoc(0,1),1.0));
            REQUIRE(moris::equal_to(tLoc(0,2),-6.0)); // extrapolation add a throw in the function
        }
