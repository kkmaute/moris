/*
 * cl_Geometry.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: gleim
 */

#include <catch.hpp>
#include "algorithms.hpp"
#include "cl_Geometry.hpp" // TOL/src/
#include "cl_Mat.hpp" // LNA/src

TEST_CASE("moris::tools::Geometry",
          "[moris],[tools],[Geometry]")
        {
            // Volume of a tetrahedron
            moris::Mat<moris::real> tCoord = {{0,0,0},{1,-1,0},{2,0,0},{1,0,1}};  // can be interpreted as x,y,z coords

            moris::real volume = moris::Geometry::vol_tetrahedron(tCoord);
        }






