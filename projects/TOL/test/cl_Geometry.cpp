/*
 * cl_Geometry.cpp
 *
 *  Created on: Apr 7, 2017
 *      Author: gleim
 */

#include <catch.hpp>
#include "algorithms.hpp"
#include "cl_Geometry.hpp" // TOL/src/
#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"

namespace moris
{
TEST_CASE("moris::tools::Geometry",
          "[moris],[tools],[Geometry]")
        {
            // Volume of a tetrahedron
            Matrix< DDRMat > tCoord = {{0,0,0},{1,-1,0},{2,0,0},{1,0,1}};  // can be interpreted as x,y,z coords

            real volume = moris::Geometry::vol_tetrahedron(tCoord);
        }
}





