/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Debug.cpp
 *
 */

#include <catch.hpp>
#include "algorithms.hpp"
#include "cl_Debug.hpp" // TOL/src/
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
TEST_CASE("moris::tools::Debug",
          "[moris],[tools],[Debug]")
{
    // Checking for duplicates in a Coord list
    Matrix< DDRMat > tCoord = {{0,0,0},{1,1,0},{2,3,0},{3,4,0},{4,4,0},{2,3,0},{2.99,4,0}};  // can be interpreted as x,y,z coords

    Matrix< DDUMat > duplicate_list = moris::Debug::duplicate_row_check(tCoord);
}
}

