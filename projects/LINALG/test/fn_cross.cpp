/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cross.cpp
 *
 */

#include <catch.hpp>
#include "op_equal_equal.hpp" //ALG
#include "fn_all_true.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_cross.hpp"

namespace moris
{
TEST_CASE(
        "cross",
        "[linalgebra],[cross]" )
        {
            Matrix<F31RMat> tVec1 = {{1.0},{2.0},{3.0}};
            Matrix<F31RMat> tVec2 = {{3.0},{2.0},{1.4}};

            Matrix<F31RMat> tGoldCrossProd = {{-3.2}, {7.6}, {-4.0}};

            Matrix<F31RMat> tCrossProd = cross(tVec1,tVec2);

            CHECK(all_true(tGoldCrossProd == tGoldCrossProd));

            tCrossProd = cross(tVec1,tVec2);
            CHECK(all_true(tGoldCrossProd == tGoldCrossProd));

            tCrossProd = cross(tVec1,tVec2);
            CHECK(all_true(tGoldCrossProd == tGoldCrossProd));

            tCrossProd = cross(tVec1,tVec2);
            CHECK(all_true(tGoldCrossProd == tGoldCrossProd));

        }
}

