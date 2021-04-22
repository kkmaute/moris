/*
 * fn_cross.cpp
 *
 *  Created on: Jan 16, 2019
 *      Author: doble
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
