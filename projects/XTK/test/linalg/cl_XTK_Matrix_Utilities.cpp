/*
 * cl_XTK_Matrix_Utilities.cpp
 *
 *  Created on: Jan 16, 2018
 *      Author: doble
 */

#include "catch.hpp"

// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"


namespace xtk
{

TEST_CASE("Has Duplicate Rows","[LINALG][DUP_ROWS]")
        {
    Mat<size_t,Default_Matrix_Integer> tTestMat({{1,2,3},{4,5,6},{7,8,9},{10,11,12}});
    CHECK(check_for_duplicate_rows(tTestMat, false));
        }
}
