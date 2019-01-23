/*
 * fn_bubble_sort.hpp
 *
 *  Created on: Jan 10, 2018
 *      Author: doble
 */



#include "tools/fn_bubble_sort.hpp"

#include "catch.hpp"

// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"


namespace xtk
{
TEST_CASE("Row Bubble Sort","[ROW_BUBBLE]")
{
    moris::Matrix< moris::DDSTMat > tIntMat({{3,2,1,4},{4,3,2,99},{9,5,10,10}});
    row_bubble_sort(tIntMat);
    moris::Matrix< moris::DDSTMat > tIntExpectedMat({{1,2,3,4},{2,3,4,99},{5,9,10,10}});

    CHECK(xtk::equal_to(tIntMat, tIntExpectedMat));

    moris::Matrix< moris::DDRMat > tRealMat({{3.0,2.0,1.0,4.0},{2.2,3.1,4.2,99.0},{9.1,5.1,10.3,10.2}});
    row_bubble_sort(tRealMat);
    moris::Matrix< moris::DDRMat > tRealExpectedMat({{1.0,2.0,3.0,4.0},{2.2,3.1,4.2,99.0},{5.1,9.1,10.2,10.3}});

    CHECK(xtk::equal_to(tRealMat, tRealExpectedMat));
}
}
