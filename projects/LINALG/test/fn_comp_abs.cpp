/*
 * fn_comp_abs.cpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */


#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_comp_abs.hpp"


TEST_CASE("moris::comp_abs",
          "[linalgebra],[comp_abs]")
{
    moris::Matrix< moris::DDRMat > tMat1(3,3);
    tMat1(0,0) = -1.0;         tMat1(0,1) = 1.0;        tMat1(0,2) = -1.0;
    tMat1(1,0) = -11.0;        tMat1(1,1) = 1.0;        tMat1(1,2) = -3.0;
    tMat1(2,0) = -13.0;        tMat1(2,1) = 3.0;        tMat1(2,2) = -1.0;

    moris::Matrix< moris::DDRMat > tMatAbs1 = moris::comp_abs(tMat1);

    CHECK(moris::equal_to( tMatAbs1(0,0), 1.0));
    CHECK(moris::equal_to( tMatAbs1(0,1), 1.0));
    CHECK(moris::equal_to( tMatAbs1(0,2), 1.0));
    CHECK(moris::equal_to( tMatAbs1(1,0), 11.0));
    CHECK(moris::equal_to( tMatAbs1(1,1), 1.0));
    CHECK(moris::equal_to( tMatAbs1(1,2), 3.0));
    CHECK(moris::equal_to( tMatAbs1(2,0), 13.0));
    CHECK(moris::equal_to( tMatAbs1(2,1), 3.0));
    CHECK(moris::equal_to( tMatAbs1(2,2), 1.0));

}
