/*
 * UT_Bounding_Box.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */

#include "catch.hpp"
#include "cl_Matrix.hpp"

#include "fn_bounding_box.hpp"

namespace moris
{
namespace ge
{

TEST_CASE("functions_test_01","[GE],[bounding_box_test]")
{
    moris::Matrix<DDRMat> tPoint(1,2);
    tPoint(0)=1.0;    tPoint(1)=1.0;

    moris::Matrix<DDRMat> tPointA(1,2);
    tPointA(0)=3.0; tPointA(1)=3.0;

    moris::Matrix<DDRMat> tPointB(1,2);
    tPointB(0)=0.0; tPointB(1)=0.0;

    bool tTest2D = check_if_in_bounding_box( tPoint,tPointA,tPointB,5 );

    CHECK(!tTest2D);

    tPointA(0)=-10.0; tPointA(1)=-10.0;
    tPointB(0)=20.0; tPointB(1)=20.0;

    tTest2D = check_if_in_bounding_box( tPoint,tPointA,tPointB,5 );

    CHECK(tTest2D);
//-----------------------------------------------------------------------------
    moris::Matrix<DDRMat> tPoint2(1,3);
    tPoint2(0)=0.0;    tPoint2(1)=0.0;    tPoint2(2)=0.0;

    moris::Matrix<DDRMat> tPointA2(1,3);
    tPointA2(0)=3.0; tPointA2(1)=3.0; tPointA2(2)=3.0;

    moris::Matrix<DDRMat> tPointB2(1,3);
    tPointB2(0)=-2.0; tPointB2(1)=-2.0;   tPointB2(2)=-2.0;

    bool tTest3D = check_if_in_bounding_box( tPoint2,tPointA2,tPointB2,5 );

    CHECK(!tTest3D);

    tPointA2(0)=-10.0; tPointA2(1)=-10.0;   tPointA2(2)=-10.0;
    tPointB2(0)=20.0; tPointB2(1)=20.0;     tPointB2(2)=20.0;

    tTest3D = check_if_in_bounding_box( tPoint2,tPointA2,tPointB2,5 );

    CHECK(tTest2D);
}

}   // end ge namespace
}   // end moris namespace




