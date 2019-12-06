/*
 * UT_Bounding_Box.cpp
 *
 *  Created on: Oct 22, 2019
 *      Author: sonne
 */

#include "catch.hpp"
#include "cl_Matrix.hpp"

#include "../src/geometry/fn_bounding_box.hpp"

namespace moris
{
namespace ge
{

TEST_CASE( "bounding_box_test_01","[GE],[bounding_box_test_01]" )
{
    moris::Matrix<DDRMat> tPoint = {{ 1.0, 1.0 }};

    moris::Matrix<DDRMat> tPointA = {{ 3.0, 2.0 }};

    moris::Matrix<DDRMat> tPointB = {{ 0.0, -1.0}};

    bool tTest2D = outside_bounding_box_check( tPoint,tPointA,tPointB,5 );

    CHECK(!tTest2D);

    tPointA(0)=-10.0; tPointA(1)= 0.0;
    tPointB(0)= 20.0; tPointB(1)= 1.0;

    tTest2D = outside_bounding_box_check( tPoint,tPointA,tPointB,5 );

    CHECK(tTest2D);
//-----------------------------------------------------------------------------
    moris::Matrix<DDRMat> tPoint2 = {{ 0.0, 0.0, 0.0 }};

    moris::Matrix<DDRMat> tPointA2 = {{ 3.0, 0.0, 4.0 }};

    moris::Matrix<DDRMat> tPointB2 = {{ -2.0, 1.0, -3.0 }};

    bool tTest3D = outside_bounding_box_check( tPoint2,tPointA2,tPointB2,5 );

    CHECK(!tTest3D);

    tPointA2(0)=-10.0; tPointA2(1)= 0.0;   tPointA2(2)= 15.0;
    tPointB2(0)= 4.0; tPointB2(1)= 20.0;     tPointB2(2)=-20.0;

    tTest3D = outside_bounding_box_check( tPoint2,tPointA2,tPointB2,5 );

    CHECK(tTest3D);
}
//-----------------------------------------------------------------------------

}   // end ge namespace
}   // end moris namespace




