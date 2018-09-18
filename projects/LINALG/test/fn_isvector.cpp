/*
 * fn_isvector.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

// Third-party header files.
#include <catch.hpp>
#include "fn_equal_to.hpp" // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "fn_isvector.hpp"

namespace moris
{
TEST_CASE( "moris::isvector", "[linalgebra],[isvector]" )
    {
    Matrix< DDRMat > a( 3, 3 );
    Matrix< DDRMat > b( 1, 3 );
    Matrix< DDRMat > c( 1, 1 );
    Matrix< DDRMat > d;

    bool tIsVector_1 = isvector( a );
    bool tIsVector_2 = isvector( b );
    bool tIsVector_3 = isvector( b );
    bool tIsVector_4 = isvector( d );

    CHECK( equal_to( tIsVector_1, false ) );
    CHECK( equal_to( tIsVector_2, true ) );
    CHECK( equal_to( tIsVector_3, true ) );
    CHECK( equal_to( tIsVector_4, false ) );
    }
}
