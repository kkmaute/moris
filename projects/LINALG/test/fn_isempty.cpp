/*
 * fn_isempty.cpp
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
#include "fn_isempty.hpp"

namespace moris
{
TEST_CASE( "moris::isempty", "[linalgebra],[isempty]" )
    {
    Matrix< real, DDRMat > a( 3, 3 );
    Matrix< real, DDRMat > b;
    Matrix< real, DDRMat > c( 0, 3 );

    bool tIsEmpty_1 = isempty( a );
    bool tIsEmpty_2 = isempty( b );
    bool tIsEmpty_3 = isempty( b );

    CHECK( equal_to( tIsEmpty_1, false ) );
    CHECK( equal_to( tIsEmpty_2, true ) );
    CHECK( equal_to( tIsEmpty_3, true ) );
    }
}
