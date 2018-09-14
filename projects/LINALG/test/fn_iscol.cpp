/*
 * fn_iscol.cpp
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
#include "fn_iscol.hpp"

namespace moris
{
TEST_CASE( "moris::iscol", "[linalgebra],[iscol]" )
    {
    Matrix< DDRMat > a( 3, 3 );
    Matrix< DDRMat > b( 3, 1 );
    Matrix< DDRMat > c( 1, 1 );
    Matrix< DDRMat > d;

    bool tIsCol_1 = iscol( a );
    bool tIsCol_2 = iscol( b );
    bool tIsCol_3 = iscol( c );
    bool tIsCol_4 = iscol( d );

    CHECK( equal_to( tIsCol_1, false ) );
    CHECK( equal_to( tIsCol_2, true ) );
    CHECK( equal_to( tIsCol_3, true ) );
    CHECK( equal_to( tIsCol_4, false ) );
    }
}
