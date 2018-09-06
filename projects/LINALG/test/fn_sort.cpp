/*
 * fn_sort.cpp
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
#include "fn_sort.hpp"
#include "fn_trans.hpp"
#include "op_times.hpp"
namespace moris
{
TEST_CASE( "moris::sort", "[linalgebra],[sort]" )
    {
    Matrix< real, DDRMat > a( 3, 1 );
    Matrix< real, DDRMat > b;

    a( 0, 0 ) = 1.0;
    a( 1, 0 ) = 0.0;
    a( 2, 0 ) = 5.0;

    sort( a ,b);

    CHECK( equal_to( b( 0, 0 ), 0.0 ) );
    CHECK( equal_to( b( 1, 0 ), 1.0 ) );
    CHECK( equal_to( b( 2, 0 ), 5.0 ) );

    }
}
