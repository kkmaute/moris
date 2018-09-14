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
TEST_CASE( "sort", "[linalgebra],[sort]" )
{
    SECTION( "real col vector" )
    {
        Matrix< DDRMat > a( 3, 1 );
        Matrix< DDRMat > b;

        a( 0 ) = 1.0;
        a( 1 ) = 0.0;
        a( 2 ) = 5.0;

        sort( a ,b );

        CHECK( b( 0 ) == 0.0 );
        CHECK( b( 1 ) == 1.0 );
        CHECK( b( 2 ) == 5.0 );
    }

    SECTION( "uint col vector" )
    {
        Matrix< DDUMat > a( 3, 1 );
        Matrix< DDUMat > b;

        a( 0 ) = 1;
        a( 1 ) = 0;
        a( 2 ) = 5;


        sort( a ,b );

        CHECK( b( 0 ) == 0 );
        CHECK( b( 1 ) == 1 );
        CHECK( b( 2 ) == 5 );
    }

}
}
