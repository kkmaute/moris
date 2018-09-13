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

TEST_CASE( "moris::sort", "[linalgebra],[sort]" )
{
    SECTION( "real col vector" )
    {
        moris::Matrix< moris::real, moris::DDRMat > a( 3, 1 );
        moris::Matrix< moris::real, moris::DDRMat > b;

        a( 0 ) = 1.0;
        a( 1 ) = 0.0;
        a( 2 ) = 5.0;

        moris::sort( a ,b );

        CHECK( b( 0 ) == 0.0 );
        CHECK( b( 1 ) == 1.0 );
        CHECK( b( 2 ) == 5.0 );
    }

    SECTION( "uint col vector" )
    {
        moris::Matrix< moris::uint, moris::DDUMat > a( 3, 1 );
        moris::Matrix< moris::uint, moris::DDUMat > b;

        a( 0 ) = 1;
        a( 1 ) = 0;
        a( 2 ) = 5;

        moris::sort( a ,b );

        CHECK( b( 0 ) == 0 );
        CHECK( b( 1 ) == 1 );
        CHECK( b( 2 ) == 5 );
    }

}
