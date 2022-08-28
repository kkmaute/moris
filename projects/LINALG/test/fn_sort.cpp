/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort.cpp
 *
 */

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

    SECTION( "DDRMat, ascending, column sorted" )
    {
    	Matrix< DDRMat > a( 3, 2 );
    	Matrix< DDRMat > b;

    	a = {{1.0, 2.0},
    	     {5.0, 0.0},
			 {3.0, 1.0}};

    	sort( a, b, "ascend", 0 );

    	CHECK( b( 0, 0 ) == 1.0 );
    	CHECK( b( 1, 0 ) == 3.0 );
    	CHECK( b( 2, 0 ) == 5.0 );
    	CHECK( b( 0, 1 ) == 0.0 );
    	CHECK( b( 1, 1 ) == 1.0 );
    	CHECK( b( 2, 1 ) == 2.0 );
    }

    SECTION( "DDRMat, descending, column sorted" )
    {
    	Matrix< DDRMat > a( 3, 2 );
    	Matrix< DDRMat > b;

    	a = {{1.0, 2.0},
    	     {5.0, 0.0},
			 {3.0, 1.0}};

    	sort( a, b, "descend", 0 );

    	CHECK( b( 0, 0 ) == 5.0 );
    	CHECK( b( 1, 0 ) == 3.0 );
    	CHECK( b( 2, 0 ) == 1.0 );
    	CHECK( b( 0, 1 ) == 2.0 );
    	CHECK( b( 1, 1 ) == 1.0 );
    	CHECK( b( 2, 1 ) == 0.0 );
    }

    SECTION( "DDRMat, ascending, row sorted" )
    {
    	Matrix< DDRMat > a( 2, 3 );
    	Matrix< DDRMat > b;

    	a = {{1.0, 5.0, 3.0},
    	     {2.0, 0.0, 1.0}};

    	sort( a, b, "ascend", 1 );

    	CHECK( b( 0, 0 ) == 1.0 );
    	CHECK( b( 0, 1 ) == 3.0 );
    	CHECK( b( 0, 2 ) == 5.0 );
    	CHECK( b( 1, 0 ) == 0.0 );
    	CHECK( b( 1, 1 ) == 1.0 );
    	CHECK( b( 1, 2 ) == 2.0 );
    }

    SECTION( "DDRMat, descending, row sorted" )
    {
    	Matrix< DDRMat > a( 2, 3 );
    	Matrix< DDRMat > b;

    	a = {{1.0, 5.0, 3.0},
    	     {2.0, 0.0, 1.0}};

    	sort( a, b, "descend", 1 );

    	CHECK( b( 0, 0 ) == 5.0 );
    	CHECK( b( 0, 1 ) == 3.0 );
    	CHECK( b( 0, 2 ) == 1.0 );
    	CHECK( b( 1, 0 ) == 2.0 );
    	CHECK( b( 1, 1 ) == 1.0 );
    	CHECK( b( 1, 2 ) == 0.0 );
    }
}
}

