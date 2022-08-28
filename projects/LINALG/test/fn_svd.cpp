/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_svd.cpp
 *
 */

#include <catch.hpp>

// MORIS project header files.

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "fn_equal_to.hpp" // ALG/src
#include "linalg_typedefs.hpp"

#include "fn_svd.hpp"
#include "op_times.hpp"
#include "fn_trans.hpp"

TEST_CASE(
        "moris::svd",
        "[linalgebra],[svd]" )
{
    moris::Matrix< moris::DDRMat > tM = { { 1.0, 0.0, -1.0}, { -2.0, 1.0, 4.0 } };

    SECTION( "SVD real small" )
    {
        moris::Matrix< moris::DDRMat > tS;

        moris::svd( tS,tM );

        moris::Matrix< moris::DDRMat > tA( 2, 1);
        tA( 0 ) = std::sqrt( 0.5 * ( 23.0 + sqrt( 505.0 ) ) );
        tA( 1 ) = std::sqrt( 0.5 * ( 23.0 - sqrt( 505.0 ) ) );

        for( uint k=0; k<tS.length(); ++k )
        {
                REQUIRE( moris::equal_to( tA( k ), tS( k ) ) );
        }
    }

    SECTION( "SVD real big" )
        {
            moris::Matrix< moris::DDRMat > tU;
            moris::Matrix< moris::DDRMat > tS;
            moris::Matrix< moris::DDRMat > tV;

            moris::svd( tU, tS, tV, tM );

            moris::Matrix< moris::DDRMat > tSigma( tS.length(), tS.length(), 0.0 );
            for( uint k=0; k<tS.length(); ++k )
            {
                tSigma( k, k ) = tS( k );
            }

            moris::Matrix< moris::DDRMat > tA;
            tA = tU * ( tSigma * moris::trans( tV ) );

            for( uint j=0; j<tM.n_cols(); ++j )
            {
                for( uint i=0; i<tM.n_rows(); ++i )
                {
                    REQUIRE(  moris::equal_to( tA( i, j ), tM( i, j ) ) );
                }
            }
        }
}

