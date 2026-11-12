/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eig_gen.cpp
 *
 */

#include <catch.hpp>
#include "fn_equal_to.hpp"    // ALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include "fn_eig_gen.hpp"

namespace moris
{
    TEST_CASE( "moris::eig_gen", "[linalgebra],[eig_gen]" )
    {
        // General solvers may order eigenpairs differently and choose any phase.
        Matrix< DDCMat > tA = { { -1.0, 1.0, 3.0 }, { 1.0, 2.0, 0.0 }, { 3.0, 0.0, 2.0 } };
        Matrix< DDCMat > tVectors, tValues;
        eig_gen( tValues, tVectors, tA );

        for ( double tExpected : { -3.0, 2.0, 4.0 } )
        {
            uint tMatches = 0;
            for ( uint j = 0; j < 3; ++j )
            {
                tMatches += std::abs( tValues( j ) - tExpected ) < 1e-12;
            }
            REQUIRE( tMatches == 1 );
        }
        for ( uint j = 0; j < 3; ++j )
        {
            for ( uint i = 0; i < 3; ++i )
            {
                std::complex< double > tResidual = -tValues( j ) * tVectors( i, j );
                for ( uint k = 0; k < 3; ++k )
                {
                    tResidual += tA( i, k ) * tVectors( k, j );
                }
                REQUIRE( std::abs( tResidual ) < 1e-12 );
            }
        }
    }
}    // namespace moris
