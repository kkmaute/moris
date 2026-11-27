/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_mat_to_vec.cpp
 *
 */

#include "fn_mat_to_vec.hpp"
#include "fn_vectorize.hpp"
#include <catch.hpp>
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"

namespace moris
{
    namespace
    {
        template< typename Matrix_Type >
        void check_flattening( Matrix< Matrix_Type > aMatrix )
        {
            // Values distinguish row/column order for rectangular and fixed inputs.
            for ( uint iRow = 0; iRow < aMatrix.n_rows(); ++iRow )
            {
                for ( uint iColumn = 0; iColumn < aMatrix.n_cols(); ++iColumn )
                {
                    aMatrix( iRow, iColumn ) = 10.0 * iRow + iColumn;
                }
            }

            const auto tColumnOrder = mat_to_vec( aMatrix );
            const auto tRowOrder    = mat_to_vec( aMatrix, true );
            REQUIRE( tColumnOrder.size() == aMatrix.numel() );
            REQUIRE( tRowOrder.size() == aMatrix.numel() );

            for ( uint iRow = 0; iRow < aMatrix.n_rows(); ++iRow )
            {
                for ( uint iColumn = 0; iColumn < aMatrix.n_cols(); ++iColumn )
                {
                    CHECK( tColumnOrder( iColumn * aMatrix.n_rows() + iRow ) == Approx( 10.0 * iRow + iColumn ) );
                    CHECK( tRowOrder( iRow * aMatrix.n_cols() + iColumn ) == Approx( 10.0 * iRow + iColumn ) );
                }
            }
        }
    }    // namespace

    TEST_CASE( "mat_to_vec fixed, dynamic and empty inputs", "[linalgebra][mat_to_vec]" )
    {
        // Direct flattening must work for fixed Eigen dimensions, both orderings,
        // vectors and empty matrices without constructing a reshaped input type.
        check_flattening( Matrix< DDRMat >( 4, 2 ) );
        check_flattening( Matrix< F33RMat >( 3, 3 ) );
        check_flattening( Matrix< DDRMat >( 1, 4 ) );
        check_flattening( Matrix< DDRMat >( 4, 1 ) );
        check_flattening( Matrix< DDRMat >( 0, 0 ) );
        check_flattening( Matrix< DDRMat >( 0, 3 ) );
        check_flattening( Matrix< DDRMat >( 3, 0 ) );

        const Matrix< DDRMat > tMatrix   = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Vector< sint >   tIntegers = mat_to_vec< DDRMat, sint >( tMatrix, true );
        REQUIRE( tIntegers.size() == 4 );
        CHECK( tIntegers( 0 ) == 1 );
        CHECK( tIntegers( 1 ) == 2 );
        CHECK( tIntegers( 2 ) == 3 );
        CHECK( tIntegers( 3 ) == 4 );
    }

    TEST_CASE( "Flattening header before vectorization header", "[linalgebra][mat_to_vec][vectorize]" )
    {
        // Keep these public headers first in this translation unit. A shared
        // include guard used to suppress vectorize() in this include order.
        const Matrix< DDRMat > tMatrix     = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Matrix< DDRMat > tVectorized = vectorize( tMatrix );
        const auto             tFlattened  = mat_to_vec( tMatrix );
        REQUIRE( tVectorized.n_rows() == 4 );
        REQUIRE( tVectorized.n_cols() == 1 );

        for ( uint i = 0; i < 4; ++i )
        {
            CHECK( tVectorized( i ) == tFlattened( i ) );
        }
    }
}    // namespace moris
