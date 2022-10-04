/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_check_equal.hpp
 *
 */

#ifndef MORIS_FN_GEN_CHECK_EQUAL_HPP
#define MORIS_FN_GEN_CHECK_EQUAL_HPP

#include "cl_Matrix.hpp"
#include "catch.hpp"

namespace moris
{
    namespace ge
    {
        /**
         * Approximate checking for moris matrices being equal.
         *
         * @tparam Matrix_Type Matrix type
         * @param aMatrix1 Comparison matrix 1
         * @param aMatrix2 Comparison matrix 2
         */
        template< typename Matrix_Type >
        void
        check_equal( Matrix< Matrix_Type > aMatrix1, Matrix< Matrix_Type > aMatrix2, real aMargin = MORIS_REAL_EPS )
        {
            REQUIRE( aMatrix1.n_rows() == aMatrix2.n_rows() );
            REQUIRE( aMatrix1.n_cols() == aMatrix2.n_cols() );
            for ( uint tRowIndex = 0; tRowIndex < aMatrix1.n_rows(); tRowIndex++ )
            {
                for ( uint tColumnIndex = 0; tColumnIndex < aMatrix1.n_cols(); tColumnIndex++ )
                {
                    CHECK( aMatrix1( tRowIndex, tColumnIndex ) == Approx( aMatrix2( tRowIndex, tColumnIndex ) ).margin( aMargin ) );

                    if ( aMatrix1( tRowIndex, tColumnIndex ) != Approx( aMatrix2( tRowIndex, tColumnIndex ) ).margin( aMargin ) )
                    {
                        std::cout << aMatrix1( tRowIndex, tColumnIndex ) << " =/= " << aMatrix2( tRowIndex, tColumnIndex ) << std::endl;
                    }
                }
            }
        }
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_FN_GEN_CHECK_EQUAL_HPP
