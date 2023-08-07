/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_check_equal.hpp
 *
 */

#ifndef MORIS_FN_CHECK_EQUAL_HPP
#define MORIS_FN_CHECK_EQUAL_HPP

#include "cl_Matrix.hpp"
#include "catch.hpp"
#include "fn_equal_to.hpp"

namespace moris
{
    /**
     * Checks for moris matrices being equal using catch.
     *
     * @tparam Matrix_Type Matrix type
     * @param aMatrix1 Comparison matrix 1
     * @param aMatrix2 Comparison matrix 2
     * @param aErrorFactor Error factor for approximate floating point value comparison
     */
    template< typename Matrix_Type >
    void
    check_equal(
            Matrix< Matrix_Type > aMatrix1,
            Matrix< Matrix_Type > aMatrix2,
            real aErrorFactor = 1.0E+03 )
    {
        // Require rows and columns to be equal before checking values
        REQUIRE( aMatrix1.n_rows() == aMatrix2.n_rows() );
        REQUIRE( aMatrix1.n_cols() == aMatrix2.n_cols() );

        // Check if each entry is equal (approximately in the case of floating point numbers)
        bool tAllMatrixEntriesEqual = true;
        for ( uint iRowIndex = 0; iRowIndex < aMatrix1.n_rows(); iRowIndex++ )
        {
            for ( uint iColumnIndex = 0; iColumnIndex < aMatrix1.n_cols(); iColumnIndex++ )
            {
                if ( !equal_to( aMatrix1( iRowIndex, iColumnIndex ), aMatrix2( iRowIndex, iColumnIndex ), aErrorFactor ) )
                {
                    std::cout << "( " << iRowIndex << ", " << iColumnIndex << " ): "
                        << aMatrix1( iRowIndex, iColumnIndex ) << " =/= " << aMatrix2( iRowIndex, iColumnIndex ) << std::endl;
                    tAllMatrixEntriesEqual = false;
                }
            }
        }
        CHECK( tAllMatrixEntriesEqual );
    }
}

#endif    // MORIS_FN_CHECK_EQUAL_HPP
