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

/**
 * \def GET_FACTOR( errorfactor, ... )
 * Pretty cool trick to grab just the first variadic argument. Don't call this, it's only used for CHECK_EQUAL.
 */
#define GET_FACTOR( errorfactor, ... ) errorfactor

/**
 * \def CHECK_EQUAL( matrix1, matrix2, errorfactor = 1.0E+06 )
 * Checks for moris matrices being equal using catch, printing the file and line number of the calling test if it fails
 * Note that you must have a comma after the second matrix if you want to use the default error factor.
 * \param matrix1 Comparison matrix 1
 * \param matrix2 Comparison matrix 2
 * \param errorfactor error factor (variadic, but only first argument is used)
 */
#define CHECK_EQUAL( matrix1, matrix2, ... ) moris::check_equal( matrix1, matrix2, GET_FACTOR( __VA_ARGS__ 1.0E+06, 1.0E+06 ), #matrix1, #matrix2, __FILE__, __LINE__ )

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
            real aErrorFactor = 1.0E+06,
            std::string aMatrix1Name = "",
            std::string aMatrix2Name = "",
            std::string aFile = "",
            uint aLine = 0 )
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
                    // If this is first value check to fail, print where this function was called from
                    if ( tAllMatrixEntriesEqual )
                    {
                        std::cout << "CHECK_EQUAL() failed at " << aFile << ":" << aLine << std::endl;
                    }

                    // Print the entry location and the values of each matrix
                    std::cout << std::setprecision( 15 - ( int ) log10( aErrorFactor ) );
                    std::cout << "  ( " << iRowIndex << ", " << iColumnIndex << " ) value does not match:" << std::endl;
                    std::cout << "    " << aMatrix1( iRowIndex, iColumnIndex ) << " from " << aMatrix1Name << std::endl;
                    std::cout << "    " << aMatrix2( iRowIndex, iColumnIndex ) << " from " << aMatrix2Name << std::endl;
                    std::cout << std::setprecision( -1 );

                    // Check should fail
                    tAllMatrixEntriesEqual = false;
                }
            }
        }
        CHECK( tAllMatrixEntriesEqual );
    }
}

#endif    // MORIS_FN_CHECK_EQUAL_HPP
