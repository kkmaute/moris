/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Matrix_Base_Utilities.hpp
 *
 */

#ifndef SRC_LINALG_CL_XTK_MATRIXBASE_UTILITIES_HPP_
#define SRC_LINALG_CL_XTK_MATRIXBASE_UTILITIES_HPP_

#include <iostream>

#include "cl_Matrix.hpp"
#include "fn_bubble_sort.hpp"
#include "fn_approximate.hpp"

namespace moris::xtk
{
    template< typename Matrix_Type >
    bool equal_to(
            Matrix< Matrix_Type > const &aMatrix1,
            Matrix< Matrix_Type > const &aMatrix2,
            bool                         aCheckRowsandCols = true )
    {
        bool tFlag = true;

        size_t tNumRows1 = aMatrix1.n_rows();
        size_t tNumCols1 = aMatrix1.n_cols();

        size_t tNumRows2 = aMatrix2.n_rows();
        size_t tNumCols2 = aMatrix2.n_cols();

        if ( aCheckRowsandCols )
        {
            if ( tNumRows1 != tNumRows2 )
            {
                std::cout << "There are a different number of rows in the provided matrices";
                tFlag = false;
            }

            if ( tNumCols1 != tNumCols2 )
            {
                std::cout << "There are a different number of columns in the provided matrices";
                tFlag = false;
            }
        }

        for ( size_t r = 0; r < tNumRows1; r++ )
        {
            for ( size_t c = 0; c < tNumCols1; c++ )
            {
                if ( !approximate( aMatrix2( r, c ), aMatrix1( r, c ) ) )
                {
                    tFlag = false;
                    break;
                }
            }
        }
        return tFlag;
    }

    template< typename Matrix_Type >
    bool row_equal(
            size_t const                &aRowIndex1,
            Matrix< Matrix_Type > const &aMatrix1,
            size_t const                &aRowIndex2,
            Matrix< Matrix_Type > const &aMatrix2 )
    {
        bool   tEqual    = true;
        size_t tNumCols1 = aMatrix1.n_cols();

        if ( tNumCols1 > 0 )
        {
            for ( size_t c = 0; c < tNumCols1; c++ )
            {
                if ( !approximate( aMatrix1( aRowIndex1, c ), aMatrix2( aRowIndex2, c ) ) )
                {
                    tEqual = false;
                    break;
                }
            }
        }

        return tEqual;
    }

    /*
     *
     * True if there are no duplicate faces
     * Order in the row does not matter (sorts in ascending order)
     */
    template< typename Matrix_Type >
    bool
    check_for_duplicate_rows(
            Matrix< Matrix_Type > &aMat,
            bool                   aOrderMatters = true )
    {
        bool                  tSame    = true;
        size_t                tNumCols = aMat.n_cols();
        Matrix< Matrix_Type > tRow1( 1, tNumCols );
        Matrix< Matrix_Type > tRow2( 1, tNumCols );

        // Sort the rows in ascending order if the order matters
        if ( !aOrderMatters ) { xtk::row_bubble_sort( aMat ); }

        for ( size_t i = 0; i < aMat.n_rows(); i++ )
        {
            tRow1 = aMat.get_row( i );

            for ( size_t j = i + 1; j < aMat.n_rows(); j++ )
            {
                tRow2 = aMat.get_row( j );

                tSame = xtk::equal_to( tRow1, tRow2 );

                if ( tSame )
                {
                    std::cout << "i = " << i << " j = " << j << std::endl;
                    return false;
                }
            }
        }

        return true;
    }

    template< typename Matrix_Type >
    bool
    check_for_duplicate_columns(
            Matrix< Matrix_Type > &aMat,
            bool                   aOrderMatters = true )
    {

        bool                                     tSame    = true;
        size_t                                   tNumRows = aMat.n_rows();
        std::shared_ptr< Matrix< Matrix_Type > > tCol1    = aMat.create( tNumRows, 1 );
        Matrix< Matrix_Type >                    tCol2    = aMat.create( tNumRows, 1 );

        // Sort the rows in ascending order if the order matters
        if ( !aOrderMatters ) { xtk::row_bubble_sort( aMat ); }

        for ( size_t i = 0; i < aMat.n_cols(); i++ )
        {
            aMat.get_column( i, *tCol1 );

            for ( size_t j = i + 1; j < aMat.n_cols(); j++ )
            {
                aMat.get_column( j, *tCol2 );

                tSame = xtk::equal_to( *tCol1, *tCol2 );

                if ( tSame )
                {
                    return false;
                }
            }
        }

        return true;
    }

    /*
     * Sets aMatrix.row(aRowIndex1) as aMatrix2.row(aRowIndex2)
     */
    template< typename Matrix_Type >
    void
    replace_row(
            size_t const                &aRowIndex1,
            Matrix< Matrix_Type > const &aMatrix1,
            size_t const                &aRowIndex2,
            Matrix< Matrix_Type >       &aMatrix2,
            bool                         aCheckRowsandCols = false )
    {
        size_t tNumCols1 = aMatrix1.n_cols();

        if ( aCheckRowsandCols )
        {

            size_t tNumRows1 = aMatrix1.n_rows();
            size_t tNumRows2 = aMatrix2.n_rows();
            size_t tNumCols2 = aMatrix2.n_cols();

            MORIS_ERROR( tNumRows1 = tNumRows2, "Different number of rows" );
            MORIS_ERROR( tNumCols1 = tNumCols2, "Different number of cols" );
            MORIS_ERROR( aRowIndex1 < tNumRows1, "Row index out of bounds matrix 1" );
            MORIS_ERROR( aRowIndex2 < tNumRows2, "Row index out of bounds matrix 2" );
        }

        for ( size_t i = 0; i < tNumCols1; i++ )
        {
            aMatrix2( aRowIndex2, i ) = aMatrix1( aRowIndex1, i );
        }
    }

    template< typename Matrix_Type >
    void
    fill_row(
            typename Matrix< Matrix_Type >::Data_Type const &aFillValue,
            size_t const                                    &aRowIndex,
            Matrix< Matrix_Type >                           &aMatrix1 )
    {
        size_t tNumCols1 = aMatrix1.n_cols();

        for ( size_t i = 0; i < tNumCols1; i++ )
        {
            aMatrix1( aRowIndex, i ) = aFillValue;
        }
    }

    /*
     * Fill Matrix with another matrix without resize
     *
     * Matrix 1 goes into Matrix 2
     */
    template< typename Matrix_Type >
    void conservative_copy(
            Matrix< Matrix_Type > const &aMatrix1,
            Matrix< Matrix_Type >       &aMatrix2 )
    {
        size_t tNumRow1 = aMatrix1.n_rows();
        size_t tNumCol1 = aMatrix1.n_cols();

        for ( size_t iR = 0; iR < tNumRow1; iR++ )
        {
            for ( size_t iC = 0; iC < tNumCol1; iC++ )
            {
                aMatrix2( iR, iC ) = aMatrix1( iR, iC );
            }
        }
    }

    /*
     * Similar to MATLAB operation
     * aMatToReindex(aMatRow,:)(aIndexMap)
     *
     * I.e
     * aMatToReindex = [2,3,4,5;
     *                  1,2,3,4;
     *                  4,2,1,1];
     * aMatRow = 0;
     * aIndexMap = [0,1,2,3;
     *              3,2,1,0];
     *
     * would result in
     *
     * [2,3,4,5;
     *  5,4,3,2];
     *
     */

    template< typename Matrix_Type, typename Integer_Matrix >
    Matrix< Matrix_Type >
    reindex_matrix(
            Matrix< Integer_Matrix > const &aIndexMap,
            size_t const                   &aMatRow,
            Matrix< Matrix_Type > const    &aMatToReindex )
    {
        size_t tNumRow = aIndexMap.n_rows();
        size_t tNumCol = aIndexMap.n_cols();

        Matrix< Matrix_Type > tReindexedMat( tNumRow, tNumCol );

        for ( size_t i = 0; i < tNumRow; i++ )
        {
            for ( size_t j = 0; j < tNumCol; j++ )
            {
                tReindexedMat( i, j ) = aMatToReindex( aMatRow, aIndexMap( i, j ) );
            }
        }

        return tReindexedMat;
    }
}    // namespace moris::xtk

#endif /* SRC_LINALG_CL_XTK_MATRIXBASE_UTILITIES_HPP_ */
