/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_mat_to_vec.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_MAT_TO_VEC_HPP_
#define PROJECTS_LINALG_SRC_FN_MAT_TO_VEC_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#include "cl_Vector.hpp"

namespace moris
{
    /**
     * @brief Generate a flattened version of matrix to either column (default) or row vector and stores it as moris::Vector.
     *
     * for column vector: elements are copied column-wise, resulting in a column vector;
     *                    equivalent to concatenating all the columns of matrix
     *
     * for row vector:    elements are copied row-wise, resulting in a row vector;
     *                    equivalent to concatenating all the rows of matrix
     *
     * @param[in] aA            Matrix.
     * @param[in] aToRowVector  flag to create row vector
     *
     * @return  row or column vector as moris::Vector
     */

    template< typename Matrix_Type, typename Vector_Type = typename Matrix< Matrix_Type >::Data_Type >
    Vector< Vector_Type >
    mat_to_vec(
            const Matrix< Matrix_Type > &aA,
            bool                         aToRowVector = false )
    {
        // Copy directly: a fixed-size input cannot hold a differently shaped
        // flattened intermediate. Explicit indexing preserves the requested
        // ordering independently of the backend's native matrix layout.
        Vector< Vector_Type > tVector( aA.numel() );
        uint                  tIndex = 0;

        if ( aToRowVector )
        {
            for ( uint iRow = 0; iRow < aA.n_rows(); ++iRow )
            {
                for ( uint iColumn = 0; iColumn < aA.n_cols(); ++iColumn )
                {
                    tVector( tIndex++ ) = aA( iRow, iColumn );
                }
            }
        }
        else
        {
            for ( uint iColumn = 0; iColumn < aA.n_cols(); ++iColumn )
            {
                for ( uint iRow = 0; iRow < aA.n_rows(); ++iRow )
                {
                    tVector( tIndex++ ) = aA( iRow, iColumn );
                }
            }
        }

        return tVector;
    }

}    // namespace moris
#endif /* PROJECTS_LINALG_SRC_FN_MAT_TO_VEC_HPP_ */
