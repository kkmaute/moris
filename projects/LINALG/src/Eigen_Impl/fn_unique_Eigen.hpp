/*
 * fn_unique_Eigen.hpp
 *
 *  Created on: Sep 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_UNIQUE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_UNIQUE_EIGEN_HPP_

#include "assert.hpp"
#include "cl_Matrix.hpp"
#include "Eigen/Dense"


namespace moris
{
//-------------------------------------------------------------------------------

    template<typename ET >
    Eigen::MatrixBase<ET>
    unique( const Eigen::MatrixBase<ET> &  aMatrix )
    {
        // Eigen does not have an internal unique function
        // 1. Create copy  of input  matrix
        // 2. Sort matrix using std::sort
        // 3. Taking only the unique part and saving in ttMat
        // 4. Finding position and resize tMat

        // get number of rows from matrix implementation
        size_t n_rows = aMatrix.rows();

        // get number of cols from matrix implementation
        size_t n_cols = aMatrix.cols();

        // get length of matrix
        size_t tLength = ( n_rows < n_cols ) ? n_cols : n_rows;

        //MORIS_ASSERT( n_rows == 1 || n_cols == 1,
        //        "unique() can only be called for a "

        // create copy
        Eigen::MatrixBase<ET> tMatrix( aMatrix );

        // sort matrix
        std::sort( tMatrix.data(),  tMatrix.data()+tLength );

        /*auto tLast = std::unique(
                tMatrix.data(),
                tMatrix.data() + tLength );

        auto tPos = std::distance( tMatrix.data(), tLast );

        // check if matrix is a row matrix
        if( n_cols == 1 )
        {
            tMatrix.resize( 1, tPos );
        }
        else // matrix must be a column matrix, otherwise length() would have thrown an error
        {
            tMatrix.resize( tPos, 1 );
        }*/

        return tMatrix;
    }

//-------------------------------------------------------------------------------
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_UNIQUE_EIGEN_HPP_ */
