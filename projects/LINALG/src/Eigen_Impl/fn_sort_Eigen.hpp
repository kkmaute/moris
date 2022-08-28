/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SORT_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SORT_EIGEN_HPP_

#include <Eigen/Dense>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET, typename Matrix_Type >
    void
    sort( const Eigen::MatrixBase<ET> & aA,
            moris::Matrix< Matrix_Type > & aSorted )
    {
        typedef typename Matrix< Matrix_Type >::Data_Type Type;

        //MORIS::ERROR(false, "sort in eigen not implemented yet");
        // Eigen does not have an internal sort function
        // 1. Create copy  of input  matrix
        // 2. Sort matrix using std::sort
        aSorted.matrix_data() = aA;

        // make sure that matrix is row or col matrix
        MORIS_ASSERT( aSorted.n_rows() == 1 || aSorted.n_cols() == 1,
                "sort() can only be performed on a vector" );

        Type * tData = aSorted.data();
        moris::uint tLen = aSorted.numel();
        std::sort(tData,tData+tLen);
    }

    template< typename ET, typename Matrix_Type, typename Num_Type>
    void
    sort( const Eigen::MatrixBase<ET>    & aA,
            moris::Matrix< Matrix_Type > & aSorted,
            char const                   * aDirection,
            Num_Type                       aDimension = 0)
    {
        MORIS::ERROR(false, "sort(mat, matSorted, dir, dim) in Eigen not implemented yet");
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SORT_EIGEN_HPP_ */

