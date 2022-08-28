/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort_index_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SORT_INDEX_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SORT_INDEX_EIGEN_HPP_

#include <Eigen/Dense>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    moris::Matrix<DDUMat>
    sort_index( const Eigen::MatrixBase<ET> & aA )
    {
        MORIS::ERROR(false, "sort_index(mat) in Eigen not implemented yet");

        return {{0}};
    }

    template< typename ET, typename Matrix_Type, typename Num_Type>
    moris::Matrix<DDUMat>
    sort_index( const Eigen::MatrixBase<ET>    & aA,
            moris::Matrix< Matrix_Type > & aSorted,
            char const                   * aDirection,
            Num_Type                       aDimension = 0)
    {
        MORIS::ERROR(false, "sort_index(mat, dir, dim) in Eigen not implemented yet");

        return {{0}};
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SORT_INDEX_EIGEN_HPP_ */

