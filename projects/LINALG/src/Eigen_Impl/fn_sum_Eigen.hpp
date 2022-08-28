/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sum_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SUM_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SUM_EIGEN_HPP_

#include <Eigen/Dense>

namespace moris
{
    template< typename T >
    auto
    sum(
            Eigen::MatrixBase< T > const & aA )
    -> decltype( aA.sum() )
    {
        return aA.sum();
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SUM_EIGEN_HPP_ */

