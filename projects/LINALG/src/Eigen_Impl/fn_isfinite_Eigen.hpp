/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isfinite_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISFINITE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISFINITE_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename ET >
    auto
    isfinite( const Eigen::MatrixBase<ET> & aA )
    ->decltype( aA.allFinite() )
    {
        return aA.allFinite();
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISFINITE_EIGEN_HPP_ */

