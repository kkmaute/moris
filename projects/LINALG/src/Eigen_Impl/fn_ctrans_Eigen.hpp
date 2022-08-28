/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_ctrans_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CTRANS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CTRANS_EIGEN_HPP_

#include <Eigen/Dense>

namespace moris
{
    template< typename ET >
    auto
    ctrans( const Eigen::MatrixBase<ET> & aA )
    ->decltype( aA.adjoint() )
    {
        return aA.adjoint();
    }

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CTRANS_EIGEN_HPP_ */

