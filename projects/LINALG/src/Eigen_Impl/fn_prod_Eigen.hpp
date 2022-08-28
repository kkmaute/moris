/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_prod_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_PROD_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_PROD_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include <Eigen/Dense>

namespace moris
{

    template<typename ET >
    auto
    prod( const Eigen::MatrixBase<ET> &  aA)
    ->decltype( aA.prod() )
    {
        return  aA.prod();
    }

    template<typename ET >
    auto
    prod( Eigen::MatrixBase<ET> &  aA)
    ->decltype( aA.prod() )
    {
        return  aA.prod();
    }

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_PROD_EIGEN_HPP_ */

