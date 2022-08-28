/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isvector_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISVECTOREIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISVECTOR_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename ET >
    bool
    isvector( const Eigen::MatrixBase<ET> & aA )
    {
        if ( aA.rows() == 1 || aA.cols() == 1 )
        {
            return true;
        }

        return false;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISVECTOR_EIGEN_HPP_ */

