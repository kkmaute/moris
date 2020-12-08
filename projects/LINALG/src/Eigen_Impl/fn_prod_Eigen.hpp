/*
 * fn_prod_Eigen.hpp
 *
 *  Created on: Nov 25, 2020
 *      Author: gates
 *
 *      Fixme: currently only has total prod() coded.  Eventually, need to add
 *             rowwise() and colwise() prod variations. See:
 *             https://eigen.tuxfamily.org/dox/classEigen_1_1DenseBase.html#af119d9a4efe5a15cd83c1ccdf01b3a4f
 *             https://eigen.tuxfamily.org/dox/classEigen_1_1VectorwiseOp.html#a01bcd17504f30b55b4910ddb75598f79
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
