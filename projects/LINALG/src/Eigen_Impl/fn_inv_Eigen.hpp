/*
 * fn_inv_Eigen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_INV_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_INV_EIGEN_HPP_

#include <Eigen/Dense>

namespace moris
{
    template< typename ET >
    auto
    inv( const Eigen::MatrixBase<ET> & aA )
    ->decltype( (aA).inverse() )
    {
        return (aA).inverse();
    }


}




#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_INV_EIGEN_HPP_ */
