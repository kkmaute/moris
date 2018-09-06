/*
 * fn_issquare_Eigen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISSQUARE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISSQUARE_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename ET >
    bool
    issquare( const Eigen::MatrixBase<ET> & aA )
    {
        if ( aA.rows() ==  aA.cols() )
        {
            return true;
        }
        return false;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISSQUARE_EIGEN_HPP_ */
