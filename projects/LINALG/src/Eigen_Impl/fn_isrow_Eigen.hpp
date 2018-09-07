/*
 * fn_isrow_Eigen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISROW_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISROW_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename ET >
    bool
    isrow( const Eigen::MatrixBase<ET> & aA )
    {
         if ( aA.rows() == 1 && aA.cols() >= 1 )
         {
             return true;
         }

         return false;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISROW_EIGEN_HPP_ */
