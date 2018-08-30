/*
 * fn_comp_abs_Eigen.hpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_COMP_ABS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_COMP_ABS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template<typename ET>
auto comp_abs(Eigen::MatrixBase<ET> & aExpTemplate)
->decltype(aExpTemplate.cwiseAbs())
{
    return aExpTemplate.cwiseAbs();
}

template<typename ET>
auto comp_abs(Eigen::MatrixBase<ET> const & aExpTemplate)
->decltype(aExpTemplate.cwiseAbs())
{
    return aExpTemplate.cwiseAbs();
}

}



#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_COMP_ABS_EIGEN_HPP_ */
