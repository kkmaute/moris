/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_reshape_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_RESHAPE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_RESHAPE_EIGEN_HPP_

#include <Eigen/Dense>
#include "cl_Matrix.hpp"

namespace moris
{

template<  typename TYPE, typename INT >
Eigen::Map<Eigen::Matrix<TYPE,   Eigen::Dynamic, Eigen::Dynamic>>
reshape(
        Eigen::Matrix<TYPE,   Eigen::Dynamic, Eigen::Dynamic>  & aA,
        INT                aB,
        INT                aC)
        {
    Eigen::Map<Eigen::Matrix<TYPE,   Eigen::Dynamic, Eigen::Dynamic>> tM2(aA.data(), aB, aC);
    return tM2;
        }

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_RESHAPE_EIGEN_HPP_ */

