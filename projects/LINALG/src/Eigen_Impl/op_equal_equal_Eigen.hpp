/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_equal_equal_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_EQUAL_EQUAL_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_EQUAL_EQUAL_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template< typename Matrix_Type, typename ET >
auto
operator==( const Eigen::MatrixBase<ET> &  aA,
            const Matrix< Matrix_Type > & aB )
->decltype( aA.array() == aB.matrix_data().array() )
{
    return  aA.array() == aB.matrix_data().array();
}

template< typename Matrix_Type, typename ET >
auto
operator==(const Matrix< Matrix_Type > & aA,
           const Eigen::MatrixBase<ET> &  aB)
->decltype( aA.matrix_data().array() == aB.array() )
{
    return  aA.matrix_data().array() == aB.array();
}

template< typename ET1, typename ET2 >
auto
operator==( const Eigen::MatrixBase<ET1> &  aA,
            const Eigen::MatrixBase<ET2> &  aB)
->decltype( aA.array() == aB.array() )
{
    return  aA.array() == aB.array();
}

// scalar to matrix comparisons
// not  typename Eigen::DenseBase< ET1 >::Scalar is way to acces underlying numeric type in eigen
template< typename ET1 >
auto
operator==(Eigen::MatrixBase<ET1> const & aA,
           typename Eigen::DenseBase< ET1 >::Scalar const & aB )
->decltype( aA.array() == aB  )
{
    return aA.array() ==  aB ;
}

template< typename ET1 >
auto
operator==(typename Eigen::DenseBase< ET1 >::Scalar & aA,
           Eigen::MatrixBase<ET1> const & aB)
->decltype( aA == aB.array() )
{
    return aA == aB.array();
}
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_EQUAL_EQUAL_EIGEN_HPP_ */

