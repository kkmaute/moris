/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_plus_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_PLUS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_PLUS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

// Native expressions also occur outside namespace moris; retain scalar arithmetic through ADL.
namespace Eigen
{
    template< typename E >
    auto operator+( const MatrixBase< E >& aMatrix, const typename E::Scalar& aScalar )
    {
        return ( aMatrix.array() + aScalar ).matrix();
    }

    template< typename E >
    auto operator+( const typename E::Scalar& aScalar, const MatrixBase< E >& aMatrix )
    {
        return ( aScalar + aMatrix.array() ).matrix();
    }
}    // namespace Eigen

namespace moris
{
    template< typename Matrix_Type, typename ET >
    auto
    operator+( const Eigen::MatrixBase< ET >& aA,
            const Matrix< Matrix_Type >&      aB )
            -> decltype( aA + aB.matrix_data() )
    {
        return aA + aB.matrix_data();
    }

    template< typename Matrix_Type, typename ET >
    auto
    operator+( const Matrix< Matrix_Type >& aA,
            const Eigen::MatrixBase< ET >&  aB )
            -> decltype( aA.matrix_data() + aB )
    {
        return aA.matrix_data() + aB;
    }

    template< typename Matrix_Type >
    auto operator+( const Matrix< Matrix_Type >&             aMatrix,
            const typename Matrix< Matrix_Type >::Data_Type& aScalar )
    {
        return ( aMatrix.matrix_data().array() + aScalar ).matrix();
    }

    template< typename Matrix_Type >
    auto operator+( const typename Matrix< Matrix_Type >::Data_Type& aScalar,
            const Matrix< Matrix_Type >&                             aMatrix )
    {
        return ( aScalar + aMatrix.matrix_data().array() ).matrix();
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_PLUS_EIGEN_HPP_ */
