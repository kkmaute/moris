/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_elemwise_mult_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_MULT_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_MULT_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
    //--------------------------------------------------------------------------------

    template< typename Matrix_Type >
    auto
    operator%(
            const Matrix< Matrix_Type > &aA,
            const Matrix< Matrix_Type > &aB )
            -> decltype( aA.matrix_data().cwiseProduct( aB.matrix_data() ) )
    {
        return aA.matrix_data().cwiseProduct( aB.matrix_data() );
    }

    //--------------------------------------------------------------------------------

    template< typename Matrix_Type, typename ET >
    auto
    operator%( const Matrix< Matrix_Type > &aA,
            const Eigen::MatrixBase< ET >  &aB )
            -> decltype( aA.matrix_data().cwiseProduct( aB ) )
    {
        return aA.matrix_data().cwiseProduct( aB );
    }

    //--------------------------------------------------------------------------------

    template< typename Matrix_Type, typename ET >
    auto
    operator%( const Eigen::MatrixBase< ET > &aA,
            const Matrix< Matrix_Type >      &aB )
            -> decltype( aA.cwiseProduct( aB.matrix_data() ) )
    {
        return aA.cwiseProduct( aB.matrix_data() );
    }

    //--------------------------------------------------------------------------------


    //--------------------------------------------------------------------------------
}    // namespace moris


// Native expression operands need an ADL-visible operator after wrapper unwrapping.
namespace Eigen
{
    template< typename ETA, typename ETB >
    auto
    operator%( const Eigen::MatrixBase< ETA > &aA,
            const Eigen::MatrixBase< ETB >    &aB )
            -> decltype( aA.cwiseProduct( aB ) )
    {
        return aA.cwiseProduct( aB );
    }

}    // namespace Eigen
namespace moris
{
    using Eigen::operator%;
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_MULT_EIGEN_HPP_ */
