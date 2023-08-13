/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_minus.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_OP_MINUS_HPP_
#define PROJECTS_LINALG_SRC_OP_MINUS_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_minus_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_minus_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Subtraction of two arrays
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ij} - \mathbf{B}_{ij} @f$
     * where C is suntraction of B from A. If A is an m-by-p and B is a m-by-p matrix, then C is an m-by-p matrix.
     *
     * Example:
     * @include LNA/src/op_minus.inc
     *
     */
    template< typename Matrix_Type >
    auto
    operator-(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type > const & aB )
    ->decltype( aA.matrix_data() - aB.matrix_data() )
    {
        return aA.matrix_data() - aB.matrix_data();
    }

    template< typename Matrix_Type>
    auto
    operator-(
            typename Matrix< Matrix_Type >::Data_Type const & aScalar,
                     Matrix< Matrix_Type > const & aMatrix )
    ->decltype( scalar_subtraction(aScalar,aMatrix.matrix_data()) )
    {
        return scalar_subtraction(aScalar,aMatrix.matrix_data());
    }

    template< typename Matrix_Type >
    auto
    operator-(
            Matrix< Matrix_Type > const & aMatrix,
            typename Matrix< Matrix_Type >::Data_Type const & aScalar
    )
    ->decltype( scalar_subtraction(aMatrix.matrix_data(), aScalar) )
    {
        return scalar_subtraction(aMatrix.matrix_data(), aScalar);
    }

    template< typename Matrix_Type >
    auto
    operator-( Matrix< Matrix_Type > const & aA )
    ->decltype( - aA.matrix_data() )
    {
        return - aA.matrix_data();
    }
}

#endif /* PROJECTS_LINALG_SRC_OP_MINUS_HPP_ */

