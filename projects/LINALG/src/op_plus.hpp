/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_plus.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_OP_PLUS_HPP_
#define PROJECTS_LINALG_SRC_OP_PLUS_HPP_

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_plus_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_plus_Arma.hpp"
#endif

#include "cl_Matrix.hpp"

namespace moris
{
    /**
     * @brief Addition of two arrays
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ij} + \mathbf{B}_{ij} @f$
     * where C is sum of A and B. If A is an m-by-p and B is a m-by-p matrix, then C is an m-by-p matrix.
     *
     * Example:
     * @include LNA/src/op_plus.inc
     *
     */
    template< typename Matrix_Type >
    auto
    operator+(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type > const & aB )
    ->decltype( aA.matrix_data() + aB.matrix_data() )
    {
        return aA.matrix_data() + aB.matrix_data();
    }

    template< typename Matrix_Type >
    auto
    operator+( Matrix< Matrix_Type > const & aA )
    ->decltype( aA.matrix_data() )
    {
        return aA.matrix_data();
    }
}

#endif /* PROJECTS_LINALG_SRC_OP_PLUS_HPP_ */

