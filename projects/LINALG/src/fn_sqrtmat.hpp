/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sqrtmat.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_SQRTMAT_HPP_
#define PROJECTS_LINALG_SRC_FN_SQRTMAT_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_sqrtmat_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_sqrtmat_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Computes the square root of a matrix such that A = X*X.
     *
     * @param[in] A matrix
     *
     * @return square root of A
     *
     */
    template< typename Matrix_Type >
    auto
    sqrtmat( Matrix< Matrix_Type > const & aA )
    -> decltype( sqrtmat( aA.matrix_data() ) )
    {
        return sqrtmat( aA.matrix_data() );
    }

    /**
     * @brief Computes the square root of a matrix such that A = X*X.
     *
     * @param[in] A matrix
     * @param[out] X matrix
     *
     * @return square root of A
     *
     */
    template< typename Matrix_Type >
    void
    sqrtmat(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type >       & aX )
    {
        return sqrtmat( aA.matrix_data(), aX.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_SQRTMAT_HPP_ */

