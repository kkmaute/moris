/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sylvester.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_SYLVESTER_HPP_
#define PROJECTS_LINALG_SRC_FN_SYLVESTER_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_sylvester_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_sylvester_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Solve the Sylvester equation, i.e. A*X + X*B + C = 0, where X is unknown.
     *
     * @param[in] A matrix
     * @param[in] B matrix
     * @param[in] C matrix
     *
     * @return solution of system
     *
     */
    template< typename Matrix_Type >
    auto
    sylvester(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type > const & aB,
            Matrix< Matrix_Type > const & aC )
    -> decltype( sylvester( aA.matrix_data(), aB.matrix_data(), aC.matrix_data() ) )
    {
        return sylvester(aA.matrix_data(), aB.matrix_data(), aC.matrix_data() );
    }

    /**
     * @brief Solve the Sylvester equation, i.e. A*X + X*B + C = 0, where X is unknown.
     *
     * @param[in] A matrix
     * @param[in] B matrix
     * @param[in] C matrix
     * @param[out] X matrix
     *
     */

    template< typename Matrix_Type >
    void
    sylvester(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type > const & aB,
            Matrix< Matrix_Type > const & aC,
            Matrix< Matrix_Type >       & aX)
   {
        sylvester( aA.matrix_data(), aB.matrix_data(), aC.matrix_data(), aX.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_SYLVESTER_HPP_ */

