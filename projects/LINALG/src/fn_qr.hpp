/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_qr.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_QR_HPP_
#define PROJECTS_LINALG_SRC_FN_QR_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_qr_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_qr_Arma.hpp"
#endif

namespace moris
{
     /* @brief QR decomposition of dense matrix
     *
     * @param[out] aQ Q-matrix
     * @param[out] aR R-matrix
     * @param[in]  aA Input matrix
     *
     * The QR decomposition of a matrix aA is such that \f$ A = Q*R \f$ where
     * \f$ Q^T*Q = I \f$ and R is an upper triangular non-square matrix. If A
     * is m-by-n, then Q is m-by-m and R is m-by-n.
     *
     * Example:
     * @include LNA/src/fn_qr.inc
     *
     */
    template< typename Matrix_Type >
    void
    qr(
            Matrix< Matrix_Type >       & aQ,
            Matrix< Matrix_Type >       & aR,
            Matrix< Matrix_Type > const & aA )
    {
        qr( aQ.matrix_data(),
            aR.matrix_data(),
            aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_QR_HPP_ */

