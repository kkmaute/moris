/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_diag_mat.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_DIAG_MAT_HPP_
#define PROJECTS_LINALG_SRC_FN_DIAG_MAT_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_diag_mat_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_diag_mat_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Generate a diagonal matrix from vector or matrix
     *
     * if aA is a vector, a square matrix with the k-th diagonal containing a copy of the vector is generated;
     * all other elements are set to zero
     *
     * if aA is a matrix, a matrix with the k-th diagonal containing a copy of the k-th diagonal of X is generated;
     * all other elements are set to zero.
     *
     * @param[in] aA Matrix.
     * @param[in] ak Diagonal index.
     * The argument k is optional; by default the main diagonal
     * is extracted (k=0).\n
     * For k > 0, the k-th super-diagonal is generated
     * (top-right corner).\n
     * For k < 0, the k-th sub-diagonal is generated
     * (bottom-left corner).\n
     */

    template< typename Matrix_Type >
    auto
    diag_mat( Matrix< Matrix_Type > &aA,
            size_t const            &ak = 0 )
            -> decltype( diag_mat( aA.matrix_data(), ak ) )
    {
        return diag_mat( aA.matrix_data(), ak );
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_DIAG_MAT_HPP_ */
