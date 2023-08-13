/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cond.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_COND_HPP_
#define PROJECTS_LINALG_SRC_FN_COND_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_cond_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_cond_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Calculates the two-norm condition number of a matrix
     *
     * @return Calculates the two-norm condition number of a matrix (ratio
     * of the largest singular value to the smallest one).\n
     * Large condition numbers suggest that matrix A is nearly singular.
     *
     * @note Neither Armadillo nor Eigen compute the p-norm condition number.\n
     * Moreover, a one-norm condition number estimate function is not provided
     * either. If there is a need for a different condition number than a two-
     * norm, or a condition number estimate, it needs to be implemented.
     */
    template< typename Matrix_Type >
    auto
    cond(Matrix< Matrix_Type > const & aA)
    -> decltype( cond( aA.matrix_data()) )
    {
        return cond( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_COND_HPP_ */

