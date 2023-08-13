/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_dot.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_DOT_HPP_
#define PROJECTS_LINALG_SRC_FN_DOT_HPP_
// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_dot_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_dot_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Calculates the dot (scalar) product between two matrices or vectors of the same size
     * @f$ s=\mathbf{A}_{i,j}*\mathbf{B}_{i,j}@f$ \n
     *
     *@param[in] aA A given matrix
     *@param[in] aB A given matrix
     *
     * Example:
     * @include LNA/src/fn_dot.inc
     *
     */
    template< typename Matrix_Type >
    auto
    dot(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type > const & aB)
    -> decltype( dot( aA.matrix_data(), aB.matrix_data() ) )
    {
        return dot( aA.matrix_data(), aB.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_DOT_HPP_ */

