/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_join_horiz.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_JOIN_VERT_HPP_
#define PROJECTS_LINALG_SRC_FN_JOIN_VERT_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_join_horiz_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_join_horiz_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief vertical concatenation
     *
     * @param[in] aA The First Matrix
     * @param[in] aB The Second Matrix
     *
     * @return join the corresponding columns of the given matrices
     *
     * @note The given matrices must have the same number of rows
     */
    template< typename Matrix_Type >
    auto
    join_vert(
            Matrix< Matrix_Type > const &aA,
            Matrix< Matrix_Type > const &aB )
            -> decltype( join_vert( aA.matrix_data(), aB.matrix_data() ) )
    {
        MORIS_ASSERT( aA.n_cols() == 0 || aB.n_cols() == 0 || aA.n_cols() == aB.n_cols(),
                "moris::join_horiz - The number of rows in the two matrices must be equal." );

        return join_vert( aA.matrix_data(), aB.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_JOIN_VERT_HPP_ */
