/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linsolve.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_LINSOLVE_HPP_
#define PROJECTS_LINALG_SRC_FN_LINSOLVE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_linsolve_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_linsolve_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Solve for a linear set of equations Ax = B.
     *
     * @param[in] A The LHS Matrix
     * @param[in] B The RHS Vector
     * @param[in] aSolver Optional solver type
     *
     * @return The vector of solutions, x. Similar to B/A in Matlab
     *
     * @note If A is square, solve() is faster and more accurate than using X = inv(A)*B .
     * If A is non-square, solve() will try to provide approximate solutions to under-determined
     * as well as over-determined systems.
     * Eigen provides various options for decomposition of matrices to facilitate a linear solve.
     * We use the default option of QR decomposition with column pivoting.
     *
     */
    template< typename Matrix_Type >
    auto
    solve( Matrix< Matrix_Type > const &  aA,
            Matrix< Matrix_Type > const & aB,
            std::string const &           aSolver = "default" )
            -> decltype( solve( aA.matrix_data(), aB.matrix_data() ) )
    {
        MORIS_ASSERT( aA.n_rows() > 10, "For matrices smaller than 10x10 use inv() instead of solve().\n" );

        return solve( aA.matrix_data(), aB.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_LINSOLVE_HPP_ */
