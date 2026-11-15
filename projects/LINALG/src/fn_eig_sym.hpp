/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eig_sym.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_EIG_SYM_HPP_
#define PROJECTS_LINALG_SRC_FN_EIG_SYM_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_eig_sym_Eigen.hpp"
#elif defined( MORIS_USE_ARMA )
#include "fn_eig_sym_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Eigenvalue decomposition of dense symmetric/Hermitian Matrix A.
     *
     * @param[in] aA Matrix.
     *
     * @param[out] eigenvalues The eigenvalues of Matrix aA.
     *
     * @param[out] eigenvectors The eigenvectors of Matrix aA.
     */
    template< typename Values, typename Vectors, typename Input >
    void eig_sym( Matrix< Values > &aEigenvalues,
            Matrix< Vectors >      &aEigenvectors,
            const Matrix< Input >  &aA )
    {
        eig_sym( aEigenvalues, aEigenvectors, aA.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_EIG_SYM_HPP_ */
