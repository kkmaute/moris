/*
 * fn_eig_sym_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_EIG_SYM_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_EIG_SYM_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET1, typename ET2, typename ET3 >
    void
    eig_sym(       ET1 & aEigenValues,
                   ET2 & aEigenVectors,
             const ET3 & aA)
    {
        arma::Col< real > tTmpS; ///< arma column type required for arma::eig_sym.

        if ( ! arma::eig_sym( tTmpS, aEigenVectors, aA ) )
        {
             throw std::runtime_error ("eig_sym failed. ");
        }

        aEigenValues = tTmpS;
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_EIG_SYM_ARMA_HPP_ */
