/*
 * fn_sqrtmat_Arma.hpp
 *
 *  Created on: Feb 01, 2020
 *      Author: maute
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SQRTMAT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SQRTMAT_ARMA_HPP_

#include "cl_Matrix.hpp"

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

    template< typename ET >
    ET
    sqrtmat( ET const & aA )
    {
        // check that sqrt of matrix is real
        MORIS_ASSERT( arma::norm( arma::imag( arma::sqrtmat( aA ))) < 1e-12*arma::norm( aA ),
                "square root of matrix is not real because it has negative eigen values; this case is not implemented\n");

        return (ET) arma::real( arma::sqrtmat( aA ) );
    }
    /**
     * @brief Computes the square root of a matrix such that A = X*X.
     *
     * @param[in] A matrix
     *
     * @return square root of A
     *
     */

    template< typename ET >
    void
    sqrtmat(
            ET const & aA,
            ET       & aX )
    {
        // check that sqrt of matrix is real
        MORIS_ASSERT( arma::norm( arma::imag( arma::sqrtmat( aA ))) < 1e-12*arma::norm( aA ),
                "square root of matrix is not real because it has negative eigen values; this case is not implemented\n");

        aX = arma::real( arma::sqrtmat( aA ) );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SQRTMAT_ARMA_HPP_ */
