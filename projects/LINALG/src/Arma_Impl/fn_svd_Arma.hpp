/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_svd_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SVD_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SVD_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename T, typename S, typename ET>
    void
    svd(
            arma::Mat<T>   & aU,
            arma::Mat<S>   & aS,
            arma::Mat<T>   & aV,
            const ET       & aM )
    {
        arma::Col< S > tS; ///< arma column type required for arma::svd.

        if ( ! arma::svd( aU, tS, aV, aM ) )
        {
            throw std::runtime_error ("SVD failed.");
        }

        aS = tS;

        // aV needs to be trimmed to get the same result as in Eigen
        aV.resize( aV.n_rows, aS.n_elem );
    }

    template<  typename S, typename ET>
    void
    svd(
            arma::Mat<S>   & aS,
            const ET       & aM )
    {
        arma::Col< S > tS; ///< arma column type required for arma::svd.

        if ( ! arma::svd( tS, aM ) )
        {
            throw std::runtime_error ("SVD failed.");
        }

        aS = tS;
    }
}

namespace moris
{
    // Keep the public output interface independent of backend storage ownership.
    template< typename U, typename S, typename V, typename E >
    void svd( Matrix< U >& aU, Matrix< S >& aS, Matrix< V >& aV, const E& aA )
    {
        moris::svd( aU.matrix_data(), aS.matrix_data(), aV.matrix_data(), aA );
    }

    template< typename S, typename E >
    void svd( Matrix< S >& aS, const E& aA )
    {
        moris::svd( aS.matrix_data(), aA );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SVD_ARMA_HPP_ */

