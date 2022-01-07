/*
 * fn_svd_Arma.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
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

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SVD_ARMA_HPP_ */
