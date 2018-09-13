/*
 * fn_svd.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_FN_SVD_HPP_
#define PROJECTS_LINALG_SRC_FN_SVD_HPP_


// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_svd_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_svd_Arma.hpp"
#endif

namespace moris
{
    template< typename  T, typename S, typename M, typename N >
    void
    svd(
                  Matrix< T, M > & aU,
                  Matrix< S, N > & aS,
                  Matrix< T, M > & aV,
            const Matrix< T, M > & aM )
    {
        svd( aU.matrix_data(),
             aS.matrix_data(),
             aV.matrix_data(),
             aM.matrix_data() );
    }

    template< typename  T, typename S, typename M, typename N >
    void
    svd(          Matrix< S, N > & aS,
            const Matrix< T, M > & aM )
    {
        svd( aS.matrix_data(),
             aM.matrix_data() );
    }

}


#endif /* PROJECTS_LINALG_SRC_FN_SVD_HPP_ */
