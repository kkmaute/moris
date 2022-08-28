/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_svd.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_SVD_HPP_
#define PROJECTS_LINALG_SRC_FN_SVD_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_svd_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_svd_Arma.hpp"
#endif

namespace moris
{
    template< typename M,
              typename N >
    void
    svd(
                  Matrix< M > & aU,
                  Matrix< N > & aS,
                  Matrix< M > & aV,
            const Matrix< M > & aM )
    {
        svd( aU.matrix_data(),
             aS.matrix_data(),
             aV.matrix_data(),
             aM.matrix_data() );
    }

    template< typename M, typename N >
    void
    svd(          Matrix< N > & aS,
            const Matrix< M > & aM )
    {
        svd( aS.matrix_data(),
             aM.matrix_data() );
    }

}

#endif /* PROJECTS_LINALG_SRC_FN_SVD_HPP_ */

