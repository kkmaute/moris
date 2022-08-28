/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_norm.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_NORM_HPP_
#define PROJECTS_LINALG_SRC_FN_NORM_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_norm_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_norm_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Calculate the L2 norm of a matrix.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename Matrix_Type >
    auto
    norm( const Matrix< Matrix_Type > & aA )
        -> decltype( norm( aA.matrix_data() ) )
    {
        return norm( aA.matrix_data() );
    }

    template< typename Matrix_Type >
    auto
    norm( Matrix< Matrix_Type > & aA )
        -> decltype( norm( aA.matrix_data() ) )
    {
        return norm( aA.matrix_data() );
    }

}

#endif /* PROJECTS_LINALG_SRC_FN_NORM_HPP_ */

