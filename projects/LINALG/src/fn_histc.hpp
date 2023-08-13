/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_histc.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_HISTC_HPP_
#define PROJECTS_LINALG_SRC_FN_HISTC_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_histc_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_histc_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief find vector.
     *
     * @param[in] aMat     Vector.
     * @param[in] WhichValue    Value, which is in the vector.
     *
     * @return  Vector of items found
     */
    template< typename Matrix_Type >
    auto
    histc( Matrix< Matrix_Type > const & aA,
           Matrix< Matrix_Type > const & aB )
    -> decltype( moris::histc( aA.matrix_data(), aB.matrix_data() ) )
    {
        return moris::histc( aA.matrix_data(), aB.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_HISTC_HPP_ */

