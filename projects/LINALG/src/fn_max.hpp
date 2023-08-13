/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_max.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_MAX_HPP_
#define PROJECTS_LINALG_SRC_FN_MAX_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"


#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_max_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief returns the maximum value for each column of a matrix.
     *
     * @param[in] aMat Vector.
     *
     * @return  maximum
     */
    template< typename Matrix_Type >
    auto
    max( Matrix< Matrix_Type > const & aA )
            -> decltype( max( aA.matrix_data() ) )
    {
        return max( aA.matrix_data() );
    }

    /**
     * @brief
     *
     * @tparam Matrix_Type
     * @param aA
     * @param aB
     * @return decltype( max( aA.matrix_data() ,aB.matrix_data() ) ) element-wise max of A and B
     */

    template< typename Matrix_Type_First, typename Matrix_Type_Second >
    auto
    max( Matrix< Matrix_Type_First > const &    aA,
            Matrix< Matrix_Type_Second > const & aB )
            -> decltype( max( aA.matrix_data(), aB.matrix_data() ) )
    {
        return max( aA.matrix_data(), aB.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_MAX_HPP_ */
