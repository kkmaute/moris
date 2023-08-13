/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_min.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_MIN_HPP_
#define PROJECTS_LINALG_SRC_FN_MIN_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_min_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief returns the minimum value for each column of a matrix.
     *
     * @param[in] aMat Vector.
     *
     * @return  minimum
     */
    template< typename Matrix_Type >
    auto
    min( Matrix< Matrix_Type > const & aA )
            -> decltype( min( aA.matrix_data() ) )
    {
        return min( aA.matrix_data() );
    }


    /**
     * @brief
     *
     * @tparam Matrix_Type
     * @param aA
     * @param aB
     * @return decltype( min( aA.matrix_data() ,aB.matrix_data() ) ) element-wise min of A and B
     */

    template< typename Matrix_Type >
    auto
    min( Matrix< Matrix_Type > const &    aA,
            Matrix< Matrix_Type > const & aB )
            -> decltype( min( aA.matrix_data(), aB.matrix_data() ) )
    {
        return min( aA.matrix_data(), aB.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_MIN_HPP_ */
