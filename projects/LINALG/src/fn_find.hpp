/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_find.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_FIND_HPP_
#define PROJECTS_LINALG_SRC_FN_FIND_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_find_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_find_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Find positions of non zeros in a vector.
     *
     * @param[in] aMat     Vector.
     *
     * @return  Vector of non zero positions found
     */
    template< typename Matrix_Type >
    auto
    find( const Matrix< Matrix_Type >& aA )
            -> decltype( find( aA.matrix_data() ) )
    {
        return find( aA.matrix_data() );
    }

    /**
     * @brief Find the first aB positions of non zeros in a vector.
     *
     * @param[in] aMat  Vector.
     * @param[in] aB    Specific number of indices.
     *
     * @return Vector of non zero positions found
     */
    template< typename Matrix_Type >
    auto
    find(
            const Matrix< Matrix_Type >& aA,
            const moris::uint&           aB )
            -> decltype( find( aA.matrix_data(), aB ) )
    {

        return find( aA.matrix_data(), aB );
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_FIND_HPP_ */
