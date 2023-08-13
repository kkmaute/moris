/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_intersect.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_INTERSECT_HPP_
#define PROJECTS_LINALG_SRC_FN_INTERSECT_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_intersect_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_intersect_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Return the unique elements common to both A and B, sorted in ascending order
     *
     * @param[in] aAmat     Vector or matrix.
     * @param[in] aBmat     Vector or matrix.
     *
     * @return  Column vector containing the indices of comment elements
     */
    template< typename Matrix_Type >
    auto
    intersect(
            const Matrix< Matrix_Type >& aAmat,
            const Matrix< Matrix_Type >& aBmat )
            -> decltype( intersect( aAmat.matrix_data(), aBmat.matrix_data() ) )
    {
        return intersect( aAmat.matrix_data(), aBmat.matrix_data() );
    }

    /**
     * @brief Determine the unique elements common to both A and B, sorted in ascending order
     *        and identifies the indices Aind and Bind in A and B of the comment elements, such that
     *        A(Aind) = B(Bind)
     *
     * @param[in]  aAmat     Vector or matrix.
     * @param[in]  aBmat     Vector or matrix.
     * @param[out] aCmat     Vector or matrix.
     * @param[out] aAind     Vector or matrix.
     * @param[out] aBind     Vector or matrix.
     */

    template< typename Matrix_Type >
    void
    intersect(
            const Matrix< Matrix_Type >& aAmat,
            const Matrix< Matrix_Type >& aBmat,
            Matrix< Matrix_Type >&       aCmat,
            Matrix< DDUMat >&            aAind,
            Matrix< DDUMat >&            aBind )
    {
        intersect(
                aAmat.matrix_data(),
                aBmat.matrix_data(),
                aCmat.matrix_data(),
                aAind.matrix_data(),
                aBind.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_UNIQUE_HPP_ */
