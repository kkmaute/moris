/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort_index.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_SORT_INDEX_HPP_
#define PROJECTS_LINALG_SRC_FN_SORT_INDEX_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_sort_index_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_sort_index_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Sorts colums of a vector/matrix.
     *
     * @param[in] aMat Vector.
     *
     * @return  sorted vector
     */
    template< typename Matrix_Type >
    Matrix<DDUMat>
    sort_index( Matrix< Matrix_Type > const & aA )
    {
        MORIS_ASSERT( aA.n_cols() == 1,
                "sort_index - only works with column vector.");

        return sort_index( aA.matrix_data() );
    }

    template< typename Matrix_Type >
    Matrix<DDUMat>
    sort_index(
            Matrix< Matrix_Type > const & aA,
            char const                  * aDirection)
    {
        MORIS_ASSERT( aA.n_cols() == 1,
                "sort_index - only works with column vector.");

        return sort_index( aA.matrix_data(), aDirection );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_SORT_INDEX_HPP_ */

