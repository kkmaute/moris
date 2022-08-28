/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_SORT_HPP_
#define PROJECTS_LINALG_SRC_FN_SORT_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_sort_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_sort_Arma.hpp"
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
    void
    sort(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type >       & aSorted)
    {
        sort( aA.matrix_data(), aSorted );
    }

    template< typename Matrix_Type, typename Num_Type >
    void
    sort(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type >       & aSorted,
            char const                  * aDirection,
            Num_Type                      aDimension = 0)
    {
        sort( aA.matrix_data(), aSorted, aDirection, aDimension);
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_SORT_HPP_ */

