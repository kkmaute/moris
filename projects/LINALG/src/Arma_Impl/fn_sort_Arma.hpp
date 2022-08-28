/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET, typename Matrix_Type >
    void
    sort(
            ET const                     & aA,
            moris::Matrix< Matrix_Type > & aSorted )
    {
        aSorted = arma::sort( aA );
    }

    template< typename ET, typename Matrix_Type, typename Num_Type >
    void
    sort(
            ET const                     & aA,
            moris::Matrix< Matrix_Type > & aSorted,
            char const                   * aDirection,
            Num_Type                       aDimension = 0 )
    {
        aSorted = arma::sort( aA, aDirection, aDimension);
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_ARMA_HPP_ */

