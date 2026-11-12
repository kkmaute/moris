/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sort_index_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_INDEX_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_INDEX_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    Matrix< DDUMat > sort_index( const ET& aA, const char* aDirection = "ascend" )
    {
        const arma::uvec tSorted = arma::sort_index( aA, aDirection );
        Matrix< DDUMat > tResult( tSorted.n_rows, 1 );

        // Armadillo's uword can be wider than MORIS uint; copy with conversion
        // for both default and explicit sort directions.
        for ( uint i = 0; i < tResult.n_rows(); ++i )
        {
            tResult( i ) = tSorted( i );
        }
        return tResult;
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_INDEX_ARMA_HPP_ */
