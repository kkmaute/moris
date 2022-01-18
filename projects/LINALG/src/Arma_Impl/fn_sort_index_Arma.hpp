/*
 * fn_sort_index_Arma.hpp
 *
 *  Created on: Jan 16, 2022
 *      Author: Kurt Maute
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_INDEX_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_INDEX_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    Matrix<DDUMat>
    sort_index(ET const & aA )
    {
        arma::uvec tSorted = arma::sort_index( aA );

        return tSorted;
    }

    template< typename ET >
    Matrix<DDUMat>
    sort_index(
            ET const                     & aA,
            char const                   * aDirection )

    {
        arma::uvec tSorted = arma::sort_index( aA, aDirection );

        Matrix<DDUMat> tMat(tSorted.n_rows,1);

        // FIXME: Find better method to copy uvec into Matrix<DDUMat>
        //        Issue arma:uvec is column vector of long long uint
        for (uint i=0; i<tMat.n_rows(); ++i)
        {
            tMat(i) = tSorted(i);
        }

        return tMat;
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_INDEX_ARMA_HPP_ */
