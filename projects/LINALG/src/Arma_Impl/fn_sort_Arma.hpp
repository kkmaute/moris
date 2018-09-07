/*
 * fn_sort_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_ARMA_HPP_
#include <armadillo>

namespace moris
{
    template< typename ET, typename Type, typename Matrix_Type >
    void
    sort( ET const                                    & aA,
                   moris::Matrix< Type, Matrix_Type > & aSorted )
    {
        aSorted = sort( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SORT_ARMA_HPP_ */
