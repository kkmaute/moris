/*
 * fn_diag_vec_Arma.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_VEC_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_VEC_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET>
    auto
    diag_vec(
            ET & aA,
            size_t     const & ak = 0 )
    ->decltype(arma::diagvec( aA, ak ))
    {
        return arma::diagvec( aA, ak );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_VEC_ARMA_HPP_ */
