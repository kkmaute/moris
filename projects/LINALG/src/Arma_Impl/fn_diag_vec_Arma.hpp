/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_diag_vec_Arma.hpp
 *
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

