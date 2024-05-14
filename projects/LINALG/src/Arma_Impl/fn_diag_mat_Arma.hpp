/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_diag_mat_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_MAT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_MAT_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    diag_mat(
            ET&            aA,
            size_t const & ak = 0 )
            -> decltype( arma::diagmat( aA, ak ) )
    {
        return arma::diagmat( aA, ak );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_MAT_ARMA_HPP_ */
