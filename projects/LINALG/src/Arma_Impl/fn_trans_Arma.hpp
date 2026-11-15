/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_trans_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_TRANS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_TRANS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris::linalg_internal
{
    template< typename ET >
    auto
    trans( const ET &A )
            -> decltype( arma::strans( A ) )
    {
        return arma::strans( A );
    }
}    // namespace moris::linalg_internal

namespace moris
{
    // Preserve the native lazy expression; the transpose borrows its operands.
    template< typename Type, typename Expression >
    auto
    trans( const arma::Base< Type, Expression > &aExpression )
            -> decltype( linalg_internal::trans( aExpression.get_ref() ) )
    {
        return linalg_internal::trans( aExpression.get_ref() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_TRANS_ARMA_HPP_ */
