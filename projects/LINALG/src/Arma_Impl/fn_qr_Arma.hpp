/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_qr_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_QR_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_QR_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ETQ, typename ETR , typename ETA  >
    void
    qr(
            ETQ       & aQ,
            ETR       & aR,
            ETA const & aA )
    {
        if ( ! arma::qr( aQ, aR, aA ) )
        {
            throw std::runtime_error( "QR decomposition failed!" );
        }
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_QR_ARMA_HPP_ */

