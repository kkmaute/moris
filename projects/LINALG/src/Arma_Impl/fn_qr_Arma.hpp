/*
 * fn_qr_Arma.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: messe
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
