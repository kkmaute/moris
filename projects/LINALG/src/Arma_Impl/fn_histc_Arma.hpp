/*
 * fn_histc_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_HISTC_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_HISTC_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET1, typename ET2>
    auto
    histc( ET1 const & aA,
           ET2 const & aB )
    ->decltype( arma::histc( aA, aB ) )
    {
        return ( arma::histc( aA, aB ) ) ;
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_HISTC_ARMA_HPP_ */
