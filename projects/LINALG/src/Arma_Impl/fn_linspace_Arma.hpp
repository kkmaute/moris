/*
 * fn_linspace_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSPACE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSPACE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename T >
    auto
    linspace_base( T             const & aStart,
                   T             const & aEnd,
                   moris::size_t const & aN )
    ->decltype( arma::linspace< arma::Mat< T > >( aStart, aEnd, aN ) )
    {
        return arma::linspace< arma::Mat< T > >( aStart, aEnd, aN );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSPACE_ARMA_HPP_ */
