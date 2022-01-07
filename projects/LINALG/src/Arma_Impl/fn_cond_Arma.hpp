/*
 * fn_cond_Arma.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COND_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COND_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET>
    auto
    cond(ET const & aA)
    -> decltype( arma::cond( aA ) )
    {
        return arma::cond( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COND_ARMA_HPP_ */
