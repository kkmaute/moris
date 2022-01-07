/*
 * fn_issquare_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISSQUARE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISSQUARE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    issquare( ET const & aA )
    ->decltype( aA.is_square() )
    {
        return aA.is_square();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISSQUARE_ARMA_HPP_ */
