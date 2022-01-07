/*
 * fn_isempty_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISEMPTY_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISEMPTY_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isempty( ET const & aA )
    ->decltype( aA.is_empty() )
    {
        return aA.is_empty();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISEMPTY_ARMA_HPP_ */
