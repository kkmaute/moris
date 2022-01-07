/*
 * fn_isfinite_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCFINITE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISFINITE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isfinite( ET const & aA )
    ->decltype( aA.is_finite() )
    {
        return aA.is_finite();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISFINITE_ARMA_HPP_ */
