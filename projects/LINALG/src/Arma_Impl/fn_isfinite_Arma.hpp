/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isfinite_Arma.hpp
 *
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

