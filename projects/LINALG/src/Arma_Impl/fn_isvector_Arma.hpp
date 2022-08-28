/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isvector_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISVECTOR_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISVECTOR_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isvector( ET const & aA )
    ->decltype( aA.is_vec() )
    {
        return aA.is_vec();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISVECTOR_ARMA_HPP_ */

