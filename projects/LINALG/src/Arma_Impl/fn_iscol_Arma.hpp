/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_iscol_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCOL_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCOL_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    iscol( ET const & aA )
    ->decltype( aA.is_colvec() )
    {
        return aA.is_colvec();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCOL_ARMA_HPP_ */

