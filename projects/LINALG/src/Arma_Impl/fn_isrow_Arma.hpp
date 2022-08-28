/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isrow_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISROW_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISROW_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isrow( ET const & aA )
    ->decltype( aA.is_rowvec() )
    {
        return aA.is_rowvec();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISROW_ARMA_HPP_ */

