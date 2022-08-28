/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_trans_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_TRANS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_TRANS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    namespace linalg_internal
    {
        template<typename ET >
        auto
        trans( const ET & A)
        ->decltype( arma::strans( A ) )
        {
            return arma::strans( A );
        }
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_TRANS_ARMA_HPP_ */

