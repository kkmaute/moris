/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eye_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_EYE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_EYE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET>
    inline
    void
    eye(
            moris::size_t const   & aNumRows,
            moris::size_t const   & aNumCols,
            ET                    & aEyeMat )
    {
        aEyeMat =arma::eye( aNumRows, aNumCols );
    }

    namespace linalg_internal
    {
        inline
        auto
        eye(
                moris::size_t const   & aNumRows,
                moris::size_t const   & aNumCols )
        ->decltype( arma::eye( aNumRows, aNumCols ) )
        {
            return arma::eye( aNumRows, aNumCols );
        }
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_EYE_ARMA_HPP_ */

