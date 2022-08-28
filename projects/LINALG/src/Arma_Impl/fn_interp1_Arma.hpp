/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_interp1_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INTERP1_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INTERP1_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    interp1(ET const &  aX,
            ET const &  aY,
            ET const &  aXI,
            ET & aYI )
    ->decltype( arma::interp1(aX,aY,aXI,aYI) )
    {
        return arma::interp1( aX, aY, aXI, aYI );
    }

    template< typename Matrix_Type, typename ET >
    auto
    interp1(ET                     const &  aX,
            ET                     const &  aY,
            Matrix< Matrix_Type >  const &  aXI,
            Matrix< Matrix_Type >        &  aYI )
    ->decltype( interp1(aX,aY,aXI.matrix_data(),aYI.matrix_data()) )
    {
        return arma::interp1( aX, aY, aXI.matrix_data(), aYI.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INTERP1_ARMA_HPP_ */

