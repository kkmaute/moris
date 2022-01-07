/*
 * fn_interp1_Arma.hpp
 *
 *  Created on: Sep. 21, 2021
 *      Author: gates
 *
 * This subroutine performs the armadillo interp1 function as described in:
 * http://arma.sourceforge.net/docs.html#interp1
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
