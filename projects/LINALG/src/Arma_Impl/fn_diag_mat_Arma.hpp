/*
 * fn_diag_Arma.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_MAT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_MAT_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET>
    auto
    diag_mat(
            ET & aA,
            size_t     const & ak = 0 )
    ->decltype(arma::diagmat( aA ))
    {
        return arma::diagmat( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DIAG_MAT_ARMA_HPP_ */
