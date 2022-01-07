/*
 * fn_ctrans_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CTRANS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CTRANS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    ctrans( const ET & aA )
    ->decltype( arma::trans(aA) )
    {
        return arma::trans(aA);
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CTRANS_ARMA_HPP_ */
