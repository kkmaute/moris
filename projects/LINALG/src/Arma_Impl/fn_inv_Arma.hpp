/*
 * fn_inv_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INV_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INV_ARMA_HPP_
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

namespace moris
{
    template< typename ET >
    auto
    inv( const ET & aA )
    ->decltype( arma::inv(aA) )
    {
        return arma::inv(aA);
    }


}



#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_INV_ARMA_HPP_ */
