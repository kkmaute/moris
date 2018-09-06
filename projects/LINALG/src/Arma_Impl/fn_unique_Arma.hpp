/*
 * fn_unique_Arma.hpp
 *
 *  Created on: Sep 5, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_UNIQUE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_UNIQUE_ARMA_HPP_


#include "cl_Matrix.hpp"
#include <armadillo>

namespace moris
{

    template<typename ET >
    auto
    norm( const ET &  aMatrix )
    ->decltype( arma::unique( aMatrix ) )
    {
        return arma::unique( aMatrix );
    }

}


#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_UNIQUE_ARMA_HPP_ */
