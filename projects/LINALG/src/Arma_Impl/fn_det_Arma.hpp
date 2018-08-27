/*
 * fn_det_Arma.hpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DET_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DET_ARMA_HPP_

#include "cl_Matrix.hpp"
#include <armadillo>

template<typename ET >
auto
det( const ET &  A)
->decltype( A.determinant() )
{
    return arma::det( A );
}

template<typename ET >
auto
det( ET &  A)
->decltype( A.determinant() )
{
    return arma::det( A );
}



#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DET_ARMA_HPP_ */
