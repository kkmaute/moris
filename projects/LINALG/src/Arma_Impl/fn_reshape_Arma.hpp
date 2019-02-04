/*
 * fn_reshape_Arma.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_

#include <armadillo>
#include "cl_Matrix.hpp"

namespace moris
{

template<  typename Matrix_Type, typename INT >

arma::Mat< Matrix_Type >
reshape(
		arma::Mat< Matrix_Type >   	& aA,
		INT 			  			  aB,
		INT				  			  aC )
{
	arma::Mat< Matrix_Type > tM2 = aA;
	tM2.reshape(aB, aC);

    return tM2;
}

} /* moris namespace */

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_ */
