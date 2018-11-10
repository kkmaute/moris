/*
 * op_move.hpp
 *
 *  Created on: Nov 9, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_OP_MOVE_HPP_
#define PROJECTS_LINALG_SRC_OP_MOVE_HPP_

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_minus_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_minus_Arma.hpp"
#endif

namespace moris
{

	template< typename Matrix_Type >
	auto
	operatorMove( Matrix< Matrix_Type > & aMatrix )
	{
		Matrix< Matrix_Type > & bMatrix;
		bMatrix = std::move(aMatrix);
		return operatorMove(bMatrix);
	}

}

#endif /* PROJECTS_LINALG_SRC_OP_MOVE_HPP_ */
