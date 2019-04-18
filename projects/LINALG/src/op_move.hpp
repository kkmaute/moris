/*
 * op_move.hpp
 *
 *  Created on: Nov 9, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_OP_MOVE_HPP_
#define PROJECTS_LINALG_SRC_OP_MOVE_HPP_

#include "cl_Matrix.hpp"

//#ifdef MORIS_USE_EIGEN
//#include "Eigen_Impl/op_move_Eigen.hpp"
//#endif
//
//#ifdef MORIS_USE_ARMA
//#include "op_move_Arma.hpp"
//#endif

namespace moris
{
	template< typename Matrix_Type >
	Matrix<Matrix_Type>
	move( Matrix< Matrix_Type > & aMatrix )
	{
		Matrix< Matrix_Type > tB = std::move(aMatrix);
		aMatrix.resize(0,0);
		return tB;
	}
}

#endif /* PROJECTS_LINALG_SRC_OP_MOVE_HPP_ */
