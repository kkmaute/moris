/*
 * op_move.cpp
 *
 *  Created on: Nov 9, 2018
 *      Author: sonne
 */

#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_move.hpp"

namespace moris
{
TEST_CASE(
         "moris::op_move",
         "[linalgebra],[op_move]" )
		{
	Matrix< DDRMat > Amatrix(3,3);
	Matrix< DDRMat > Bmatrix(3,3);
	Matrix< DDRMat > Cmatrix(3,3);

	Amatrix(0,0)=1.0; Amatrix(0,1)=2.0; Amatrix(0,2)=3.0;
	Amatrix(1,0)=4.0; Amatrix(1,1)=5.0; Amatrix(2,2)=6.0;

	Bmatrix = Amatrix;
	Cmatrix = move(Bmatrix);

	REQUIRE( Cmatrix == Amatrix );
	}
}
