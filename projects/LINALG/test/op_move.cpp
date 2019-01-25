/*
 * op_move.cpp
 *
 *  Created on: Nov 9, 2018
 *      Author: sonne
 */

#include <catch.hpp>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_move.hpp"
#include "op_equal_equal.hpp" //ALG
#include "fn_all_true.hpp"
#include "fn_isempty.hpp"
#include "fn_print.hpp"

namespace moris
{
TEST_CASE(
         "moris::op_move",
         "[linalgebra],[op_move]" )
		{
	Matrix< DDRMat > Amatrix(300,300,0.0);
	Matrix< DDRMat > Bmatrix(300,300,0.0);
	Matrix< DDRMat > Cmatrix(300,300);

	Amatrix(0,0)=1.0; Amatrix(0,1)=2.0; Amatrix(0,2)=3.0;
	Amatrix(1,0)=4.0; Amatrix(1,1)=5.0; Amatrix(1,2)=6.0;
	Amatrix(2,0)=7.0; Amatrix(2,1)=8.0; Amatrix(2,2)=9.0;

	Bmatrix = Amatrix.copy();

	Cmatrix = move(Bmatrix);

	REQUIRE( all_true(Cmatrix == Amatrix) );

	std::cout<<"BMatrix cols = "<<Bmatrix.n_rows();
	std::cout<<"BMatrix rows = "<<Bmatrix.n_cols();

	REQUIRE(isempty(Bmatrix));

	}
}
