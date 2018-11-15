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


//	std::clock_t    start;

//	start = std::clock();
	Bmatrix = Amatrix.copy();
//	std::cout << "Time Copy: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000000) << " ms" << std::endl;

//	start = std::clock();
	Cmatrix = move(Bmatrix);
//	std::cout << "Time Move: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000000) << " ms" << std::endl;


	arma::Mat<real> tBArma(300,300);
	tBArma.fill(0.0);

//	start = std::clock();
	arma::Mat<real> tCArma = std::move(tBArma);
//	std::cout << "Time Move: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000000) << " ms" << std::endl;





	REQUIRE( all_true(Cmatrix == Amatrix) );

	REQUIRE(isempty(Bmatrix));

	}
}
