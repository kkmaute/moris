// Third-party header files.
#include <catch.hpp>

// C++ header files
#include "iostream"

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::eig_gen",
		"[linalgebra],[eig_gen],[Mat]" )
{
	SECTION( "eig_gen( moris::Mat )" )
	{
		moris::Mat< moris::real > A( 3, 3 );

		A( 0, 0 ) = -1.0; A( 0, 1 ) = 1.0; A( 0, 2 ) = 3.0;
		A( 1, 0 ) =  1.0; A( 1, 1 ) = 2.0; A( 1, 2 ) = 0.0;
		A( 2, 0 ) =  3.0; A( 2, 1 ) = 0.0; A( 2, 2 ) = 2.0;

		moris::Mat< moris::cplx > eig_vec;
		moris::Mat< moris::cplx > eig_val;

		moris::eig_gen( eig_val, eig_vec, A );
		// NOTE: This test evaluates the values, not the signs. To test/return certain signs and get consistent results,
		// one could require that the first element in each eigenvector be of a certain sign. Then, if you get an
		// eigenvector that does not follow this rule, you multiply it by -1.

		// NOTE: Eigen and Armadillo both return the eigenvalues in the same order, but Matlab returns them in a different order. Does this matter?

		REQUIRE( moris::equal_to( std::abs( eig_val( 0, 0 ) ), 3.000000000000000 ) );
		REQUIRE( moris::equal_to( std::abs( eig_val( 1, 0 ) ), 4.000000000000000 ) );
		REQUIRE( moris::equal_to( std::abs( eig_val( 2, 0 ) ), 2.000000000000000 ) );

		REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 0 ) ), 0.8451542547285165 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 2 ) ), 0.0000000000000000 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 1 ) ), 0.5345224838248487 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 0 ) ), 0.1690308509457034 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 2 ) ), 0.9486832980505139 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 1 ) ), 0.2672612419124244 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 0 ) ), 0.5070925528371100 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 2 ) ), 0.3162277660168378 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 1 ) ), 0.8017837257372730 ) );
	}
}
