// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::eig_sym",
		"[linalgebra],[eig_sym],[Mat]" )
{
	SECTION( "eig_sym( moris::Mat )" )
	{
		moris::Mat< moris::real > A( 3, 3 );

		A( 0, 0 ) = -1.0; A( 0, 1 ) = 3.0; A( 0, 2 ) = 0.0;
		A( 1, 0 ) =  3.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 2.0;
		A( 2, 0 ) =  0.0; A( 2, 1 ) = 2.0; A( 2, 2 ) = 4.0;

		moris::Mat< moris::real > eig_vec;
		moris::Mat< moris::real > eig_val;

		moris::eig_sym( eig_val, eig_vec, A );

		// NOTE: This test evalues the values, not the signs. To test/return certain signs and get consistent results,
		// one could require that the first element in each eigenvector be of a certain sign. Then, if you get an
		// eigenvector that does not follow this rule, you multiply it by -1.

		REQUIRE( moris::equal_to( std::abs( eig_val( 0, 0 ) ), 2.341200659105694 ) );
		REQUIRE( moris::equal_to( std::abs( eig_val( 1, 0 ) ), 3.043564371011120 ) );
		REQUIRE( moris::equal_to( std::abs( eig_val( 2, 0 ) ), 7.297636288094574 ) );

		REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 0 ) ), 0.905449906553528 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 1 ) ), 0.304846459001434 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 0, 2 ) ), 0.295345734955654 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 0 ) ), 0.404796670485593 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 1 ) ), 0.410888760082367 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 1, 2 ) ), 0.816890495967332 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 0 ) ), 0.127671932256022 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 1 ) ), 0.859208393390255 ) );
		REQUIRE( moris::equal_to( std::abs( eig_vec( 2, 2 ) ), 0.495440021033577 ) );
	}
}
