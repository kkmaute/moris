// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::inv",
		"[linalgebra],[inv],[Mat]" )
{
	SECTION( "inv( moris::Mat )" )
	{
		moris::Mat< moris::real > A( 3, 3 );

		A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
		A( 1, 0 ) = 0.0; A( 1, 1 ) = 1.0; A( 1, 2 ) = 4.0;
		A( 2, 0 ) = 5.0; A( 2, 1 ) = 6.0; A( 2, 2 ) = 0.0;

		moris::Mat< moris::real > B = moris::inv( A );

		REQUIRE( moris::equal_to( B( 0, 0 ), -24.0 ) );
		REQUIRE( moris::equal_to( B( 0, 1 ),  18.0 ) );
		REQUIRE( moris::equal_to( B( 0, 2 ),  5.0 ) );

		REQUIRE( moris::equal_to( B( 1, 0 ),  20.0 ) );
		REQUIRE( moris::equal_to( B( 1, 1 ), -15.0 ) );
		REQUIRE( moris::equal_to( B( 1, 2 ), -4.0 ) );

		REQUIRE( moris::equal_to( B( 2, 0 ),  -5.0 ) );
		REQUIRE( moris::equal_to( B( 2, 1 ),   4.0 ) );
		REQUIRE( moris::equal_to( B( 2, 2 ),  1.0 ) );
	}
}
