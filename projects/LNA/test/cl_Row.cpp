// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
		"moris::Row",
		"[linalgebra],[Row]" )
{
	moris::Mat< moris::real > a( 1, 3 );

	a( 0, 0 ) = 1.0;
	a( 0, 1 ) = 4.0;
	a( 0, 2 ) = 9.0;


	// Size.
	REQUIRE( moris::equal_to( a( 0, 0 ), 1.0 ) );
	REQUIRE( moris::equal_to( a( 0, 1 ), 4.0 ) );
	REQUIRE( moris::equal_to( a( 0, 2 ), 9.0 ) );

	REQUIRE( moris::equal_to( a( 0 ), 1.0 ) );
	REQUIRE( moris::equal_to( a( 1 ), 4.0 ) );
	REQUIRE( moris::equal_to( a( 2 ), 9.0 ) );


	// Out of range.
	REQUIRE_THROWS( a( 0, 3 ) );
	REQUIRE_THROWS( a( 1, 0 ) );

	// Out of range single-index
	REQUIRE_THROWS( a( 3 ) );

}
