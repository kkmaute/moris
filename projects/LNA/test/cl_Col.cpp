// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE( "moris::Col", "[linalgebra], [Col]" )
{
	moris::Mat< moris::real > a( 3, 1 );

	a( 0, 0 ) = 1.0;
	a( 1, 0 ) = 4.0;
	a( 2, 0 ) = 9.0;

	// Size.
	REQUIRE( moris::equal_to( a( 0, 0 ), 1.0 ) );
	REQUIRE( moris::equal_to( a( 1, 0 ), 4.0 ) );
	REQUIRE( moris::equal_to( a( 2, 0 ), 9.0 ) );

	REQUIRE( moris::equal_to( a( 0 ), 1.0 ) );
	REQUIRE( moris::equal_to( a( 1 ), 4.0 ) );
	REQUIRE( moris::equal_to( a( 2 ), 9.0 ) );

	// Out of range.
	REQUIRE_THROWS( a( 3, 0 ) );
	REQUIRE_THROWS( a( 0, 1 ) );
	REQUIRE_THROWS( a( 3 ) );
}
