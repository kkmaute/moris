// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::op_less_equal",
        "[linalgebra],[Mat],[Col],[op_less_equal]" )
{
    SECTION( "moris::Col <= moris::Col" )
    {
        moris::Mat< moris::real > A( 3, 1 );
        moris::Mat< moris::real > B( 3, 1 );

        A( 0, 0 ) = 1.0;
        A( 1, 0 ) = 0.0;
        A( 2, 0 ) = 3.0;

        B( 0, 0 ) = 1.0;
        B( 1, 0 ) = 2.0;
        B( 2, 0 ) = 1.0;

        moris::Mat< moris::uint > C = ( A <= B );

        REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
        REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
    }

    SECTION( "moris::Mat <= moris::Mat" )
    {
        #include "linalg/op_less_equal.inc"
        REQUIRE( C( 0, 0 ) == 1 );
        REQUIRE( C( 0, 1 ) == 1 );
        REQUIRE( C( 0, 2 ) == 0 );

        REQUIRE( C( 1, 0 ) == 1 );
        REQUIRE( C( 1, 1 ) == 1 );
        REQUIRE( C( 1, 2 ) == 0 );

        REQUIRE( C( 2, 0 ) == 1 );
        REQUIRE( C( 2, 1 ) == 1 );
        REQUIRE( C( 2, 2 ) == 0 );
    }

    SECTION( "moris::Col <= scalar" )
     {
         moris::Mat< moris::real > A( 3, 1 );
         moris::real B;

         A( 0, 0 ) = 1.0;
         A( 1, 0 ) = 2.0;
         A( 2, 0 ) = 3.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( A <= B );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 1 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
     }

     SECTION( "moris::Mat <= scalar" )
     {
         moris::Mat< moris::real > A( 3, 2 );
         moris::real B;

         A( 0, 0 ) = 1.0;    A( 0, 1 ) = 4.0;
         A( 1, 0 ) = 2.0;    A( 1, 1 ) = 0.0;
         A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( A <= B );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 2 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 0, 1 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 1 ), 1 ) );
     }

     SECTION( "scalar <= moris::Col" )
     {
         moris::Mat< moris::real > A( 3, 1 );
         moris::real B;

         A( 0, 0 ) = 1.0;
         A( 1, 0 ) = 2.0;
         A( 2, 0 ) = 3.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( B <= A );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 1 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
     }

     SECTION( "scalar <= moris::Mat " )
     {
         moris::Mat< moris::real > A( 3, 2 );
         moris::real B;

         A( 0, 0 ) = 1.0;    A( 0, 1 ) = 1.0;
         A( 1, 0 ) = 2.0;    A( 1, 1 ) = 3.0;
         A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( B <= A );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 2 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 0, 1 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 1 ), 1 ) );
     }
}
