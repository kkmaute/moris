// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::op_less",
        "[linalgebra],[Mat],[op_less]" )
{

#include "linalg/op_less.inc"

    SECTION(
            "moris::Mat < moris::Mat" )
    {
        REQUIRE( Cm( 0, 0 ) == 1 );
        REQUIRE( Cm( 0, 1 ) == 0 );
        REQUIRE( Cm( 0, 2 ) == 0 );

        REQUIRE( Cm( 1, 0 ) == 1 );
        REQUIRE( Cm( 1, 1 ) == 0 );
        REQUIRE( Cm( 1, 2 ) == 0 );

        REQUIRE( Cm( 2, 0 ) == 1 );
        REQUIRE( Cm( 2, 1 ) == 1 );
        REQUIRE( Cm( 2, 2 ) == 0 );
    }

    SECTION(
            "moris::Col < moris::Col" )
    {
        REQUIRE( Cc( 0, 0 ) == 0 );
        REQUIRE( Cc( 1, 0 ) == 1 );
        REQUIRE( Cc( 2, 0 ) == 0 );
    }

    SECTION( "moris::Col < scalar" )
     {
         moris::Mat< moris::real > A( 3, 1 );
         moris::real B;

         A( 0, 0 ) = 1.0;
         A( 1, 0 ) = 2.0;
         A( 2, 0 ) = 3.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( A < B );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 1 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
     }

     SECTION( "moris::Mat < scalar" )
     {
         moris::Mat< moris::real > A( 3, 2 );
         moris::real B;

         A( 0, 0 ) = 1.0;    A( 0, 1 ) = 4.0;
         A( 1, 0 ) = 2.0;    A( 1, 1 ) = 0.0;
         A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( A < B );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 2 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 0, 1 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 1 ), 0 ) );
     }

     SECTION( "scalar < moris::Col" )
     {
         moris::Mat< moris::real > A( 3, 1 );
         moris::real B;

         A( 0, 0 ) = 1.0;
         A( 1, 0 ) = 2.0;
         A( 2, 0 ) = 3.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( B < A );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 1 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
     }

     SECTION( "scalar < moris::Mat " )
     {
         moris::Mat< moris::real > A( 3, 2 );
         moris::real B;

         A( 0, 0 ) = 1.0;    A( 0, 1 ) = 1.0;
         A( 1, 0 ) = 2.0;    A( 1, 1 ) = 3.0;
         A( 2, 0 ) = 3.0;    A( 2, 1 ) = 2.0;

         B = 2.0;

         moris::Mat< moris::uint > C = ( B < A );

         REQUIRE( moris::equal_to( C.size(0), 3 ) );
         REQUIRE( moris::equal_to( C.size(1), 2 ) );

         REQUIRE( moris::equal_to( C( 0, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 0 ), 0 ) );
         REQUIRE( moris::equal_to( C( 2, 0 ), 1 ) );
         REQUIRE( moris::equal_to( C( 0, 1 ), 0 ) );
         REQUIRE( moris::equal_to( C( 1, 1 ), 1 ) );
         REQUIRE( moris::equal_to( C( 2, 1 ), 0 ) );
     }
}
