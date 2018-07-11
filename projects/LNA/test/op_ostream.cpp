//// Third-party header files.
//#include <catch.hpp>
//
//// MORIS project header files.
//#include "linalg.hpp"
//
//// ----------------------------------------------------------------------------
//
//TEST_CASE(
//         "moris::op_ostream",
//         "[linalgebra],[op_ostream]" )
//{
//
//#include "linalg/op_ostream.inc"
//
//    SECTION( "moris::Mat" )
//    {
//        REQUIRE( Am( 0, 0 ) == 1.0 );
//        REQUIRE( Am( 0, 1 ) == 1.0 );
//        REQUIRE( Am( 0, 2 ) == 1.0 );
//
//        REQUIRE( Am( 1, 0 ) == 2.0 );
//        REQUIRE( Am( 1, 1 ) == 2.0 );
//        REQUIRE( Am( 1, 2 ) == 2.0 );
//
//        REQUIRE( Am( 2, 0 ) == 3.0 );
//        REQUIRE( Am( 2, 1 ) == 3.0 );
//        REQUIRE( Am( 2, 2 ) == 3.0 );
//    }
//
//    SECTION( "moris::Sp_Mat" )
//    {
//        REQUIRE( As( 0, 0 ) == 0.0 );
//        REQUIRE( As( 0, 1 ) == 1.0 );
//        REQUIRE( As( 0, 2 ) == 1.0 );
//
//        REQUIRE( As( 1, 0 ) == 2.0 );
//        REQUIRE( As( 1, 1 ) == 0.0 );
//        REQUIRE( As( 1, 2 ) == 2.0 );
//
//        REQUIRE( As( 2, 0 ) == 3.0 );
//        REQUIRE( As( 2, 1 ) == 3.0 );
//        REQUIRE( As( 2, 2 ) == 0.0 );
//    }
//}
