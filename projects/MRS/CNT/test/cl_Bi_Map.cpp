/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Bi_Map.cpp
 *
 */

//#include <catch.hpp>
//
//// MORIS project header files.
//#include "core.hpp"
//#include "algorithms.hpp"
//#include "cl_Bi_Map.hpp" // CON/src
//
//// ----------------------------------------------------------------------------
//
//TEST_CASE( "moris::Bi_Map" )
//{
//
//    moris::Bi_Map<std::string, moris::real > testMap;
//
//    testMap.insert("one",   1.0);
//    testMap.insert("two",   2.0);
//    testMap.insert("three", 3.0);
//
//    SECTION("moris::Dist_Map", "right access")
//    {
//
//        REQUIRE( testMap.getRight("one")   == 1.0 );
//        REQUIRE( testMap.getRight("two")   == 2.0 );
//        REQUIRE( testMap.getRight("three") == 3.0 );
//    }
//
//    SECTION("moris::Bi_Map", "left access")
//    {
//        moris::Bi_Map<int, std::string > tMap;
//
//        tMap.insert(1, "first");
//        tMap.insert(4, "fourth");
//        tMap.insert(9, "ninth");
//
//        REQUIRE( tMap.getLeft("first")  == 1 );
//        REQUIRE( tMap.getLeft("fourth") == 4 );
//        REQUIRE( tMap.getLeft("ninth")  == 9 );
//    }
//
//    SECTION("moris::Bi_Map", "fail access")
//    {
//        REQUIRE_THROWS( testMap.getRight("four") );
//
//        REQUIRE_THROWS( testMap.getLeft(4.0) );
//    }
//
//}

