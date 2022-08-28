/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Dist_Map.cpp
 *
 */

//#include <catch.hpp>
//#include<iostream>
//
//// MORIS project header files.
//#include "core.hpp"
//#include "algorithms.hpp"
//#include "cl_Dist_Map.hpp" // CON/src
//
//// ----------------------------------------------------------------------------
//
//TEST_CASE( "moris::Dist_Map" )
//{
//
//    std::vector<moris::lint> data = {1,200,3,400,5};
//
//    moris::Dist_Map testMap(data);
//
//    SECTION("moris::Dist_Map", "access")
//    {
//        REQUIRE( testMap.glbId(0) == 1   );
//        REQUIRE( testMap.glbId(1) == 200 );
//        REQUIRE( testMap.glbId(2) == 3   );
//        REQUIRE( testMap.glbId(3) == 400 );
//        REQUIRE( testMap.glbId(4) == 5   );
//
//        REQUIRE_THROWS( testMap.glbId(5) );
//    }
//
//    SECTION("moris::Dist_Map", "find")
//    {
//        REQUIRE( testMap.locId(1)   == 0 );
//        REQUIRE( testMap.locId(200) == 1 );
//        REQUIRE( testMap.locId(3)   == 2 );
//        REQUIRE( testMap.locId(400) == 3 );
//        REQUIRE( testMap.locId(5)   == 4 );
//
//        REQUIRE_THROWS( testMap.locId( 10 ) );
//    }
//
//    SECTION("moris::Dist_Map", "insert")
//    {
//        testMap.insert(5, 600);
//
//        REQUIRE( testMap.glbId(5) == 600 );
//        REQUIRE( testMap.locId(600) == 5 );
//
//        REQUIRE_THROWS( testMap.glbId(6) );
//    }
//
//    SECTION("moris::DistMap", "copy")
//    {
//        moris::Dist_Map secondMap;
//
//        secondMap = testMap;
//
//        REQUIRE( secondMap.glbId(0) == 1   );
//        REQUIRE( secondMap.glbId(1) == 200 );
//        REQUIRE( secondMap.glbId(2) == 3   );
//        REQUIRE( secondMap.glbId(3) == 400 );
//        REQUIRE( secondMap.glbId(4) == 5   );
//
//        REQUIRE_THROWS( secondMap.glbId(5) );
//
//        REQUIRE( secondMap.locId(1)   == 0 );
//        REQUIRE( secondMap.locId(200) == 1 );
//        REQUIRE( secondMap.locId(3)   == 2 );
//        REQUIRE( secondMap.locId(400) == 3 );
//        REQUIRE( secondMap.locId(5)   == 4 );
//
//        REQUIRE_THROWS( secondMap.locId(100) );
//    }
//
//    SECTION("moris::DistMap", "copy")
//    {
//        moris::Dist_Map thirdMap(testMap);
//
//
//        REQUIRE( thirdMap.glbId(0) == 1   );
//        REQUIRE( thirdMap.glbId(1) == 200 );
//        REQUIRE( thirdMap.glbId(2) == 3   );
//        REQUIRE( thirdMap.glbId(3) == 400 );
//        REQUIRE( thirdMap.glbId(4) == 5   );
//
//        REQUIRE_THROWS( thirdMap.glbId(5) );
//
//        REQUIRE( thirdMap.locId(1)   == 0 );
//        REQUIRE( thirdMap.locId(200) == 1 );
//        REQUIRE( thirdMap.locId(3)   == 2 );
//        REQUIRE( thirdMap.locId(400) == 3 );
//        REQUIRE( thirdMap.locId(5)   == 4 );
//
//        REQUIRE_THROWS( thirdMap.locId(100) );
//    }
//
//}

