///*
// * cl_GeometryObject.cpp
// *
// *  Created on: Jan 4, 2017
// *      Author: doble
// */
//#include <catch.hpp>
//#include <iostream>
//
////MORIS Specific Headers
//#include "algorithms.hpp"
//#include "cl_GeometryObject.hpp" // GEN/src
//
//TEST_CASE("moris::geomeng::GeometryObject",
//          "[moris], [geomeng], [GeometryObject]")
//          {
//
//    SECTION(" Set and Get Functions")
//        {
//        // Setting/getting parent element and mesh tests
//        moris::geomeng::GeometryObject tGeoObj1; // initialize
//        tGeoObj1.set_global_parent_element_id(15); // Set parent element id to 2
//        tGeoObj1.set_local_parent_element_index(2);
//        tGeoObj1.set_parent_mesh_id(3);    // Set parent mesh id to 3
//        REQUIRE(moris::equal_to(tGeoObj1.get_global_parent_element_id(),15));
//        REQUIRE(moris::equal_to(tGeoObj1.get_local_parent_element_index(),2));
//        REQUIRE(moris::equal_to(tGeoObj1.get_parent_mesh_id(),3));
//
//        //Setting/getting node request
//        // Set node request 1
//        moris::Mat<moris::real> tReq1(1,5,0);
//        tReq1(0,0) = 1;   tReq1(0,1) = 2; tReq1(0,2) = 0.0; tReq1(0,3) = 1.0; tReq1(0,4) = -1.0;
//        tGeoObj1.set_node_request(tReq1); // set one request
//        REQUIRE(moris::equal_to(tGeoObj1.get_num_requests(),1)); // Make sure there is only one request
//
//        //Access node request 1
//        moris::Mat<moris::real> tReqCheck(1,5,0);
//        tGeoObj1.get_node_request(0,tReqCheck); // Retrieve request 1 from geometry object
//
//        // Check the node request 1
//        REQUIRE(moris::equal_to(tReqCheck(0,0),1));
//        REQUIRE(moris::equal_to(tReqCheck(0,4),-1.0));
//
//        // Set node request 2
//        moris::Mat<moris::real> tReq2(1,5,0);
//        tReq2(0,0) = 2;   tReq2(0,1) = 7; tReq2(0,2) = 10.0; tReq2(0,3) = 11.0; tReq2(0,4) = -14.0;  // Request parameters
//
//        tGeoObj1.set_node_request(tReq2);                        // Set request 2 in geometry object
//        tGeoObj1.get_node_request(1,tReqCheck);                  // Retrieve request 2 from geometry object
//
//        // Test second request
//        REQUIRE(moris::equal_to(tReqCheck(0,1),7));              // Check parent entity id of request
//        REQUIRE(moris::equal_to(tReqCheck(0,3),11));             // Check local zeta_2 coordinate
//        REQUIRE(moris::equal_to(tGeoObj1.get_num_requests(),2)); // Make sure there are two requests
//
//        //Test fulfill_node_request
//        tGeoObj1.fulfill_node_request(0,48);   // give request 1 node ID 48
//        tGeoObj1.fulfill_node_request(1,5543); // give request 2 node ID 5543
//
//        // Use get_node_ids function
//        moris::Mat<moris::uint> tFulReqTest = tGeoObj1.get_node_ids();
//        REQUIRE(moris::equal_to(tFulReqTest(0,0),48));
//        REQUIRE(moris::equal_to(tFulReqTest(1,0),5543));
//        }
//          }


