///*
// * cl_Mesh_Input_Buckets.cpp
// *
// *  Created on: Sep 1, 2017
// *      Author: doble
// */
//
//
//#include "catch.hpp"
//#include "cl_Test_Algorithms.hpp"
//#include "cl_Matrix.hpp" // LNA/src
//#include "linalg_typedefs.hpp"
//#include "cl_Cell.hpp" // CON/src
//#include "cl_Mesh_Element_Bucket.hpp" // STK/src/Hierarchical
//
//
//
//namespace moris
//{
//
//TEST_CASE("MESH ELEMENT BUCKET","[MESH][BUCKET]")
//{
//    /*
//     * Initialize
//     * Note: the numbers provided in the constructor must be accurate to prevent dynamic memory allocation
//     */
//    uint tNumElements = 3;
//    uint tNumNodesPerElement = 4;
//    uint tMaxPartNameLength = 6;
//    uint tMaxNumParts = 2;
//    Element_Bucket tElementBucket(tNumElements,tNumNodesPerElement,tMaxNumParts,tMaxPartNameLength);
//
//    // Make sure it answers that there have been no elements placed in bucket
//    CHECK(!tElementBucket.has_entities());
//
//
//    //Parts
//    Cell< std::string >  tPartNames = {"Part_1","Part_2"};
//    tElementBucket.add_part_names(tPartNames);
//    Cell< std::string >  const & tPartNamesOut = tElementBucket.get_part_names();
//    CHECK(tPartNamesOut(0).compare(tPartNames(0)) == 0);
//    CHECK(tPartNamesOut(1).compare(tPartNames(1)) == 0);
//    CHECK(tElementBucket.get_num_parts() == 2);
//
//    // Element 1
//    uint tElementId1 = 3;
//    Matrix< DDUMat >  tElementNodes1 = {{1,2,3,4}};
//
//    tElementBucket.add_entity(tElementNodes1);
//    tElementBucket.add_entity_id(tElementId1);
//
//    CHECK(tElementBucket.has_entities());
//    CHECK(tElementBucket.get_num_entities_in_bucket()==1);
//    CHECK(tElementBucket.get_entity_id(0) == tElementId1);
//    CHECK(equal_to(tElementBucket.get_entity(0),tElementNodes1));
//
//
//    // Elements 2 and 3
//    uint tElementId2 = 4;
//    Matrix< DDUMat >  tElementNodes2 = {{14,56,33,1}};
//    uint tElementId3 = 7;
//    Matrix< DDUMat >  tElementNodes3 = {{1,2,14,23}};
//
//    // Test adding entities more than one at a time
//    Cell<Matrix< DDUMat > > tElements2and3 = {tElementNodes2,tElementNodes3};
//    Cell<uint> tEntity2and3Ids = {tElementId2,tElementId3};
//    tElementBucket.add_entities(tElements2and3);
//    tElementBucket.add_entity_ids(tEntity2and3Ids);
//
//    CHECK(tElementBucket.get_num_entities_in_bucket()==3);
//    CHECK(equal_to(tElementBucket.get_entity(1),tElementNodes2));
//    CHECK(equal_to(tElementBucket.get_entity(2),tElementNodes3));
//
//    CHECK(tElementBucket.get_entity_id(1) == tElementId2);
//    CHECK(tElementBucket.get_entity_id(2) == tElementId3);
//
//
//
//}
//}
//
//
