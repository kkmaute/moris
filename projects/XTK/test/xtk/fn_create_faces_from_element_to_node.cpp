/*
 * fn_create_faces_from_element_to_node.cpp
 *
 *  Created on: Jun 19, 2018
 *      Author: ktdoble
 */


#include "catch.hpp"

#include <unordered_map>

// Linear Algebra Includes
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg_typedefs.hpp"


#include "mesh/cl_Mesh_Enums.hpp"

#include "xtk/fn_create_faces_from_element_to_node.hpp"





namespace xtk
{

TEST_CASE("fn_create_faces_from_element_to_node","[CREATE_FACES]")
{


    moris::Mat_New<size_t,Default_Matrix_Integer> tElementToNode(4,4);
    tElementToNode(0,0) = 0;  tElementToNode(0,1) = 1;  tElementToNode(0,2) = 2;  tElementToNode(0,3) = 3;
    tElementToNode(1,0) = 1;  tElementToNode(1,1) = 4;  tElementToNode(1,2) = 2;  tElementToNode(1,3) = 3;
    tElementToNode(2,0) = 2;  tElementToNode(2,1) = 4;  tElementToNode(2,2) = 5;  tElementToNode(2,3) = 3;
    tElementToNode(3,0) = 0;  tElementToNode(3,1) = 2;  tElementToNode(3,2) = 6;  tElementToNode(3,3) = 3;

    enum EntityTopology tElementTopo = EntityTopology::TET_4;

    // Number of nodes
    size_t tNumNodes = 7;

    // Element to face output
    moris::Mat_New<size_t,Default_Matrix_Integer> tElementToFace(4,4);

    // Face to Node output
    moris::Mat_New<size_t,Default_Matrix_Integer> tFaceToNode(12,3);

    // Node to face output
    moris::Mat_New<size_t,Default_Matrix_Integer> tNodeToFace(7,10);

    // Face to Element output
    moris::Mat_New<size_t,Default_Matrix_Integer> tFaceToElement(12,2);

    create_faces_from_element_to_node(tElementTopo,
                                      tNumNodes,
                                      tElementToNode,
                                      tElementToFace,
                                      tFaceToNode,
                                      tNodeToFace,
                                      tFaceToElement);


    moris::Mat_New<size_t,Default_Matrix_Integer> tExpElementToFace ({{0, 1, 2, 3},
                                                           {4, 5, 1, 6},
                                                           {5, 7, 8, 9},
                                                           {2, 10, 11, 12}});


    CHECK(equal_to(tElementToFace,tExpElementToFace));


    moris::Mat_New<size_t,Default_Matrix_Integer> tExpFaceToNode({
        {0, 1, 3},
        {2, 1, 3},
        {0, 2, 3},
        {0, 2, 1},
        {1, 4, 3},
        {2, 4, 3},
        {1, 2, 4},
        {5, 4, 3},
        {2, 5, 3},
        {2, 5, 4},
        {6, 2, 3},
        {0, 6, 3},
        {0, 6, 2}});


    CHECK(equal_to(tFaceToNode,tExpFaceToNode));


    size_t tMax = std::numeric_limits<size_t>::max();
    moris::Mat_New<size_t,Default_Matrix_Integer> tExpNodeToFace(
        {{0, 2, 3, 11, 12, tMax, tMax, tMax},
         {0, 1, 3, 4, 6, tMax, tMax, tMax},
         {1, 2, 3, 5, 6, 8, 9, 10},
         {0, 1, 2, 4, 5, 7, 8, 10},
         {4, 5, 6, 7, 9, tMax, tMax, tMax},
         {7, 8, 9, tMax, tMax, tMax, tMax, tMax},
         {10, 11, 12, tMax, tMax, tMax, tMax, tMax}});



    CHECK(equal_to(tNodeToFace,tExpNodeToFace));


    moris::Mat_New<size_t,Default_Matrix_Integer> tExpFaceToElement(
    {{0, tMax},
        {0, 1},
        {0, 3},
        {0, tMax},
        {1, tMax},
        {1, 2},
        {1, tMax},
        {2, tMax},
        {2, tMax},
        {2, tMax},
        {3, tMax},
        {3, tMax},
        {3, tMax}});


    CHECK(equal_to(tFaceToElement,tExpFaceToElement));
}



}
