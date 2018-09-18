/*
 * fn_generate_element_to_element.cpp
 *
 *  Created on: Feb 8, 2018
 *      Author: ktdoble
 */


#include "catch.hpp"

// Linear Algebra Includes
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"

#include"xtk/fn_generate_element_to_element.hpp"




namespace xtk
{

/*
 * This Test case checks borth sequential and nonsequential versions of generated element to element connectivity
 */
TEST_CASE("Generate Element to Element Connectivity Sequential","[ELEM_TO_ELEM]")
{
    // This test starts from a face to element connectivity and generates the element to element connectivity
    // using the function generate_element_to_element();

    SECTION("On a Tetrahedral Mesh")
    {
        // Using the tet from cl_XTK_Cut_Mesh_Modification.cpp

        // There are 8 elements in this problem
        size_t tNumElements = 8;

        // Each element has 4 faces
        size_t tNumFacesPerElement = 4;

        // The face to element connectivity is as follow:
        size_t tMax = std::numeric_limits<size_t>::max();
        moris::Matrix< Default_Matrix_Integer > tFaceToElement(23,2);
        tFaceToElement.fill(tMax);
        tFaceToElement(0,0)  = 0;
        tFaceToElement(1,0)  = 3;
        tFaceToElement(2,0)  = 0;
        tFaceToElement(3,0)  = 0;    tFaceToElement(3,1)  = 1;
        tFaceToElement(4,0)  = 1;
        tFaceToElement(5,0)  = 6;
        tFaceToElement(6,0)  = 1;
        tFaceToElement(7,0)  = 0;    tFaceToElement(7,1)  = 2;
        tFaceToElement(8,0)  = 2;    tFaceToElement(8,1)  = 4;
        tFaceToElement(9,0)  = 2;
        tFaceToElement(10,0) = 2;
        tFaceToElement(11,0) = 3;    tFaceToElement(11,1) = 6;
        tFaceToElement(12,0) = 3;    tFaceToElement(12,1) = 4;
        tFaceToElement(13,0) = 3;
        tFaceToElement(14,0) = 4;
        tFaceToElement(15,0) = 4;    tFaceToElement(15,1) = 7;
        tFaceToElement(16,0) = 1;    tFaceToElement(16,1) = 5;
        tFaceToElement(17,0) = 5;
        tFaceToElement(18,0) = 5;    tFaceToElement(18,1) = 7;
        tFaceToElement(19,0) = 5;
        tFaceToElement(20,0) = 6;
        tFaceToElement(21,0) = 6;    tFaceToElement(21,1) = 7;
        tFaceToElement(22,0) = 7;


        // Run the generate element connectivity
        moris::Matrix< Default_Matrix_Integer > tElementToElement = generate_element_to_element(tFaceToElement,tNumElements,tNumFacesPerElement,tMax);

        // Expected result
        moris::Matrix< Default_Matrix_Integer > tExpectedElementToElement(8,4);
        tExpectedElementToElement.fill(tMax);
        // Element 0 neighbors
        tExpectedElementToElement(0,0) = 1; tExpectedElementToElement(0,1) = 2;
        // Element 1 neighbors
        tExpectedElementToElement(1,0) = 0; tExpectedElementToElement(1,1) = 5;
        // Element 2 neighbors
        tExpectedElementToElement(2,0) = 0; tExpectedElementToElement(2,1) = 4;
        // Element 3 neighbors
        tExpectedElementToElement(3,0) = 6; tExpectedElementToElement(3,1) = 4;
        // Element 4 neighbors
        tExpectedElementToElement(4,0) = 2; tExpectedElementToElement(4,1) = 3; tExpectedElementToElement(4,2) = 7;
        // Element 5 neighbors
        tExpectedElementToElement(5,0) = 1; tExpectedElementToElement(5,1) = 7;
        // Element 6 neighbors
        tExpectedElementToElement(6,0) = 3; tExpectedElementToElement(6,1) = 7;
        // Element 7 neighbors
        tExpectedElementToElement(7,0) = 4; tExpectedElementToElement(7,1) = 5; tExpectedElementToElement(7,2) = 6;

        // Make sure the neighbor relations are correct

        CHECK(equal_to(tElementToElement,tExpectedElementToElement));
    }
}

TEST_CASE("Generate Element to Element Connectivity Non Sequential","[ELEM_TO_ELEM]")
{
    // This test starts from a face to element connectivity and generates the element to element connectivity
    // using the function generate_element_to_element_nonsequential();

    // Using the tet from cl_XTK_Cut_Mesh_Modification.cpp

    // There are 8 elements in this problem
    size_t tNumElements = 8;

    // Each element has 4 faces
    size_t tNumFacesPerElement = 4;

    // Element Index vecot
    moris::Matrix< Default_Matrix_Integer > tElementVector(1,8);

   tElementVector(0,0)=  3;
   tElementVector(0,1)=  5;
   tElementVector(0,2)=  6;
   tElementVector(0,3)=  7;
   tElementVector(0,4)= 10;
   tElementVector(0,5)= 16;
   tElementVector(0,6)= 19;
   tElementVector(0,7)=  1;


   // Generate Map
   std::unordered_map<size_t,size_t> tElementMap(tNumElements);
   for(size_t iE = 0; iE<tNumElements; iE++)
   {
       tElementMap[tElementVector(0,iE)] = iE;
   }

   // The face to element connectivity is as follow:
    size_t tMax = std::numeric_limits<size_t>::max();
    moris::Matrix< Default_Matrix_Integer > tFaceToElement(23,2);
    tFaceToElement.fill(tMax);
    tFaceToElement(0,0)  = tElementVector(0,0);
    tFaceToElement(1,0)  = tElementVector(0,3);
    tFaceToElement(2,0)  = tElementVector(0,0);
    tFaceToElement(3,0)  = tElementVector(0,0);    tFaceToElement(3,1)  = tElementVector(0,1);
    tFaceToElement(4,0)  = tElementVector(0,1);
    tFaceToElement(5,0)  = tElementVector(0,6);
    tFaceToElement(6,0)  = tElementVector(0,1);
    tFaceToElement(7,0)  = tElementVector(0,0);    tFaceToElement(7,1)  = tElementVector(0,2);
    tFaceToElement(8,0)  = tElementVector(0,2);    tFaceToElement(8,1)  = tElementVector(0,4);
    tFaceToElement(9,0)  = tElementVector(0,2);
    tFaceToElement(10,0) = tElementVector(0,2);
    tFaceToElement(11,0) = tElementVector(0,3);    tFaceToElement(11,1) = tElementVector(0,6);
    tFaceToElement(12,0) = tElementVector(0,3);    tFaceToElement(12,1) = tElementVector(0,4);
    tFaceToElement(13,0) = tElementVector(0,3);
    tFaceToElement(14,0) = tElementVector(0,4);
    tFaceToElement(15,0) = tElementVector(0,4);    tFaceToElement(15,1) = tElementVector(0,7);
    tFaceToElement(16,0) = tElementVector(0,1);    tFaceToElement(16,1) = tElementVector(0,5);
    tFaceToElement(17,0) = tElementVector(0,5);
    tFaceToElement(18,0) = tElementVector(0,5);    tFaceToElement(18,1) = tElementVector(0,7);
    tFaceToElement(19,0) = tElementVector(0,5);
    tFaceToElement(20,0) = tElementVector(0,6);
    tFaceToElement(21,0) = tElementVector(0,6);    tFaceToElement(21,1) = tElementVector(0,7);
    tFaceToElement(22,0) = tElementVector(0,7);

    // Expected result
    moris::Matrix< Default_Matrix_Integer > tExpectedElementToElement(8,4);
    tExpectedElementToElement.fill(tMax);
    tExpectedElementToElement(0,0) = tElementMap[tElementVector(0,1)]; tExpectedElementToElement(0,1) = tElementMap[tElementVector(0,2)];
    tExpectedElementToElement(1,0) = tElementMap[tElementVector(0,0)]; tExpectedElementToElement(1,1) = tElementMap[tElementVector(0,5)];
    tExpectedElementToElement(2,0) = tElementMap[tElementVector(0,0)]; tExpectedElementToElement(2,1) = tElementMap[tElementVector(0,4)];
    tExpectedElementToElement(3,0) = tElementMap[tElementVector(0,6)]; tExpectedElementToElement(3,1) = tElementMap[tElementVector(0,4)];
    tExpectedElementToElement(4,0) = tElementMap[tElementVector(0,2)]; tExpectedElementToElement(4,1) = tElementMap[tElementVector(0,3)]; tExpectedElementToElement(4,2) = tElementMap[tElementVector(0,7)];
    tExpectedElementToElement(5,0) = tElementMap[tElementVector(0,1)]; tExpectedElementToElement(5,1) = tElementMap[tElementVector(0,7)];
    tExpectedElementToElement(6,0) = tElementMap[tElementVector(0,3)]; tExpectedElementToElement(6,1) = tElementMap[tElementVector(0,7)];
    tExpectedElementToElement(7,0) = tElementMap[tElementVector(0,4)]; tExpectedElementToElement(7,1) = tElementMap[tElementVector(0,5)]; tExpectedElementToElement(7,2) = tElementMap[tElementVector(0,6)];


    // Run the generate element connectivity
    moris::Matrix< Default_Matrix_Integer > tElementToElement =
            generate_element_to_element_nonsequential(tFaceToElement, tElementMap, tNumElements, tNumFacesPerElement, tMax);


    // Make sure the neighbor relations are correct
    CHECK(equal_to(tElementToElement,tExpectedElementToElement));

}
}

