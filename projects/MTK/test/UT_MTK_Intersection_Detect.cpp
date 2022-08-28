/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Intersection_Detect.cpp
 *
 */

#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
//#//include "cl_XTK_Periodic_Boundary_Condition_Helper.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Set.hpp" //MTK/src
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "catch.hpp"
#include "paths.hpp"

// implementations to test
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Intersection_Detect.cpp"

namespace moris
{
    namespace mtk
    {
        TEST_CASE("MTK Intersection Polygon Clipping","[MTK],[MTK_Intersect_Polygon_Clipping]")
                                                                            {
            if(par_size() ==1)
            {
                //generate a cubic mesh
                std::string tInterpString = "generated:2x1x1|sideset:yY";

                //create interpolation integration mesh
                mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tInterpString );
                mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh );

                // create a mesh manager and register the mesh pair
                auto tMeshManager = std::make_shared< mtk::Mesh_Manager >();
                tMeshManager->register_mesh_pair( tInterpMesh, tIntegMesh );

                //parameter list input for the surfaces that are going to be periodic
                moris::ParameterList tParameterList;
                tParameterList.insert( "periodic_side_set_pair", "surface_1,surface_2");

                //construct the intersection detect
                mtk::Intersection_Detect tIsDetetc( tMeshManager, 0, tParameterList, 2) ;

                // Test "EdgeIntersect"
                // Test same coordinates

                // Assign coords
                moris::Matrix < DDRMat >   tFirstTriCoordsSame = { { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
                moris::Matrix < DDRMat >   tSecondTriCoordsSame = { { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };

                // Initialize the outputs
                moris::Matrix < DDRMat >  tIntersectedPoints;
                moris::Matrix < DDUMat >  tIntersectVec;

                //use "EdgeIntersect" method
                tIsDetetc.EdgeIntersect( tFirstTriCoordsSame, tSecondTriCoordsSame, tIntersectedPoints, tIntersectVec);

                moris::Matrix < DDRMat > tExpectedTriCoords = { { 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 }, { 0.0, 0.0,0.0, 1.0,0.0, 1.0 } };

                REQUIRE( norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001 );

                // Test different coordinates

                // Assign coords
                moris::Matrix < DDRMat >  tFirstTriCoords = { { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
                moris::Matrix < DDRMat > tSecondTriCoords = { { 0.5, 1.0, 0.0 }, { 0.0, 0.5, 0.5 } };

                //use "EdgeIntersect" method
                tIsDetetc.EdgeIntersect( tFirstTriCoords, tSecondTriCoords, tIntersectedPoints, tIntersectVec);

                // Expected result
                tExpectedTriCoords = { { 0.5, 0.5, 1.0, 1.0, 0.5, 0.25 }, { 0.0, 0.0, 0.5, 0.5, 0.5, 0.25 } };

                // Check if the expected result is the same as the output
                REQUIRE( norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001 );

                // Test "PointsXInY"

                // Initialize outputs
                moris::Matrix < DDRMat > tPointsXInY;
                moris::Matrix < DDRMat > tPointsYInX;

                // Test points of Tri 1 in Tri 2
                tIsDetetc.PointsXInY( tFirstTriCoords , tSecondTriCoords, tPointsXInY);

                // Check if the expected result is the same as the output
                REQUIRE(  tPointsXInY.numel() == 0 );

                // Test points of Tri 2 in Tri 1
                tIsDetetc.PointsXInY(  tSecondTriCoords, tFirstTriCoords, tPointsYInX);

                tExpectedTriCoords = { { 0.5, 1.0 }, { 0.0, 0.5 } };

                // Check if the expected result is the same as the output
                REQUIRE( norm( tPointsYInX - tExpectedTriCoords ) < 0.00000001 );

                // Test "SortandRemove"

                // Apply the method
                tIsDetetc.SortandRemove(tIntersectedPoints);

                // Expected results
                tExpectedTriCoords = { { 0.25, 0.5, 1.0, 0.5 }, { 0.25, 0.0, 0.5, 0.5 } };

                // Check if the expected result is the same as the output
                REQUIRE(norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001  );

                //Test "Intersect"
                tIsDetetc.Intersect(tFirstTriCoords, tSecondTriCoords, tIntersectedPoints, tIntersectVec);

                tExpectedTriCoords = { { 0.25, 0.5, 1.0, 0.5 }, { 0.25, 0.0, 0.5, 0.5 } };

                // Check if the expected result is the same as the output
                REQUIRE(norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001  );

                std::cout<<"Hello World"<<std::endl;

                delete tInterpMesh;
                delete tIntegMesh;
            }
                                                                            }

        TEST_CASE("MTK Intersection Intersection","[MTK],[MTK_Interface]")
        {

            if(par_size() ==1)
            {

                //generate a cubic mesh
                std::string tInterpString = "generated:2x1x1|sideset:yY";

                //create interpolation integration mesh
                mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tInterpString );
                mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh );

                // create a mesh manager and register the mesh pair
                auto tMeshManager = std::make_shared< mtk::Mesh_Manager >();
                tMeshManager->register_mesh_pair( tInterpMesh, tIntegMesh );

                //parameter list input for the surfaces that are going to be periodic
                moris::ParameterList tParameterList;
                tParameterList.insert( "periodic_side_set_pair", "surface_1,surface_2");

                //construct the intersection detect
                mtk::Intersection_Detect tIsDetetc( tMeshManager, 0, tParameterList, 2) ;

                //Test Interface Matrix

                moris::Matrix < DDRMat > tFirstMesh1 = { { 1.0, -1.0, -1.0},{ -1.0, -1.0, 1.0} };
                moris::Matrix < DDRMat > tFirstMesh2 = { { 1.0, -1.0, +1.0},{ -1.0, 1.0, 1.0} };

                moris::Matrix < DDRMat > tSecondMesh1 = { { -1.0, +1.0, +1.0},{ -1.0, +1.0, -1.0} };
                moris::Matrix < DDRMat > tSecondMesh2 = { { -1.0, -1.0, +1.0},{ -1.0, 1.0, 1.0} };

                moris::Cell< moris::Matrix < DDRMat > >  tFirstMesh = { tFirstMesh1,  tFirstMesh2};
                moris::Cell< moris::Matrix < DDRMat > >  tSecondMesh = { tSecondMesh1,  tSecondMesh2};

                Matrix<IndexMat> tFirstMeshIdentifier = { {0,1} };
                Matrix<IndexMat> tSecondMeshIdentifier = { {0,1} };

                moris::Cell < moris::Matrix <DDRMat> >          tCutTriangles;
                moris::Matrix< moris::IndexMat>                tCutTrianglesIdentifier;

                tIsDetetc.elementwise_bruteforce_search(tFirstMesh,
                        tFirstMeshIdentifier,
                        tSecondMesh,
                        tSecondMeshIdentifier,
                        tCutTriangles,
                        tCutTrianglesIdentifier);

                for(uint i = 0 ; i < tCutTriangles.size(); i++ )
                {
                    print(tCutTriangles(i), "tCutTriangles");
                }
                // Check if the expected result is the same as the output
                moris::Matrix < DDRMat > tIntsersectionArea00 = { {-1.0, 1.0 , 0.0}, {-1.0, -1.0 ,0.0 } };

                moris::Matrix < DDRMat > tIntsersectionArea01 = { {-1.0, -1.0 , 0.0}, {-1.0, +1.0 ,0.0 } };

                moris::Matrix < DDRMat > tIntsersectionArea10 = { {0.0 , 1.0 , 1.0}, {0.0, +1.0 ,-1.0 } };

                moris::Matrix < DDRMat > tIntsersectionArea11 = { {-1.0, +1.0 , 0.0}, {+1.0, +1.0 ,0.0 } };

                REQUIRE( norm( tCutTriangles( 0 ) ) - norm( tIntsersectionArea00 ) < 0.00000001);
                REQUIRE( norm( tCutTriangles( 1 ) ) - norm( tIntsersectionArea01 ) < 0.00000001);
                REQUIRE( norm( tCutTriangles( 2 ) ) - norm( tIntsersectionArea10 ) < 0.00000001);
                REQUIRE( norm( tCutTriangles( 3 ) ) - norm( tIntsersectionArea11 ) < 0.00000001);

                //

                delete tInterpMesh;
                delete tIntegMesh;
            }
        }

    }
}

