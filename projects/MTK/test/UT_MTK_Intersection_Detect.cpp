/*
 * UT_MTK_Intersection_Detect.cpp
 *
 *  Created on: Jul 20, 2021
 *      Author: momo
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
                print( tIntersectedPoints, "tIntersectedPointsSame");

                moris::Matrix < DDRMat > tExpectedTriCoords = { { 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 }, { 0.0, 0.0,0.0, 1.0,0.0, 1.0 } };

                REQUIRE( norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001 );


                // Test different coordinates

                // Assign coords
                moris::Matrix < DDRMat >  tFirstTriCoords = { { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
                moris::Matrix < DDRMat > tSecondTriCoords = { { 0.5, 1.0, 0.0 }, { 0.0, 0.5, 0.5 } };

                //use "EdgeIntersect" method
                tIsDetetc.EdgeIntersect( tFirstTriCoords, tSecondTriCoords, tIntersectedPoints, tIntersectVec);

                print( tIntersectedPoints, "tIntersectedPoints");

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

                print( tPointsXInY, "tPointsXInY");

                // Check if the expected result is the same as the output
                REQUIRE(  tPointsXInY.numel() == 0 );

                // Test points of Tri 2 in Tri 1
                tIsDetetc.PointsXInY(  tSecondTriCoords, tFirstTriCoords, tPointsYInX);

                print( tPointsYInX, "tPointsYInX");

                tExpectedTriCoords = { { 0.5, 1.0 }, { 0.0, 0.5 } };

                // Check if the expected result is the same as the output
                REQUIRE( norm( tPointsYInX - tExpectedTriCoords ) < 0.00000001 );


                // Test "SortandRemove"

                // Apply the method
                tIsDetetc.SortandRemove(tIntersectedPoints);

                print( tIntersectedPoints, "tIntersectedPoints");

                // Expected results
                tExpectedTriCoords = { { 0.25, 0.5, 1.0, 0.5 }, { 0.25, 0.0, 0.5, 0.5 } };

                // Check if the expected result is the same as the output
                REQUIRE(norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001  );

                //Test "Intersect"
                tIsDetetc.Intersect(tFirstTriCoords, tSecondTriCoords, tIntersectedPoints, tIntersectVec);

                tExpectedTriCoords = { { 0.25, 0.5, 1.0, 0.5 }, { 0.25, 0.0, 0.5, 0.5 } };

                // Check if the expected result is the same as the output
                REQUIRE(norm( tIntersectedPoints - tExpectedTriCoords ) < 0.00000001  );


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

                // Coordinates of the two meshes
                moris::Matrix < DDRMat >  tFirstTriCoords = { {-1.0, 1.0, -1.0, 0.0, 1.0, 0.6115, 0.1089 } , { -1.0, -1.0 , 1.0, 0.0, 0.4301, 0.6115, 1.0}};

                moris::Matrix < DDRMat >  tSecondTriCoords = { {-1.0, 1.0, -1.0, 0.0, 0.6115, 1.0 , 0.1089 } , {-1.0, -1.0 , 1.0, 0.0, 0.6115, 0.4301, 1.0}};

                // Local indexing of the nodes Tri 1
                moris::Matrix<moris::IdMat> tFisrtTriNodeIndex = { {1, 4, 2},
                        {1, 3, 4},
                        {2, 6, 5},
                        {4, 6, 2},
                        {4, 7, 6},
                        {3, 7, 4} };

                // Local indexing of the nodes Tri 2
                moris::Matrix<moris::IdMat> tSecondNodeIndex = { {1, 2, 4},
                        {1, 4, 3},
                        {4, 6, 5},
                        {2, 6, 4},
                        {4, 5, 7},
                        {3, 4, 7} };

                moris::Cell< moris::Matrix < DDRMat > > tIntersectionPoints = tIsDetetc.elementwise_bruteforce_search(tFirstTriCoords,
                        tSecondTriCoords,
                        tFisrtTriNodeIndex,
                        tSecondNodeIndex);

                for( uint  i = 0 ; i < tIntersectionPoints.size(); i++ )
                {
                    std::cout<<"i is: "<<i<<std::endl;
                    print(tIntersectionPoints(i), "tIntersectionPoints");
                }

                // Check if the expected result is the same as the output
                moris::Matrix < DDRMat > tIntsersectionArea0 = { {-1.000000000000000e+00,  +1.000000000000000e+00,  +0.000000000000000e+00},
                        {-1.000000000000000e+00,  -1.000000000000000e+00,  +0.000000000000000e+00 } };

                moris::Matrix < DDRMat > tIntsersectionArea1 = { {-1.000000000000000e+00,  +0.000000000000000e+00,  -1.000000000000000e+00},
                        {-1.000000000000000e+00,  +0.000000000000000e+00,  +1.000000000000000e+00} };

                moris::Matrix < DDRMat > tIntsersectionArea2 = { {+6.876218536345440e-01,  +1.000000000000000e+00,  +6.115000000000002e-01},
                        {+2.957461592482173e-01,  +4.301000000000000e-01,  +6.115000000000000e-01} };

                moris::Matrix < DDRMat > tIntsersectionArea3 = { {+0.000000000000000e+00,  +6.876218536345440e-01,  +6.115000000000002e-01},
                        {+0.000000000000000e+00,  +2.957461592482173e-01,  +6.115000000000000e-01} };

                REQUIRE( norm( tIntersectionPoints(0) - tIntsersectionArea0 ) < 0.00000001);
                REQUIRE( norm( tIntersectionPoints(1) - tIntsersectionArea1 ) < 0.00000001);
                REQUIRE( norm( tIntersectionPoints(2) - tIntsersectionArea2 ) < 0.00000001);
                REQUIRE( norm( tIntersectionPoints(3) - tIntsersectionArea3 ) < 0.00000001);


                //

                tFirstTriCoords = { {1.0, 0.4300, 0.6115, 1.0 } , { 1.0, 1.0 , 0.6115, 0.108}};

                tSecondTriCoords = { {-1.0, - 1.0, 1.0, 0.0, 0.4300, 1.0 , 0.6115 } , {-1.0, 1.0 , 1.0, 0.0, -1.0, -0.108, -0.6115}};

                tFisrtTriNodeIndex = { {1, 2, 3},
                        {1, 3, 4}};

                tSecondNodeIndex = { {1, 2, 4},
                        {2, 3, 4},
                        {3, 6, 7},
                        {4, 3, 7},
                        {5, 1, 7},
                        {1, 4, 7} };

                tIntersectionPoints = tIsDetetc.elementwise_bruteforce_search(tFirstTriCoords,
                        tSecondTriCoords,
                        tFisrtTriNodeIndex,
                        tSecondNodeIndex);

                for( uint  i = 0 ; i < tIntersectionPoints.size(); i++ )
                {
                    std::cout<<"i is: "<<i<<std::endl;
                    print(tIntersectionPoints(i), "tIntersectionPoints");
                }


                // Check if the expected result is the same as the output
                tIntsersectionArea0 = { {+6.114999999999999e-01,  +1.000000000000000e+00,  +4.299999999999999e-01},
                        {+6.114999999999999e-01,  +1.000000000000000e+00,  +1.000000000000000e+00 } };

                tIntsersectionArea1 = { {+8.361503546099290e-01,  +1.000000000000000e+00,  +1.000000000000000e+00},
                        {+3.203508274231677e-01,  +1.080000000000000e-01,  +1.000000000000000e+00} };

                tIntsersectionArea2 = { {+6.115000000000000e-01,  +8.361503546099291e-01,  +1.000000000000000e+00},
                        {+6.115000000000000e-01,  +3.203508274231680e-01,  +1.000000000000000e+00} };

                REQUIRE( norm( tIntersectionPoints(0) - tIntsersectionArea0 ) < 0.00000001);
                REQUIRE( norm( tIntersectionPoints(1) - tIntsersectionArea1 ) < 0.00000001);
                REQUIRE( norm( tIntersectionPoints(2) - tIntsersectionArea2 ) < 0.00000001);

                delete tInterpMesh;
                delete tIntegMesh;
            }
        }


    }
}



