/*
 * cl_Geometry_Engine.cpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

// Unit Test Suite
#include "catch.hpp"
#include <memory>

#include "ios/cl_Logger.hpp"
#include "core/xtk_typedefs.hpp"


//XTKL: Linear Algebra
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp" // For print

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Geometry  Includes
#include "geometry/cl_Geometry.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "geomeng/cl_MGE_Geometry_Object.hpp"
#include "geometry/cl_Sphere.hpp"
#include "geometry/cl_Discrete_Level_Set.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"



#include "xtk/cl_XTK_Phase_Table.hpp"

TEST_CASE("Function", "[FUNCTION]")
{
    SECTION("Analytical Levelset Sphere","[SPHERE]"){
// Analytical levelset sphere function parameters
    real tRadius = 1;
    real tXCenter = 0.5;
    real tYCenter = 0.25;
    real tZCenter = 0.3;

// Create an analytical levelset sphere function directly
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    // Set test coordinates
    Mat<real, Default_Matrix_Real> tTestCoordinate(1, 3);
    tTestCoordinate(0, 0) = 0;
    tTestCoordinate(0, 1) = 2;
    tTestCoordinate(0, 2) = 3;

    // Check the value of the function
    real tFieldValue = tLevelsetSphere.evaluate_field_value_with_coordinate(0,tTestCoordinate);
    REQUIRE(tFieldValue==Approx(9.6025));

    // Check the sensitivity at the point
//    std::shared_ptr<Matrix_Base<real,Default_Matrix_Real>> tSensitivityInformation = tLevelsetSphere.evaluate_sensitivity_dx_dp_with_coordinate(tTestCoordinate);

    // Create a matrix of known expected values
    Mat<real,Default_Matrix_Real> tExpectedSensitivity({{-0.326991257596857, 0.3910309435028875, -0.6575959492214292}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
//    REQUIRE(equal_to(tSensitivityInformation, tExpectedSensitivity));


}
}

TEST_CASE("Geometry Engine is intersected analytic level set sphere","[GEOMENG][ANALYTIC]")
{
    // Problem Setup
    // Analytical levelset sphere function parameters
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 0;

    // Create an analytical levelset sphere function directly
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    // Create Geometry Engine
    // Create an analytical levelset sphere function directly
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);


    //        (2)
    //     3--x----2            Element with intersected sides as marked by x
    //     |       |            Edge numbers denoted by ()
    // (3) |       x (1)
    //     |       |
    //     0-------1
    //        (0)

    // Coordinates
    Mat<real, Default_Matrix_Real> tNodeCoordinate(4, 3);

    // Node 0 Coordinates
    (tNodeCoordinate)(0, 0) = 0.0;
    (tNodeCoordinate)(0, 1) = 0.0;
    (tNodeCoordinate)(0, 2) = 0.0;

    // Node 1 Coordinates
    (tNodeCoordinate)(1, 0) = 1.0;
    (tNodeCoordinate)(1, 1) = 0.0;
    (tNodeCoordinate)(1, 2) = 0.0;

    // Node 2 Coordinates
    (tNodeCoordinate)(2, 0) = 1.0;
    (tNodeCoordinate)(2, 1) = 1.0;
    (tNodeCoordinate)(2, 2) = 0.0;

    // Node 3 Coordinates
    (tNodeCoordinate)(3, 0) = 1.0;
    (tNodeCoordinate)(3, 1) = 0.0;
    (tNodeCoordinate)(3, 2) = 0.0;

    tGeometryEngine.create_geometry_objects_for_background_mesh_nodes(tNodeCoordinate);

    SECTION("Intersection check with no interface information_ and an analytic level set sphere"){
    // Node to element connectivity
    Mat<size_t, Default_Matrix_Integer> tNodetoElementConnectivity(1, 4);
    (tNodetoElementConnectivity)(0, 0) = 0;
    (tNodetoElementConnectivity)(0, 1) = 1;
    (tNodetoElementConnectivity)(0, 2) = 2;
    (tNodetoElementConnectivity)(0, 3) = 3;

    // Perform intersection check without interface information
    Cell<Geometry_Object<real, size_t, Default_Matrix_Real,Default_Matrix_Integer>> tGeometryObjects;
    tGeometryEngine.is_intersected(tNodeCoordinate, tNodetoElementConnectivity, (size_t) 0,tGeometryObjects);

    REQUIRE(tGeometryObjects.size() == 1);
    REQUIRE(tGeometryObjects(0).get_parent_entity_index() == 0);
}

    SECTION("Intersection check with interface information on an analytic level set sphere"){
    // For interface information the connectivity required by the geometry engine is node to edge rather than node to element
    // Set node to edge connectivity
    Mat<size_t, Default_Matrix_Integer> tNodetoEdgeConnectivity(4, 2);

    // Edge 0
    (tNodetoEdgeConnectivity)(0, 0) = 0;
    (tNodetoEdgeConnectivity)(0, 1) = 1;

    // Edge 1
    (tNodetoEdgeConnectivity)(1, 0) = 1;
    (tNodetoEdgeConnectivity)(1, 1) = 2;

    // Edge 2
    (tNodetoEdgeConnectivity)(2, 0) = 2;
    (tNodetoEdgeConnectivity)(2, 1) = 3;

    // Edge 3
    (tNodetoEdgeConnectivity)(3, 0) = 3;
    (tNodetoEdgeConnectivity)(3, 1) = 0;

    // Perform intersection check with information about the boundary
    Cell<Geometry_Object<real, size_t, Default_Matrix_Real,Default_Matrix_Integer>> tGeometryObjects;
    tGeometryEngine.is_intersected(tNodeCoordinate, tNodetoEdgeConnectivity, (size_t) 1,tGeometryObjects);

    // Get information from geometry object
    size_t tNumGeometryObjects = tGeometryObjects.size();
    size_t tParentEdgeIndex1 = tGeometryObjects(0).get_parent_entity_index();
    size_t tParentEdgeIndex2 = tGeometryObjects(1).get_parent_entity_index();

    real tLocalCoordinate1 = tGeometryObjects(0).get_interface_lcl_coord();
    real tLocalCoordinate2 = tGeometryObjects(1).get_interface_lcl_coord();

    REQUIRE(tNumGeometryObjects == 2);
    REQUIRE(tParentEdgeIndex1 == 1);
    REQUIRE(tParentEdgeIndex2 == 2);
    REQUIRE(tLocalCoordinate1 == Approx(0.875));
    REQUIRE(tLocalCoordinate2 == Approx(-0.875));
}
}

TEST_CASE("Geometry Engine is intersected discrete level set sphere","[GEOMENG][DISCRETE_2D]")
{
    // Problem Setup
    // Analytical levelset sphere function parameters
    real tRadius = 0.25;
    real tXCenter = 0.5;
    real tYCenter = 0.5;
    real tZCenter = 0;

    // Create an analytical levelset sphere function directly
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    // Discretize the Levelset Sphere as a nodal value field on a mesh
    std::string tPrefix = std::getenv("XTKROOT");
    std::string tMeshFileName = tPrefix + "/TestExoFiles/test_mesh_2d_single_element.e";
    Cell<std::string> tScalarFieldNames = {"LEVEL_SET_SPHERE"};
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;

    Cell<Geometry<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>*> tLevelSetFunctions = {&tLevelSetSphere};
    Discrete_Level_Set<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetMesh(tLevelSetFunctions,tMeshFileName,tScalarFieldNames,tMeshBuilder);

    // Create Geometry Engine
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelSetMesh,tPhaseTable);

    //        (2)
    //     3--x----2            Element with intersected sides as marked by x
    //     |       |            Edge numbers denoted by ()
    // (3) |       x (1)
    //     |       |
    //     0-------1
    //        (0)

    // Coordinates
    Mat<real, Default_Matrix_Real> tNodeCoordinate = tLevelSetMesh.get_level_set_mesh()->get_all_node_coordinates_loc_inds();
    tGeometryEngine.create_geometry_objects_for_background_mesh_nodes(tNodeCoordinate);


    SECTION("Intersection check with no interface information_ and an analytic level set sphere"){
    // Node to element connectivity
    Mat<size_t, Default_Matrix_Integer> tNodetoElementConnectivity(1, 4);
    (tNodetoElementConnectivity)(0, 0) = 0;
    (tNodetoElementConnectivity)(0, 1) = 1;
    (tNodetoElementConnectivity)(0, 2) = 2;
    (tNodetoElementConnectivity)(0, 3) = 3;

    // Perform intersection check without interface information
    Cell<Geometry_Object<real, size_t, Default_Matrix_Real,Default_Matrix_Integer>> tGeometryObjects;
    tGeometryEngine.is_intersected(tNodeCoordinate, tNodetoElementConnectivity, (size_t) 0,tGeometryObjects);

    REQUIRE(tGeometryObjects.size() == 1);
    REQUIRE(tGeometryObjects(0).get_parent_entity_index() == 0);
}

    SECTION("Intersection check with interface information on a discrete level set sphere"){
    // For interface information the connectivity required by the geometry engine is node to edge rather than node to element
    // Set node to edge connectivity
    Mat<size_t, Default_Matrix_Integer> tNodetoEdgeConnectivity(4, 2);

    // Get level set mesh
    std::shared_ptr<mesh::Mesh_Data<real,size_t, Default_Matrix_Real, Default_Matrix_Integer>> tLevelSetMeshData = tLevelSetMesh.get_level_set_mesh();

    // Edge 0
    (tNodetoEdgeConnectivity)(0, 0) = 0;
    (tNodetoEdgeConnectivity)(0, 1) = 1;

    // Edge 1
    (tNodetoEdgeConnectivity)(1, 0) = 1;
    (tNodetoEdgeConnectivity)(1, 1) = 2;

    // Edge 2
    (tNodetoEdgeConnectivity)(2, 0) = 2;
    (tNodetoEdgeConnectivity)(2, 1) = 3;

    // Edge 3
    (tNodetoEdgeConnectivity)(3, 0) = 3;
    (tNodetoEdgeConnectivity)(3, 1) = 0;

    // Perform intersection check with information about the boundary
    Cell<Geometry_Object<real, size_t, Default_Matrix_Real,Default_Matrix_Integer>> tGeometryObjects;
    tGeometryEngine.is_intersected(tNodeCoordinate, tNodetoEdgeConnectivity, (size_t) 1,tGeometryObjects);

    // Get information from geometry object
    size_t tNumGeometryObjects = tGeometryObjects.size();
    size_t tParentEdgeIndex1 = tGeometryObjects(0).get_parent_entity_index();
    size_t tParentEdgeIndex2 = tGeometryObjects(1).get_parent_entity_index();

    real tLocalCoordinate1 = tGeometryObjects(0).get_interface_lcl_coord();
    real tLocalCoordinate2 = tGeometryObjects(1).get_interface_lcl_coord();

    REQUIRE(tNumGeometryObjects == 2);
    // Note these are slightly different from above because of my hard-coded numbering scheme
    // Node 0 in this case is = to node 3 in the analytic test case
    REQUIRE(tParentEdgeIndex1 == 0);
    REQUIRE(tParentEdgeIndex2 == 3);
    REQUIRE(tLocalCoordinate1 == Approx(-0.875));
    REQUIRE(tLocalCoordinate2 == Approx(0.875));
}
}

TEST_CASE("Geometry Engine is intersected discrete level set sphere on a 3d mesh","[GEOMENG][DISCRETE_3D]")
{
    // Problem Setup
    // Analytical levelset sphere function parameters
    real tRadius  = 1;
    real tXCenter = 2.5;
    real tYCenter = 2.5;
    real tZCenter = 2.5;

    // Create an analytical levelset sphere function directly
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    // Discretize the Levelset Sphere as a nodal value field on a mesh
    std::string tPrefix = std::getenv("XTKROOT");
    std::string tMeshFileName = tPrefix + "/TestExoFiles/cube_5x5x5_not_unit.e";
    Cell<std::string> tScalarFieldNames = {"LEVEL_SET_SPHERE"};
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    Cell<Geometry<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>*> tLevelSetFunctions = {&tLevelSetSphere};
    Discrete_Level_Set<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelSetMesh(tLevelSetFunctions,tMeshFileName,tScalarFieldNames,tMeshBuilder);

    // Create Geometry Engine
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelSetMesh,tPhaseTable);
    Mat<real, Default_Matrix_Real> tNodeCoordinate = tLevelSetMesh.get_level_set_mesh()->get_all_node_coordinates_loc_inds();

    tGeometryEngine.create_geometry_objects_for_background_mesh_nodes(tNodeCoordinate);


    SECTION("Intersection check with no interface information_ and an analytic level set sphere on a 3d mesh"){

        std::shared_ptr<mesh::Mesh_Data<real,size_t, Default_Matrix_Real, Default_Matrix_Integer>> tLevelSetMeshData = tLevelSetMesh.get_level_set_mesh();
        size_t tNumEdge = tLevelSetMeshData->get_num_entities(EntityRank::ELEMENT);
        Mat<real, Default_Matrix_Real> tNodeCoordinate = tLevelSetMeshData->get_all_node_coordinates_loc_inds();
        Mat<size_t, Default_Matrix_Integer> tNodetoElemConnectivity(tNumEdge, 8);
        Mat<size_t, Default_Matrix_Integer> tNodetoElemConnRow (1, 8);
        // Get level set mesh


        // Package up node to element connectivity
        // TODO: Nest this in a mesh function
        for (size_t i = 0; i < tNumEdge; i++)
        {
            tNodetoElemConnRow = tLevelSetMeshData->get_entity_connected_to_entity_loc_inds(i, EntityRank::ELEMENT, EntityRank::NODE);
            tNodetoElemConnectivity.set_row(i, tNodetoElemConnRow);
        }

        // Perform intersection check without interface information
        Cell<Geometry_Object<real, size_t, Default_Matrix_Real,Default_Matrix_Integer>> tGeometryObjects;
        tGeometryEngine.is_intersected(tNodeCoordinate, tNodetoElemConnectivity, (size_t) 0,tGeometryObjects);


        REQUIRE(tGeometryObjects.size() == 4);
    }
}

TEST_CASE("Compute dx/dp","[DxDp]")
{
    // Set derivative of phi with respect to an arbitrary design variable
    Mat<real, Default_Matrix_Real> tDPhiADp(2,3);
    tDPhiADp.fill(10.0);
    tDPhiADp(0,0) = 0.0;
    tDPhiADp(1,0) = 3.0;

    Mat<real, Default_Matrix_Real> tDPhiBDp(2,3);
    tDPhiBDp.fill(10.0);
    tDPhiBDp(0,0) = 3.0;
    tDPhiBDp(1,0) = 0.0;

   // Set edge node coordinates
    Mat<real, Default_Matrix_Real> tEdgeNodeCoordinates(2,3);

    (tEdgeNodeCoordinates)(0,0) = 0;  (tEdgeNodeCoordinates)(0,1) =   0;  (tEdgeNodeCoordinates)(0,2) = 0;
    (tEdgeNodeCoordinates)(1,0) = 1;  (tEdgeNodeCoordinates)(1,1) = 1.5;  (tEdgeNodeCoordinates)(1,2) = 2;

    // Set edge node coordinates
     Mat<real, Default_Matrix_Real> tEdgeNodePhi(2,1);
     (tEdgeNodePhi)(0,0) = -2;  (tEdgeNodePhi)(1,0) = 1;
     Mat<real, Default_Matrix_Real> tDxDp(2,3);
     //compute_dx_dp_with_linear_basis(tDPhiADp, tDPhiBDp, tEdgeNodeCoordinates, tEdgeNodePhi,tDxDp);

}

