/*
 * cl_XTK_Mesh.cpp
 *
 *  Created on: Jun 21, 2017
 *      Author: ktdoble
 */

#include <memory>   // Shared ptrs
#include "catch.hpp"

// XTKL: Logging and Assertion Includes
#include "ios/cl_Logger.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"


#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "geometry/cl_Gyroid.hpp"
#include "geometry/cl_Sphere.hpp"
#include "catch.hpp"
#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "linalg_typedefs.hpp"
#include "tools/fn_tet_volume.hpp"

//XTK  Includes
#include "xtk/cl_XTK_Cut_Mesh.hpp"

#include "xtk/cl_XTK_Child_Mesh.hpp"
#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"

// XTKL: Geometry  Include
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Tools.hpp"

namespace xtk
{
size_t find_edge_in_base_mesh(moris::Matrix<size_t, Default_Matrix_Integer> const & aBaseEdges,
                              moris::Matrix<size_t, Default_Matrix_Integer> const & aEdgeToFind)
{
    size_t tNumEdges = aBaseEdges.n_rows();
    for(size_t i = 0; i<tNumEdges; i++)
    {
        if(aBaseEdges(i,0) == aEdgeToFind(0,0) || aBaseEdges(i,0) == aEdgeToFind(0,1))
        {
            if(aBaseEdges(i,1) == aEdgeToFind(0,0) || aBaseEdges(i,1) == aEdgeToFind(0,1))
            {
                return i;
            }
        }
    }

    return 10000;

}

TEST_CASE("Simple Mesh Testing","[XTK][CUT_MESH]"){
    // Functionality tested in this case
    // - get_num_simple_meshes
    // - set/get_node_inds/ids
    // TODO Add add_index function test

    //       Indices
    //    9----10----11
    //   /|    /|    /|    // Element 1 Indices - 0,1,4,3,6,7,10,9
    //  6-----7-----8 |    // Element 2 Indices - 1,2,5,4,7,8,11,10
    //  | |   | |   | |
    //  | 3---|-4---|-5
    //  |/    |/    |/
    //  0-----1-----2
    //

    //         Ids
    //    2----10----11
    //   /|    /|    /|
    //  7-----36----8 |    // Element 1 Ids - 14,5,18,4,7,36,10,2
    //  | |   | |   | |    // Element 2 Ids - 5,3,8,18,36,8,11,10
    //  | 4---|18---|-8
    //  |/    |/    |/
    //  14----5-----3

    // Initialize an Cut Mesh with 2 simple meshes that are 3d
    size_t tModelDim = 3;
    size_t tNumSimpleMesh = 2;
    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tCutMesh(tNumSimpleMesh, tModelDim);

    // Make sure the correct number of simple meshes have been created
    REQUIRE(tCutMesh.get_num_simple_meshes() == 2);

    // Add node Indices then node ids for each element
    Child_Mesh_Test<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCM1 = tCutMesh.get_child_mesh(0);// Index of element 1
    Child_Mesh_Test<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCM2 = tCutMesh.get_child_mesh(1);// Index of element 2
    moris::Matrix<size_t,Default_Matrix_Integer> tInds1({{0, 1, 4, 3, 6, 7, 10, 9}}); // Indices of element 1
    moris::Matrix<size_t,Default_Matrix_Integer> tIds1({{14, 5, 18, 4, 7, 36, 10, 2}}); // Ids of element 1
    moris::Matrix<size_t,Default_Matrix_Integer> tInds2({{1, 2, 5, 4, 7, 8, 11, 10}}); // Indices of element 2
    moris::Matrix<size_t,Default_Matrix_Integer> tIds2({{   5, 3, 8, 18, 36, 8, 11, 10}}); // Ids of element 2

    // Set node indices
    tCM1.add_node_indices(tInds1);
    tCM2.add_node_indices(tInds2);

    // Set node Ids
    tCM1.add_node_ids(tInds1);
    tCM2.add_node_ids(tInds2);

    //Test node indices for element 1
    moris::Matrix<size_t,Default_Matrix_Integer> const & tNodeInds1 = tCM1.get_node_indices();
    moris::Matrix<size_t,Default_Matrix_Integer> const & tNodeInds2 = tCM2.get_node_indices();
    moris::Matrix<size_t,Default_Matrix_Integer> const & tNodeIds1 = tCM1.get_node_ids();
    moris::Matrix<size_t,Default_Matrix_Integer> const & tNodeIds2 = tCM2.get_node_ids();

    REQUIRE(tNodeInds1(0, 0) == Approx(0));
    REQUIRE(tNodeInds1(0, 1) == Approx(1));
    REQUIRE(tNodeInds1(0, 2) == Approx(4));
    REQUIRE(tNodeInds1(0, 3) == Approx(3));
    REQUIRE(tNodeInds1(0, 4) == Approx(6));
    REQUIRE(tNodeInds1(0, 5) == Approx(7));
    REQUIRE(tNodeInds1(0, 6) == Approx(10));
    REQUIRE(tNodeInds1(0, 7) == Approx(9));

    //Test node indices for element 2
    REQUIRE(tNodeInds2(0, 0) == 1);
    REQUIRE(tNodeInds2(0, 1) == 2);
    REQUIRE(tNodeInds2(0, 2) == 5);
    REQUIRE(tNodeInds2(0, 3) == 4);
    REQUIRE(tNodeInds2(0, 4) == 7);
    REQUIRE(tNodeInds2(0, 5) == 8);
    REQUIRE(tNodeInds2(0, 6) == 11);
    REQUIRE(tNodeInds2(0, 7) == 10);

    //Test node ids for element 1
    REQUIRE(tNodeIds1(0, 0) == 0);
    REQUIRE(tNodeIds1(0, 1) == 1);
    REQUIRE(tNodeIds1(0, 2) == 4);
    REQUIRE(tNodeIds1(0, 3) == 3);
    REQUIRE(tNodeIds1(0, 4) == 6);
    REQUIRE(tNodeIds1(0, 5) == 7);
    REQUIRE(tNodeIds1(0, 6) == 10);
    REQUIRE(tNodeIds1(0, 7) == 9);

    //Test node ids for element 2
    REQUIRE(tNodeIds2(0, 0) == 1);
    REQUIRE(tNodeIds2(0, 1) == 2);
    REQUIRE(tNodeIds2(0, 2) == 5);
    REQUIRE(tNodeIds2(0, 3) == 4);
    REQUIRE(tNodeIds2(0, 4) == 7);
    REQUIRE(tNodeIds2(0, 5) == 8);
    REQUIRE(tNodeIds2(0, 6) == 11);
    REQUIRE(tNodeIds2(0, 7) == 10);

}




/*
 * Tests the following:
 * - Signed volume of regular subdivision
 *
 * Will test: signed surface area to make sure it's 0
 */
TEST_CASE("Regular Subdivision Geometry Check","[VOLUME_CHECK]")
{
   // Geometry Engine Setup -----------------------
    // Using a Levelset Sphere as the Geometry
//    real tRadius = 3.15;
//    real tXCenter = 5.0;
//    real tYCenter = 4.0;
//    real tZCenter = 3.0;
//    Analytic_Level_Set_Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    Gyroid<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetGyroid;

    xtk::Phase_Table<xtk::size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetGyroid,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:10x10x10";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);
    moris::Matrix<size_t,Default_Matrix_Integer> tParent = tMeshData->get_entity_connected_to_entity_loc_inds(0,EntityRank::ELEMENT,EntityRank::NODE);


    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    /*
     * Get the output mesh and write to exodus file
     */

    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCutMesh =  tXTKModel.get_cut_mesh();

    Child_Mesh_Test<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCHM = tCutMesh.get_child_mesh(129);

    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mInternalUseFlag = true;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);
    std::string tPrefix = std::getenv("XTKOUTPUT");;
    std::string tMeshOutputFile = tPrefix + "/geometry_check_regular_subdivision.e";
//
    tCutMeshData->write_output_mesh(tMeshOutputFile);


    // Compute volumes of beginning and cut mesh
    xtk::Cell<real> tVolumes;

    mesh::Mesh_Helper::compute_mesh_volume_by_buckets(*tCutMeshData,tVolumes);

    real tTotalOutputVolume = 0;
    size_t tNumBuckets = tVolumes.size();
    for(size_t iBucket = 0; iBucket<tNumBuckets; iBucket++)
    {
        tTotalOutputVolume += tVolumes(iBucket);
    }

    mesh::Mesh_Helper::compute_mesh_volume_by_buckets(*tMeshData,tVolumes);

    real tTotalInputVolume = 0;
    tNumBuckets = tVolumes.size();
    for(size_t iBucket = 0; iBucket<tNumBuckets; iBucket++)
    {
        tTotalInputVolume += tVolumes(iBucket);
    }


    CHECK(tTotalInputVolume==Approx(tTotalOutputVolume));
    /*
     * Check surface area
     */

}



TEST_CASE("Node Hierarchy Geometry Check","[REGULAR_SUBDIVISION][TEMPLATE]")
{
    // Geometry Engine Setup -----------------------
    // Using a Levelset Sphere as the Geometry
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    xtk::Phase_Table<xtk::size_t, Default_Matrix_Integer> tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);



    /*
     * Specify Mesh Inputs
     */
    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:1x1x1";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);
    moris::Matrix<size_t,Default_Matrix_Integer> tParent = tMeshData->get_entity_connected_to_entity_loc_inds(0,EntityRank::ELEMENT,EntityRank::NODE);




    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
    tXTKModel.decompose(tDecompositionMethods);

    /*
     * Get the output mesh and write to exodus file
     */

    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCutMesh =  tXTKModel.get_cut_mesh();


    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mInternalUseFlag = true;
    tOutputOptions.mAddSideSets = true;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder, tOutputOptions);

//    mesh::Mesh_Helper::write_mesh_to_template(tMeshData,"indices.inc","offsets.inc");
    std::string tPrefix = std::getenv("XTKOUTPUT");;
    std::string tMeshOutputFile = tPrefix + "/geometry_check_regular_subdivision.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);

    size_t tNumElements = tCutMeshData->get_num_entities(EntityRank::ELEMENT);

    real tTotalVolume = 0;
    for(size_t i = 0; i<24; i++)
    {
        moris::Matrix<size_t,Default_Matrix_Integer> tElementToNodes = tCutMeshData->get_entity_connected_to_entity_loc_inds(i,EntityRank::ELEMENT,EntityRank::NODE);
        moris::Matrix<real, Default_Matrix_Real> tCoordinates = tCutMeshData->get_selected_node_coordinates_loc_inds(tElementToNodes);

        real tVolume = xtk::vol_tetrahedron(tCoordinates);
        tTotalVolume +=tVolume;
    }

    CHECK(tTotalVolume == Approx(1.0));

    /*
     * Check surface area
     */


    // Compute volume
    xtk::Cell<real> tVolumes;

    mesh::Mesh_Helper::compute_mesh_volume_by_buckets(*tCutMeshData,tVolumes);

    real tTotalSTKVolume = 0;
    size_t tNumBuckets = tVolumes.size();
    for(size_t iBucket = 0; iBucket<tNumBuckets; iBucket++)
    {
        tTotalSTKVolume += tVolumes(iBucket);
    }

    CHECK(tTotalSTKVolume == Approx(1.0));

}


TEST_CASE("Regular Subdivision Base Data","[BASE_REG_SUB]")
{
    bool tOn = false;
    if(tOn)
    {
    /*
     * Specify Mesh Inputs
     */
    // Spatial Dimension
    xtk::size_t tSpatialDimension = 3;

    // Intialize STK Mesh Builder
    mesh::Mesh_Builder_Stk<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;

    moris::Matrix<size_t,Default_Matrix_Integer> tTetElementConnectivity(
            { {0,  8,  1,  14},
              {1,  8,  5,  14},
              {4,  5,  8,  14},
              {0,  4,  8,  14},
              {1,  9,  2,  14},
              {2,  9,  6,  14},
              {5,  6,  9,  14},
              {1,  5,  9,  14},
              {2,  10,  3, 14},
              {2,  6,  10, 14},
              {6,  7,  10, 14},
              {3,  10,  7, 14},
              {0,  3,  11, 14},
              {3,  7,  11, 14},
              {4,  11,  7, 14},
              {0,  11,  4, 14},
              {0,  1,  12, 14},
              {1,  2,  12, 14},
              {2,  3,  12, 14},
              {0,  12,  3, 14},
              {4,  13,  5, 14},
              {5,  13,  6, 14},
              {6,  13,  7, 14},
              {4,  7,  13, 14}});

    /*
     * Node Coordinates (Index corresponds to coordinate indexed in local to global map)
     */
    moris::Matrix<real, Default_Matrix_Real>  tNodeCoords(15,3,10);
    (tNodeCoords)(0 ,0) = 1.0;  (tNodeCoords)(0 ,1) = 0.0;  (tNodeCoords)(0 ,2) = 0.0;
    (tNodeCoords)(1 ,0) = 1.0;  (tNodeCoords)(1 ,1) = 1.0;  (tNodeCoords)(1 ,2) = 0.0;
    (tNodeCoords)(2 ,0) = 0.0;  (tNodeCoords)(2 ,1) = 1.0;  (tNodeCoords)(2 ,2) = 0.0;
    (tNodeCoords)(3 ,0) = 0.0;  (tNodeCoords)(3 ,1) = 0.0;  (tNodeCoords)(3 ,2) = 0.0;
    (tNodeCoords)(4 ,0) = 1.0;  (tNodeCoords)(4 ,1) = 0.0;  (tNodeCoords)(4 ,2) = 1.0;
    (tNodeCoords)(5 ,0) = 1.0;  (tNodeCoords)(5 ,1) = 1.0;  (tNodeCoords)(5 ,2) = 1.0;
    (tNodeCoords)(6 ,0) = 0.0;  (tNodeCoords)(6 ,1) = 1.0;  (tNodeCoords)(6 ,2) = 1.0;
    (tNodeCoords)(7 ,0) = 0.0;  (tNodeCoords)(7 ,1) = 0.0;  (tNodeCoords)(7 ,2) = 1.0;

    (tNodeCoords)(8 ,0) = 1.0;  (tNodeCoords)(8 ,1) = 0.5;  (tNodeCoords)(8 ,2) = 0.5;
    (tNodeCoords)(9 ,0) = 0.5;  (tNodeCoords)(9 ,1) = 1.0;  (tNodeCoords)(9 ,2) = 0.5;
    (tNodeCoords)(10,0) = 0.0;  (tNodeCoords)(10,1) = 0.5;  (tNodeCoords)(10,2) = 0.5;
    (tNodeCoords)(11,0) = 0.5;  (tNodeCoords)(11,1) = 0.0;  (tNodeCoords)(11,2) = 0.5;
    (tNodeCoords)(12,0) = 0.5;  (tNodeCoords)(12,1) = 0.5;  (tNodeCoords)(12,2) = 0.0;
    (tNodeCoords)(13,0) = 0.5;  (tNodeCoords)(13,1) = 0.5;  (tNodeCoords)(13,2) = 1.0;
    (tNodeCoords)(14,0) = 0.5;  (tNodeCoords)(14,1) = 0.5;  (tNodeCoords)(14,2) = 0.5;

    moris::Matrix<real, Default_Matrix_Real> tCoords(tTetElementConnectivity.n_cols(),3);
    moris::Matrix<real, Default_Matrix_Real> tCoordRow(1,3);

    real tTotalChildVol = 0;
    real tChildVol;
    size_t k = 0;

    for(size_t c = 0; c<4; c++)
    {
        for (size_t r = 0; r<24; r++)
        {
           (tTetElementConnectivity)(r,c) = (tTetElementConnectivity)(r,c) +1;
        }
    }

    xtk::Cell<moris::Matrix<size_t,Default_Matrix_Integer>> tConnectivity({tTetElementConnectivity});

    xtk::print(tTetElementConnectivity,"Tets");

    /*
     * Using different ordering than a typical ordinal of hex 8 to check robustness
     */
    moris::Matrix<size_t,Default_Matrix_Integer> tNodeLocaltoGlobal({{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}});



    /*
     * Part Names
     */
     xtk::Cell< xtk::Cell<std::string> >  tPartNames(1);
     tPartNames.push_back({"block_1"});


     /*
      * Construct Mesh From Data
      */
//     std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_data( tSpatialDimension, tConnectivity, tNodeCoords, tNodeLocaltoGlobal, tPartNames, true);
//
//     std::string tMeshOutputFile = "../TestExoFiles/Outputs/base_reg_sub.e";

//     tMeshData->write_output_mesh(tMeshOutputFile);

    }
}

}
