/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Cluster.cpp
 *
 */

#include "catch.hpp"
#include "cl_MTK_Cell_Cluster_Proxy.hpp"
#include "cl_MTK_Cell_Cluster_Input.hpp"
#include "cl_MTK_Vertex_Proxy.hpp"
#include "cl_MTK_Cell_Proxy.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

namespace moris
{

namespace mtk
{

TEST_CASE("Cell Cluster Proxy","[MTK_CLUSTER_PROXY]")
        {

    // nodes and node coordinates (used for both meshes)
    Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};
    Matrix<DDRMat> tNodeCoordinates ={{0, 0, 0},{1, 0, 0},{0, 1, 0},{1, 1, 0},{0, 0, 1},{1, 0, 1},{0, 1, 1},{1, 1, 1},{0, 0, 2},{1, 0, 2},{0, 1, 2},{1, 1, 2},{0, 0, 3},{1, 0, 3},{0, 1, 3},{1, 1, 3},{0, 0, 4},{1, 0, 4},{0, 1, 4},{1, 1, 4},{0.5, 0, 3.5},{1, 0.5, 3.5},{0.5, 1, 3.5},{0, 0.5, 3.5},{0.5, 0.5, 3},{0.5, 0.5, 4},{0.5, 0.5, 3.5},{0.1, 0, 3.1},{0.9, 0, 3.1},{0.1, 0.1, 3.1},{0.9, 0.1, 3.1},{1, 0, 3.1},{0, 0, 3.1},{1, 0.1, 3.1},{1, 0.9, 3.1},{0.9, 0.9, 3.1},{1, 1, 3.1},{0.9, 1, 3.1},{0.1, 1, 3.1},{0.1, 0.9, 3.1},{0, 1, 3.1},{0, 0.9, 3.1},{0, 0.1, 3.1},{0.5, 0.5, 3.1}};

    // setup interpolation mesh
    // interpolation cell information
    Matrix<IndexMat> tInterpCellIds     = {{4}};
    Matrix<IndexMat> tInterpCellVertices = {{13, 14, 16, 15, 17, 18, 20, 19}};

    // interpolation vertex setup
    moris::Cell<mtk::Vertex_Proxy> tInterpVertices(8);
    for(moris::uint  i = 0; i<tInterpVertices.size(); i++)
    {
        tInterpVertices(i).mVertexId  = tInterpCellVertices(i);
        tInterpVertices(i).mVertexInd = i;
        tInterpVertices(i).mVertexCoord = tNodeCoordinates.get_row(tInterpCellVertices(i)-1);
    }

    // interpolation cell setup
    mtk::Cell_Proxy tInterpCell;
    tInterpCell.mId    = tInterpCellIds(0);
    tInterpCell.mIndex = 3;
    for(moris::uint  i = 0; i<tInterpVertices.size(); i++)
    {
        tInterpCell.mVertices.push_back(&tInterpVertices(i));
    }

    // setup integration mesh
    // Cells and cell topology in material phase 0
    // Tetrathedral cells in material phase 1
    Matrix<IndexMat> tCellIdsPhase0    = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};
    Matrix<IndexMat> tCellToNodePhase0 = {{29, 14, 31, 32},{13, 28, 30, 33},{35, 16, 36, 37},{14, 31, 32, 34},{16, 36, 37, 38},{39, 15, 40, 41},{13, 30, 40, 43},{15, 40, 41, 42},{30, 13, 33, 43},{13, 14, 31, 29},{28, 13, 31, 29},{13, 28, 31, 30},{16, 14, 35, 36},{14, 31, 34, 36},{14, 34, 35, 36},{15, 16, 39, 40},{16, 36, 38, 40},{16, 38, 39, 40},{15, 13, 40, 43},{15, 40, 42, 43},{13, 30, 31, 44},{14, 13, 31, 44},{25, 13, 14, 44},{14, 31, 36, 44},{16, 14, 36, 44},{25, 14, 16, 44},{16, 36, 40, 44},{15, 16, 40, 44},{25, 16, 15, 44},{30, 13, 40, 44},{13, 15, 40, 44},{13, 25, 15, 44}};

    // Tetrathedral cells in material phase 0 void
    Matrix<IndexMat> tCellToNodeGhost0 = {{21, 27, 31, 30},{17, 18, 21, 27},{31, 27, 34, 36},{18, 20, 22, 27},{36, 27, 38, 40},{20, 19, 23, 27},{17, 24, 19, 27},{30, 27, 31, 44},{31, 27, 36, 44},{36, 27, 40, 44},{27, 30, 40, 44},{17, 26, 18, 27},{18, 26, 20, 27},{20, 26, 19, 27},{17, 19, 26, 27},{21, 28, 30, 31},{21, 28, 31, 29},{21, 29, 31, 32},{27, 21, 31, 32},{18, 21, 27, 32},{28, 21, 30, 33},{21, 27, 30, 33},{21, 17, 27, 33},{27, 22, 34, 36},{34, 22, 35, 36},{22, 35, 36, 37},{27, 22, 36, 37},{20, 22, 27, 37},{31, 27, 32, 34},{27, 18, 32, 34},{27, 22, 18, 34},{27, 23, 38, 40},{38, 23, 39, 40},{36, 27, 37, 38},{27, 20, 37, 38},{27, 23, 20, 38},{23, 39, 40, 41},{27, 23, 40, 41},{19, 23, 27, 41},{27, 24, 42, 43},{30, 27, 40, 43},{40, 27, 42, 43},{40, 27, 41, 42},{27, 19, 41, 42},{27, 24, 19, 42},{27, 30, 33, 43},{17, 27, 33, 43},{24, 27, 17, 43}};
    Matrix<IndexMat> tCellIdsGhost0    = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};

    // local coordinates
    // element level parameter
    Matrix<IndexMat> tVertexIDsInCluster          = {{13, 14, 16, 15, 17, 18, 20, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};
    Matrix<DDRMat> tLocalCoordinatesWrtInterpCell = {{-1, -1, -1},{ 1, -1, -1},{ 1, 1, -1},{-1, 1, -1},{-1, -1, 1},{ 1, -1, 1},{ 1, 1, 1},{-1, 1, 1},{ 0, -1, 0},{ 1, 0, 0},{ 0, 1, 0},{ -1, 0, 0},{0, 0, -1},{0, 0, 1},{0, 0, 0},{-0.8, -1, -0.8},{0.8, -1, -0.8},{-0.8, -0.8, -0.8},{0.8, -0.8, -0.8},{1, -1, -0.8},{-1, -1, -0.8},{1, -0.8, -0.8},{1, 0.8, -0.8},{0.8, 0.8, -0.8},{1, 1, -0.8},{0.8, 1, -0.8},{-0.8, 1, -0.8},{-0.8, 0.8, -0.8},{-1, 1, -0.8},{-1, 0.8, -0.8},{-1, -0.8, -0.8},{0, 0, -0.8}};

    // setup vertices
    uint tNumVertices = tLocalToGlobalNodeMap.numel();
    moris::Cell<Vertex_Proxy> tVertices(tNumVertices);

    for(moris::uint  i = 0; i < tNumVertices; i++)
    {
        tVertices(i).mVertexId    = tLocalToGlobalNodeMap(i);
        tVertices(i).mVertexInd   = (moris_index) i;
        tVertices(i).mVertexCoord = tNodeCoordinates.get_row(i);
    }

    // setup cells and cell cluster
    moris::Cell<Cell_Proxy> tPrimaryCells(tCellIdsPhase0.numel());
    Cell_Info_Tet4 tTet4CellInfo;
    Cell_Info_Hex8 tHex8CellInfo;

    Cell_Cluster_Proxy tCellCluster;

    moris_index  tCount = 0;
    // for cells in the material
    for(moris::uint  i = 0; i < tCellIdsPhase0.numel(); i++)
    {
        tPrimaryCells(i).mId = tCellIdsPhase0(i);
        tPrimaryCells(i).mIndex = tCount;
        tPrimaryCells(i).mCellInfo = &tTet4CellInfo;
        // add vertices
        for(moris::uint j = 0; j< tCellToNodeGhost0.n_cols(); j++)
        {
            moris_index tVertId = tCellToNodePhase0(i,j);
            tPrimaryCells( i).mVertices.push_back(&tVertices(tVertId-1));
        }

        tCellCluster.mPrimaryIntegrationCells.push_back(&tPrimaryCells(i));
        tCount++;
    }

    // for cells in the ghost
    moris::Cell<Cell_Proxy> tVoidCells( tCellIdsGhost0.numel() );

    for(moris::uint  i = 0; i < tCellIdsGhost0.numel(); i++)
    {
        tVoidCells(i).mId = tCellIdsGhost0(i);
        tVoidCells(i).mIndex = tCount;
        tVoidCells(i).mCellInfo = &tTet4CellInfo;

        // add vertices
        for(moris::uint j = 0; j< tCellToNodeGhost0.n_cols(); j++)
        {
            moris_index tVertId = tCellToNodeGhost0(i,j);
            tVoidCells(i).mVertices.push_back(&tVertices(tVertId-1));
        }
        tCellCluster.mVoidIntegrationCells.push_back(&tVoidCells(i));
        tCount++;
    }

    // set interpolation cell
    tCellCluster.mInterpolationCell = &tInterpCell;
    CHECK( tInterpCell.mId == 4 );
    CHECK( tInterpCell.mIndex == 3 );

    // set vertex in cluster
    moris::Cell<moris::mtk::Vertex const *> tVerticesInCluster(tVertexIDsInCluster.numel());
    for(moris::uint  i = 0; i <tVertexIDsInCluster.numel(); i++)
    {
        tVerticesInCluster(i) = &tVertices(tVertexIDsInCluster(i)-1);
    }

    tCellCluster.mVerticesInCluster = tVerticesInCluster;

    // set parametric coordinate relative to the interpolation cell
    tCellCluster.mVertexParamCoords = tLocalCoordinatesWrtInterpCell;
    moris::Cell<moris::mtk::Cell  const *> const & tPrimaryCells2 = tCellCluster.get_primary_cells_in_cluster();

    for(moris::uint  i = 0; i <tPrimaryCells2.size(); i++)
    {
        CHECK(tPrimaryCells2(i)->get_id() == tCellIdsPhase0(i));
    }

    moris::Cell<moris::mtk::Cell const *> const & tVoidCells2    = tCellCluster.get_void_cells_in_cluster();

    for(moris::uint  i = 0; i <tVoidCells2.size(); i++)
    {
        CHECK(tVoidCells2(i)->get_id() == tCellIdsGhost0(i));
    }

    CHECK(all_true(tLocalCoordinatesWrtInterpCell == tCellCluster.get_vertices_local_coordinates_wrt_interp_cell()));

    moris::Cell<moris::mtk::Vertex const *> const & tVerticesInCluster2 = tCellCluster.get_vertices_in_cluster();
    for(moris::uint  i = 0; i <tVerticesInCluster2.size(); i++)
    {
        CHECK(tVerticesInCluster2(i)->get_id() == tVertexIDsInCluster(i));
    }

    CHECK(tCellCluster.compute_cluster_cell_measure(Primary_Void::PRIMARY,Leader_Follower::LEADER));
    CHECK(tCellCluster.compute_cluster_cell_measure(Primary_Void::VOID,Leader_Follower::LEADER));
        }

TEST_CASE(" Same Interpolation and Integration Mesh + Cluster Input ","[MTK_MESH_CLUSTER]")
{

    uint p_size = moris::par_size();

    if( p_size == 1 ) // specify it is a serial test only
    {
        // setup the interpolation mesh
        std::string tInterpString = "generated:1x1x4";
        Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tInterpString );

        // setup the integration mesh
        // Define the Integration Mesh (from data from xtk)
        Matrix<DDRMat> tNodeCoordinates ={{0, 0, 0},
                                          {1, 0, 0},
                                          {0, 1, 0},
                                          {1, 1, 0},
                                          {0, 0, 1},
                                          {1, 0, 1},
                                          {0, 1, 1},
                                          {1, 1, 1},
                                        {0, 0, 2},
                                        {1, 0, 2},
                                        {0, 1, 2},
                                        {1, 1, 2},
                                        {0, 0, 3},
                                        {1, 0, 3},
                                        {0, 1, 3},
                                        {1, 1, 3},
                                        {0, 0, 4},
                                        {1, 0, 4},
                                        {0, 1, 4},
                                        {1, 1, 4},
                                        {0.5, 0, 3.5},
                                        {1, 0.5, 3.5},
                                        {0.5, 1, 3.5},
                                        {0, 0.5, 3.5},
                                        {0.5, 0.5, 3},
                                        {0.5, 0.5, 4},
                                        {0.5, 0.5, 3.5},
                                        {0.1, 0, 3.1},
                                        {0.9, 0, 3.1},
                                        {0.1, 0.1, 3.1},
                                        {0.9, 0.1, 3.1},
                                        {1, 0, 3.1},
                                        {0, 0, 3.1},
                                        {1, 0.1, 3.1},
                                        {1, 0.9, 3.1},
                                        {0.9, 0.9, 3.1},
                                        {1, 1, 3.1},
                                        {0.9, 1, 3.1},
                                        {0.1, 1, 3.1},
                                        {0.1, 0.9, 3.1},
                                        {0, 1, 3.1},
                                        {0, 0.9, 3.1},
                                        {0, 0.1, 3.1},
                                        {0.5, 0.5, 3.1}};

        Matrix<IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};

        Matrix<IndexMat> tInterpElemsAsIntegCellIds     = {{1,2,3,4}};
        Matrix<IndexMat> tInterpElemsAsIntegCellToNodes = {{1, 2, 4, 3, 5, 6, 8, 7},
                                                           {5, 6, 8, 7, 9, 10, 12, 11},
                                                           {9, 10, 12, 11, 13, 14, 16, 15},
                                                           {13, 14, 16, 15, 17, 18, 20, 19}};

        // Tetrathedral cells in material phase 1
        Matrix<IndexMat> tCellIdsPhase0    = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};
        Matrix<IndexMat> tCellToNodePhase0 = {{29, 14, 31, 32},
                                              {13, 28, 30, 33},
                                              {35, 16, 36, 37},
                                              {14, 31, 32, 34},
                                              {16, 36, 37, 38},
                                              {39, 15, 40, 41},
                                              {13, 30, 40, 43},
                                              {15, 40, 41, 42},
                                              {30, 13, 33, 43},
                                              {13, 14, 31, 29},
                                              {28, 13, 31, 29},
                                              {13, 28, 31, 30},
                                              {16, 14, 35, 36},
                                              {14, 31, 34, 36},
                                              {14, 34, 35, 36},
                                              {15, 16, 39, 40},
                                              {16, 36, 38, 40},
                                              {16, 38, 39, 40},
                                              {15, 13, 40, 43},
                                              {15, 40, 42, 43},
                                              {13, 30, 31, 44},
                                              {14, 13, 31, 44},
                                              {25, 13, 14, 44},
                                              {14, 31, 36, 44},
                                              {16, 14, 36, 44},
                                              {25, 14, 16, 44},
                                              {16, 36, 40, 44},
                                              {15, 16, 40, 44},
                                              {25, 16, 15, 44},
                                              {30, 13, 40, 44},
                                              {13, 15, 40, 44},
                                              {13, 25, 15, 44}};

        // Tetrathedral cells in material phase 1
        Matrix<IndexMat> tCellToNodeGhost0 = {{21, 27, 31, 30},
                                                    {17, 18, 21, 27},
                                                    {31, 27, 34, 36},
                                                    {18, 20, 22, 27},
                                                    {36, 27, 38, 40},
                                                    {20, 19, 23, 27},
                                                    {17, 24, 19, 27},
                                                    {30, 27, 31, 44},
                                                    {31, 27, 36, 44},
                                                    {36, 27, 40, 44},
                                                    {27, 30, 40, 44},
                                                    {17, 26, 18, 27},
                                                    {18, 26, 20, 27},
                                                    {20, 26, 19, 27},
                                                    {17, 19, 26, 27},
                                                    {21, 28, 30, 31},
                                                    {21, 28, 31, 29},
                                                    {21, 29, 31, 32},
                                                    {27, 21, 31, 32},
                                                    {18, 21, 27, 32},
                                                    {28, 21, 30, 33},
                                                    {21, 27, 30, 33},
                                                    {21, 17, 27, 33},
                                                    {27, 22, 34, 36},
                                                    {34, 22, 35, 36},
                                                    {22, 35, 36, 37},
                                                    {27, 22, 36, 37},
                                                    {20, 22, 27, 37},
                                                    {31, 27, 32, 34},
                                                    {27, 18, 32, 34},
                                                    {27, 22, 18, 34},
                                                    {27, 23, 38, 40},
                                                    {38, 23, 39, 40},
                                                    {36, 27, 37, 38},
                                                    {27, 20, 37, 38},
                                                    {27, 23, 20, 38},
                                                    {23, 39, 40, 41},
                                                    {27, 23, 40, 41},
                                                    {19, 23, 27, 41},
                                                    {27, 24, 42, 43},
                                                    {30, 27, 40, 43},
                                                    {40, 27, 42, 43},
                                                    {40, 27, 41, 42},
                                                    {27, 19, 41, 42},
                                                    {27, 24, 19, 42},
                                                    {27, 30, 33, 43},
                                                    {17, 27, 33, 43},
                                                    {24, 27, 17, 43}};

        Matrix<IndexMat> tCellIdsGhost0 = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};

        moris::mtk::MtkSetsInfo tMtkMeshSets;
        // Define side sets on the integration mesh (i.e. fixed bc, interface and ghost)
        // Fixed bc
        moris::mtk::MtkSideSetInfo tFixed;
        Matrix<IndexMat> tFixedElementsAndOrds = {{1, 4}};
        tFixed.mElemIdsAndSideOrds = &tFixedElementsAndOrds;
        tFixed.mSideSetName        = "fixed" ;
        tMtkMeshSets.add_side_set(&tFixed);

        // Fixed bc
        moris::mtk::MtkSideSetInfo tGhost;
        Matrix<IndexMat> tGhostCellAndOrds = {{3,5},{4,4}};
        tGhost.mElemIdsAndSideOrds = &tGhostCellAndOrds;
        tGhost.mSideSetName        = "ghost_facets" ;
        tMtkMeshSets.add_side_set(&tGhost);

        // add block sets
        // Tet Cells in Omega 0
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
        tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
        tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
        tOmega0BlockSetTet.mBlockSetTopo = CellTopology::TET4;
        tMtkMeshSets.add_block_set(&tOmega0BlockSetTet);

        // Hex Cells in Omega 0
        Matrix<IdMat> tOmega0HexCellIds = {{1,2,3}};
        moris::mtk::MtkBlockSetInfo tOmega0BlockSetHex;
        tOmega0BlockSetHex.mCellIdsInSet = &tOmega0HexCellIds;
        tOmega0BlockSetHex.mBlockSetName = "Omega_0_hex";
        tOmega0BlockSetHex.mBlockSetTopo = CellTopology::HEX8;
        tMtkMeshSets.add_block_set(&tOmega0BlockSetHex);

        // Cells in the ghost domain of omega 1
        moris::mtk::MtkBlockSetInfo tOmega0GhostBlockSetTet;
        tOmega0GhostBlockSetTet.mCellIdsInSet = &tCellIdsGhost0;
        tOmega0GhostBlockSetTet.mBlockSetName = "Omega_0_Ghost";
        tOmega0GhostBlockSetTet.mBlockSetTopo = CellTopology::TET4;
        tMtkMeshSets.add_block_set(&tOmega0GhostBlockSetTet);

        // Integration Cells for Ghost penalization only
        Matrix<IdMat> tGhostCellIds = {{3,4}};
        moris::mtk::MtkBlockSetInfo tCellsForGhost;
        tCellsForGhost.mCellIdsInSet = &tGhostCellIds;
        tCellsForGhost.mBlockSetName = "Ghost_Cells_0";
        tCellsForGhost.mBlockSetTopo = CellTopology::HEX8;
        tMtkMeshSets.add_block_set(&tCellsForGhost);

        // Mesh data input structure
        moris::mtk::MtkMeshData tMeshDataInput(3);

        moris::uint tSpatialDim   = 3;
        Matrix<IdMat> tNodeOwner(1,tNodeCoordinates.n_rows(),moris::par_rank());
        tMeshDataInput.ElemConn(0)             = &tInterpElemsAsIntegCellToNodes;
        tMeshDataInput.ElemConn(1)             = &tCellToNodePhase0;
        tMeshDataInput.ElemConn(2)             = &tCellToNodeGhost0;
        tMeshDataInput.LocaltoGlobalElemMap(0) = (&tInterpElemsAsIntegCellIds);
        tMeshDataInput.LocaltoGlobalElemMap(1) = (&tCellIdsPhase0);
        tMeshDataInput.LocaltoGlobalElemMap(2) = (&tCellIdsGhost0);
        tMeshDataInput.CreateAllEdgesAndFaces  = true;
        tMeshDataInput.Verbose                 = false;
        tMeshDataInput.SpatialDim              = &tSpatialDim;
        tMeshDataInput.NodeCoords              = &tNodeCoordinates;
        tMeshDataInput.NodeProcOwner           = &tNodeOwner;
        tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
        tMeshDataInput.SetsInfo                = &tMtkMeshSets;
        tMeshDataInput.MarkNoBlockForIO        = false;

        // Get the 4th cell from interpolation mesh
        moris::moris_id tInterpCellIndex = 3;
        moris::mtk::Cell* tInterpCell = &tInterpMesh1->get_mtk_cell(tInterpCellIndex);

        // setup integration mesh
        // Cells and cell topology in material phase 0
        // Tetrathedral cells in material phase 1
        Matrix<IndexMat> tCellIdsCluster1Material = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};

        // Tetrathedral cells in material phase 0 void
        Matrix<IndexMat> tCellIdsCluster1Void = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};

        // local coordinates
        // element level parameter
        Matrix<IdMat> tVertexIDsInCluster  = {{13, 14, 16, 15, 17, 18, 20, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};
        Matrix<DDRMat> tLocalCoordinatesWrtInterpCell = {{-1, -1, -1},{ 1, -1, -1},{ 1, 1, -1},{-1, 1, -1},{-1, -1, 1},{ 1, -1, 1},{ 1, 1, 1},{-1, 1, 1},{ 0, -1, 0},{ 1, 0, 0},{ 0, 1, 0},{ -1, 0, 0},{0, 0, -1},{0, 0, 1},{0, 0, 0},{-0.8, -1, -0.8},{0.8, -1, -0.8},{-0.8, -0.8, -0.8},{0.8, -0.8, -0.8},{1, -1, -0.8},{-1, -1, -0.8},{1, -0.8, -0.8},{1, 0.8, -0.8},{0.8, 0.8, -0.8},{1, 1, -0.8},{0.8, 1, -0.8},{-0.8, 1, -0.8},{-0.8, 0.8, -0.8},{-1, 1, -0.8},{-1, 0.8, -0.8},{-1, -0.8, -0.8},{0, 0, -0.8}};

        // setup cluster data
        Cell_Cluster_Input tCellClusterInput;
        tCellClusterInput.add_cluster_data(tInterpCell,&tCellIdsCluster1Material,&tCellIdsCluster1Void,&tVertexIDsInCluster,&tLocalCoordinatesWrtInterpCell);

        tMeshDataInput.CellClusterInput = & tCellClusterInput;

        Integration_Mesh* tIntegMesh1  = create_integration_mesh(MeshType::STK,tMeshDataInput,tInterpMesh1);

        // check the integration mesh cluster of interpolation cell id 4
        Cell_Cluster const & tCellClusterIndex3 = tIntegMesh1->get_cell_cluster(*tInterpCell);

        // Check Ids in cluster
        moris::Matrix<moris::IdMat> tPrimaryCellIdsInCluster = tCellClusterIndex3.get_primary_cell_ids_in_cluster();
        moris::Matrix<moris::IdMat> tVoidCellIdsInCluster = tCellClusterIndex3.get_void_cell_ids_in_cluster();

        CHECK(all_true(tPrimaryCellIdsInCluster == tCellIdsCluster1Material));
        CHECK(all_true(tVoidCellIdsInCluster == tCellIdsCluster1Void));

        // Check Indices in cluster
        moris::Matrix<moris::IndexMat> tPrimaryCellIndsInCluster = tCellClusterIndex3.get_primary_cell_indices_in_cluster();
        moris::Matrix<moris::IndexMat> tVoidCellIndsInCluster = tCellClusterIndex3.get_void_cell_indices_in_cluster();

        // convert to ids
        moris::Matrix<moris::IdMat> tConvertedPrimaryCellIds = convert_entity_indices_to_ids(tPrimaryCellIndsInCluster, EntityRank::ELEMENT, tIntegMesh1);
        moris::Matrix<moris::IdMat> tConvertedVoidCellIds = convert_entity_indices_to_ids(tVoidCellIndsInCluster, EntityRank::ELEMENT, tIntegMesh1);

        CHECK(all_true(tPrimaryCellIdsInCluster == tConvertedPrimaryCellIds));
        CHECK(all_true(tVoidCellIdsInCluster == tConvertedVoidCellIds));

        // Check vertices
        moris::Matrix<moris::IdMat> tVertexIds = tCellClusterIndex3.get_vertex_ids_in_cluster();
        CHECK(all_true(tVertexIDsInCluster == tVertexIds));

        moris::Matrix<moris::IndexMat> tVertexInds = tCellClusterIndex3.get_vertex_indices_in_cluster();
        moris::Matrix<moris::IdMat> tConvertedVertexIds = convert_entity_indices_to_ids(tVertexInds, EntityRank::NODE, tIntegMesh1);
        CHECK(all_true(tConvertedVertexIds == tVertexIds));

        // Check local coordinates wrt interpolation cell
        moris::Matrix<moris::DDRMat> tLocalCoords = tCellClusterIndex3.get_vertices_local_coordinates_wrt_interp_cell();

        CHECK(all_true(tLocalCoordinatesWrtInterpCell == tLocalCoords));

        // check the local coordinates we receive from primary cells
        moris::Cell<moris::mtk::Cell const *> const & tPrimaryCells = tCellClusterIndex3.get_primary_cells_in_cluster();
        for(moris::uint i = 0; i < tPrimaryCells.size(); i++)
        {
            moris::Matrix<moris::DDRMat> tCellParametricCoords = tCellClusterIndex3.get_primary_cell_local_coords_on_side_wrt_interp_cell(i);

            moris::Cell<moris::mtk::Vertex*> tVertsOnSide = tPrimaryCells(i)->get_vertex_pointers();

            moris::Matrix<moris::DDRMat> tGoldParamCoords(4,3);
            for(moris::uint j = 0; j<tVertsOnSide.size(); j++)
            {
                tGoldParamCoords.get_row(j) = tCellClusterIndex3.get_vertex_local_coordinate_wrt_interp_cell(tVertsOnSide(j)).get_row(0);
            }

            CHECK(all_true(tGoldParamCoords == tCellParametricCoords));

        }

        // check the local coordinates we receive from void cells
        moris::Cell<moris::mtk::Cell const *> const & tVoidCells = tCellClusterIndex3.get_void_cells_in_cluster();

        for(moris::uint i = 0; i < tVoidCells.size(); i++)
        {
            moris::Matrix<moris::DDRMat> tCellParametricCoords = tCellClusterIndex3.get_void_cell_local_coords_on_side_wrt_interp_cell(i);

            moris::Cell<moris::mtk::Vertex*> tVertsOnSide = tVoidCells(i)->get_vertex_pointers();

            moris::Matrix<moris::DDRMat> tGoldParamCoords(4,3);
            for(moris::uint j = 0; j<tVertsOnSide.size(); j++)
            {
                tGoldParamCoords.get_row(j) = tCellClusterIndex3.get_vertex_local_coordinate_wrt_interp_cell(tVertsOnSide(j)).get_row(0);
            }

            CHECK(all_true(tGoldParamCoords == tCellParametricCoords));

        }

        // check the integration mesh cluster of interpolation cell index 0
        moris_index tInterpCellIndex0 = 0;
        moris_index tInterpCellId0 = 1;
        mtk::Cell const & tInterpCell0   = tInterpMesh1->get_mtk_cell(tInterpCellIndex0);
        Cell_Cluster const & tCellCluster0 = tIntegMesh1->get_cell_cluster(tInterpCell0);

        CHECK(tCellCluster0.is_trivial());

        // Check Ids in cluster
        tPrimaryCellIdsInCluster = tCellCluster0.get_primary_cell_ids_in_cluster();

        CHECK(tInterpCellId0 == tPrimaryCellIdsInCluster(0));

        // Check Indices in cluster
        tPrimaryCellIndsInCluster = tCellCluster0.get_primary_cell_indices_in_cluster();

        // convert to ids
        tConvertedPrimaryCellIds = convert_entity_indices_to_ids(tPrimaryCellIndsInCluster, EntityRank::ELEMENT, tIntegMesh1);

        CHECK(tInterpCellId0 == tConvertedPrimaryCellIds(0));

        // test block set access
        moris::Cell<std::string> tBlockSetNames = tIntegMesh1->get_block_set_names();

        // the sets that have no primary integration cells do not show up here
        CHECK(tBlockSetNames.size() == 7);
        CHECK(tBlockSetNames(3).compare("Omega_0_tets")==0);
        CHECK(tBlockSetNames(4).compare("Omega_0_hex")==0);
        CHECK(tBlockSetNames(6).compare("Ghost_Cells_0")==0);

        moris::Cell<Cluster const *> tClustersInBlock0 = tIntegMesh1->get_cell_clusters_in_set(3);
        moris::Cell<Cluster const *> tClustersInBlock1 = tIntegMesh1->get_cell_clusters_in_set(4);
        moris::Cell<Cluster const *> tClustersInBlock2 = tIntegMesh1->get_cell_clusters_in_set(6);

        CHECK(tClustersInBlock0.size() == 1);
        CHECK(all_true(tClustersInBlock0(0)->get_primary_cell_ids_in_cluster() == tCellIdsCluster1Material));
        CHECK(tClustersInBlock1.size() == 3);

        // cleanup
        delete tInterpMesh1;
        delete tIntegMesh1;
    }

}
}
}

