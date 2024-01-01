/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Writer_Exodus_Test.cpp
 *
 */

#include <MTK/src/cl_MTK_Writer_Exodus.hpp>
#include <MTK/src/cl_MTK_Reader_Exodus.hpp>
#include <MTK/src/cl_MTK_Integration_Mesh.hpp>
#include "catch.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_MTK_Mesh.hpp" // MTK/src
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris
{
    namespace mtk
    {
        TEST_CASE("MTK Exodus Writer","[WRITE_EXO]")
        {
            Matrix <DDRMat> tNodeCoordinates = {{0,   0,   0},
                                                {1,   0,   0},
                                                {0,   1,   0},
                                                {1,   1,   0},
                                                {0,   0,   1},
                                                {1,   0,   1},
                                                {0,   1,   1},
                                                {1,   1,   1},
                                                {0,   0,   2},
                                                {1,   0,   2},
                                                {0,   1,   2},
                                                {1,   1,   2},
                                                {0,   0,   3},
                                                {1,   0,   3},
                                                {0,   1,   3},
                                                {1,   1,   3},
                                                {0,   0,   4},
                                                {1,   0,   4},
                                                {0,   1,   4},
                                                {1,   1,   4},
                                                {0.5, 0,   3.5},
                                                {1,   0.5, 3.5},
                                                {0.5, 1,   3.5},
                                                {0,   0.5, 3.5},
                                                {0.5, 0.5, 3},
                                                {0.5, 0.5, 4},
                                                {0.5, 0.5, 3.5},
                                                {0.1, 0,   3.1},
                                                {0.9, 0,   3.1},
                                                {0.1, 0.1, 3.1},
                                                {0.9, 0.1, 3.1},
                                                {1,   0,   3.1},
                                                {0,   0,   3.1},
                                                {1,   0.1, 3.1},
                                                {1,   0.9, 3.1},
                                                {0.9, 0.9, 3.1},
                                                {1,   1,   3.1},
                                                {0.9, 1,   3.1},
                                                {0.1, 1,   3.1},
                                                {0.1, 0.9, 3.1},
                                                {0,   1,   3.1},
                                                {0,   0.9, 3.1},
                                                {0,   0.1, 3.1},
                                                {0.5, 0.5, 3.1}};

            Matrix <IndexMat> tLocalToGlobalNodeMap = {{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44}};

            Matrix <IndexMat> tInterpElemsAsIntegCellIds = {{1, 2, 3, 4}};
            Matrix <IndexMat> tInterpElemsAsIntegCellToNodes = {{1,  2,  4,  3,  5,  6,  8,  7},
                                                                {5,  6,  8,  7,  9,  10, 12, 11},
                                                                {9,  10, 12, 11, 13, 14, 16, 15},
                                                                {13, 14, 16, 15, 17, 18, 20, 19}};

            // Tetrathedral cells in material phase 1
            Matrix <IndexMat> tCellIdsPhase0 = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};
            Matrix <IndexMat> tCellToNodePhase0 = {{29, 14, 31, 32},
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
            Matrix <IndexMat> tCellToNodeGhost0 = {{21, 27, 31, 30},
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

            Matrix <IndexMat> tCellIdsGhost0 = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};

            moris::mtk::MtkSetsInfo tMtkMeshSets;
            // Define side sets on the integration mesh (i.e. fixed bc, interface and ghost)

            // interface side set
            moris::mtk::MtkSideSetInfo tInterfaceSideSet;
            Matrix <IndexMat> tInterfaceElemIdandSideOrd = {{6,  2},
                                                            {8,  1},
                                                            {9,  2},
                                                            {10, 2},
                                                            {12, 1},
                                                            {13, 2},
                                                            {14, 1},
                                                            {16, 2},
                                                            {17, 1},
                                                            {18, 1},
                                                            {20, 2},
                                                            {21, 2},
                                                            {22, 2},
                                                            {23, 2},
                                                            {24, 1},
                                                            {29, 1},
                                                            {30, 1},
                                                            {32, 2},
                                                            {33, 1},
                                                            {34, 1},
                                                            {37, 2},
                                                            {41, 2},
                                                            {43, 1},
                                                            {44, 1},
                                                            {45, 1},
                                                            {48, 2},
                                                            {52, 2},
                                                            {54, 1},
                                                            {55, 1},
                                                            {56, 2},
                                                            {59, 1},
                                                            {63, 1},
                                                            {65, 2},
                                                            {66, 2},
                                                            {67, 2},
                                                            {70, 1},
                                                            {73, 1},
                                                            {76, 1},
                                                            {79, 1},
                                                            {82, 2}};
            tInterfaceSideSet.mElemIdsAndSideOrds = &tInterfaceElemIdandSideOrd;
            tInterfaceSideSet.mSideSetName = "iside";
            tMtkMeshSets.add_side_set(&tInterfaceSideSet);

            // Fixed bc
            moris::mtk::MtkSideSetInfo tFixed;
            Matrix <IndexMat> tFixedElementsAndOrds = {{1, 4}};
            tFixed.mElemIdsAndSideOrds = &tFixedElementsAndOrds;
            tFixed.mSideSetName = "fixed";
            tMtkMeshSets.add_side_set(&tFixed);

            // Fixed bc
            moris::mtk::MtkSideSetInfo tGhost;
            Matrix <IndexMat> tGhostCellAndOrds = {{3, 5},
                                                   {4, 4}};
            tGhost.mElemIdsAndSideOrds = &tGhostCellAndOrds;
            tGhost.mSideSetName = "ghost_facets";
            tMtkMeshSets.add_side_set(&tGhost);

            // add block sets
            // Tet Cells in Omega 0
            moris::mtk::MtkBlockSetInfo tOmega0BlockSetTet;
            tOmega0BlockSetTet.mCellIdsInSet = &tCellIdsPhase0;
            tOmega0BlockSetTet.mBlockSetName = "Omega_0_tets";
            tOmega0BlockSetTet.mBlockSetTopo = CellTopology::TET4;
            tMtkMeshSets.add_block_set(&tOmega0BlockSetTet);

            // Hex Cells in Omega 0
            Matrix <IdMat> tOmega0HexCellIds = {{1, 2, 3}};
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
            Matrix <IdMat> tGhostCellIds = {{3, 4}};
            moris::mtk::MtkBlockSetInfo tCellsForGhost;
            tCellsForGhost.mCellIdsInSet = &tGhostCellIds;
            tCellsForGhost.mBlockSetName = "Ghost_Cells_0";
            tCellsForGhost.mBlockSetTopo = CellTopology::HEX8;
            tMtkMeshSets.add_block_set(&tCellsForGhost);

            // Mesh data input structure
            moris::mtk::MtkMeshData tMeshDataInput(3);

            moris::uint tSpatialDim = 3;
            Matrix <IdMat> tNodeOwner(1, tNodeCoordinates.n_rows(), moris::par_rank());
            tMeshDataInput.ElemConn(0) = &tInterpElemsAsIntegCellToNodes;
            tMeshDataInput.ElemConn(1) = &tCellToNodePhase0;
            tMeshDataInput.ElemConn(2) = &tCellToNodeGhost0;
            tMeshDataInput.LocaltoGlobalElemMap(0) = (&tInterpElemsAsIntegCellIds);
            tMeshDataInput.LocaltoGlobalElemMap(1) = (&tCellIdsPhase0);
            tMeshDataInput.LocaltoGlobalElemMap(2) = (&tCellIdsGhost0);
            tMeshDataInput.CreateAllEdgesAndFaces = true;
            tMeshDataInput.Verbose = false;
            tMeshDataInput.SpatialDim = &tSpatialDim;
            tMeshDataInput.NodeCoords = &tNodeCoordinates;
            tMeshDataInput.NodeProcOwner = &tNodeOwner;
            tMeshDataInput.LocaltoGlobalNodeMap = &tLocalToGlobalNodeMap;
            tMeshDataInput.SetsInfo = &tMtkMeshSets;
            tMeshDataInput.MarkNoBlockForIO = false;

            Integration_Mesh *tIntegMeshData = moris::mtk::create_integration_mesh(MeshType::STK, tMeshDataInput);

            // Write mesh
            Writer_Exodus writer(tIntegMeshData);
            writer.write_mesh("", "test_write.exo", "", "test_temp.exo");

            // Nodal fields
            moris::Vector<std::string> tNodalFieldNames(3);
            tNodalFieldNames(0) = "ux";
            tNodalFieldNames(1) = "uy";
            tNodalFieldNames(2) = "uz";
            moris::Matrix<moris::DDRMat> xField(44, 1, 1);
            moris::Matrix<moris::DDRMat> yField(44, 1, 2);
            moris::Matrix<moris::DDRMat> zField(44, 1, 3);
            for (int i = 0; i < 44; i++)
            {
                xField(i) = moris::real(44 - i);
                yField(i) = moris::real(i);
                zField(i) = moris::real(i * i);
            }

            // Elemental fields
            moris::Vector<std::string> tElementalFieldNames(1);
            tElementalFieldNames(0) = "pressure";
            moris::Matrix<moris::DDRMat> tetField(32, 1, 4);
            moris::Matrix<moris::DDRMat> hexField(3, 1, 5);

            // Global Fields
            moris::Vector<std::string> tGlobalVariableNames(2);
            tGlobalVariableNames(0) = "mass";
            tGlobalVariableNames(1) = "time";
            moris::Matrix<DDRMat> tGlobalVariableValues(2,1);
            tGlobalVariableValues(0) = 99;
            tGlobalVariableValues(1) = 1;

            // Write the fields
            writer.set_nodal_fields(tNodalFieldNames);
            writer.set_elemental_fields(tElementalFieldNames);
            writer.set_global_variables(tGlobalVariableNames);
            writer.set_time(0.0);
            writer.write_nodal_field("ux", xField);
            writer.write_nodal_field("uy", yField);
            writer.write_nodal_field("uz", zField);
            writer.write_elemental_field("Omega_0_tets", "pressure", tetField);
            writer.write_elemental_field("Omega_0_hex", "pressure", hexField);
            writer.write_global_variables(tGlobalVariableNames, tGlobalVariableValues);
            writer.close_file();

            delete tIntegMeshData;
        }
    }
}

