/*
 * UT_MTK_Cell_Shape.cpp
 *
 *  Created on: Apr. 4, 2021
 *      Author: gates
 */

#include "catch.hpp"
#include "cl_Communication_Tools.hpp"

// base class
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Cell_STK.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Vertex_STK.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell_Info_Quad9.hpp"
#include "cl_MTK_Cell_Info_Quad16.hpp"
#include "cl_MTK_Cell_Info_Tri3.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Space_Interpolator.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_reindex_mat.hpp"

namespace moris
{
    namespace mtk
    {
        TEST_CASE("MTK Cell Shape Hex8 Rectangular", "[MTK],[Cell_Shape],[Shape_Hex8_Rectangular]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0, 0.0},
                                          {1.0, 0.0, 0.0},
                                          {1.0, 1.0, 0.0},
                                          {0.0, 1.0, 0.0},
                                          {0.0, 0.0, 1.0},
                                          {1.0, 0.0, 1.0},
                                          {1.0, 1.0, 1.0},
                                          {0.0, 1.0, 1.0}};
                Matrix<IdMat> tNodeLocalToGlobal = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX8;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                // std::string tFileOutput = "./mtk_hex8_cell_shape_ut.exo";
                // tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 8);
                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::RECTANGULAR);
            }
        }

        TEST_CASE("MTK Cell Shape Hex8 Straight", "[MTK],[Cell_Shape],[Shape_Hex8_Straight]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0, 0.0},
                                          {1.0, 0.0, 0.0},
                                          {1.0, 1.0, 0.0},
                                          {0.0, 1.0, 0.0},
                                          {0.0, 0.0, 1.0},
                                          {2.0, 0.0, 1.0},
                                          {2.0, 1.0, 1.0},
                                          {0.0, 1.0, 1.0}};
                Matrix<IdMat> tNodeLocalToGlobal = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX8;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                // std::string tFileOutput = "./mtk_hex8_cell_shape_ut.exo";
                // tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 8);
                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::STRAIGHT);
            }
        }

        TEST_CASE("MTK Cell Shape Hex8 General", "[MTK],[Cell_Shape],[Shape_Hex8_General]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0, 0.0},
                                          {1.0, 0.0, 0.0},
                                          {1.0, 1.0, 0.0},
                                          {0.0, 1.0, 0.0},
                                          {0.0, 0.0, 1.0},
                                          {2.0, 0.0, 1.0},
                                          {2.0, 1.0, 1.0},
                                          {0.0, 1.0, 1.5}};
                Matrix<IdMat> tNodeLocalToGlobal = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX8;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                // std::string tFileOutput = "./mtk_hex8_cell_shape_ut.exo";
                // tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 8);
                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::GENERAL);
            }
        }

        TEST_CASE("MTK Cell Shape Hex27 Rectangular", "[MTK],[Cell_Shape],[Shape_Hex27_Rectangular]")
        {
            if (par_size() <= 1)
            {
                //create a 3D MORIS mesh of hex27's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{ -1.0, -1.0, -1.0 },
                                          { +1.0, -1.0, -1.0 },
                                          { +1.0, +1.0, -1.0 },
                                          { -1.0, +1.0, -1.0 },
                                          { -1.0, -1.0, +1.5 },
                                          { +1.0, -1.0, +1.5 },
                                          { +1.0, +1.0, +1.5 },
                                          { -1.0, +1.0, +1.5 },
                                          { +0.0, -1.0, -1.0 },
                                          { +1.0,  0.0, -1.0 },
                                          {  0.0, +1.0, -1.0 },
                                          { -1.0,  0.0, -1.0 },
                                          { -1.0, -1.0,  0.0 },
                                          { +1.0, -1.0,  0.0 },
                                          { +1.0, +1.0,  0.0 },
                                          { -1.0, +1.0,  0.0 },
                                          {  0.0, -1.0, +1.5 },
                                          { +1.0,  0.0, +1.5 },
                                          {  0.0, +1.0, +1.5 },
                                          { -1.0,  0.0, +1.5 },
                                          {  0.0,  0.0,  0.0 },
                                          {  0.0,  0.0, -1.0 },
                                          {  0.0,  0.0, +1.5 },
                                          { -1.0,  0.0,  0.0 },
                                          { +1.0,  0.0,  0.0 },
                                          {  0.0, -1.0,  0.0 },
                                          {  0.0, +1.0,  0.0 }};

                Matrix<IdMat> tNodeLocalToGlobal = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26}};
                Matrix<IndexMat> tNodeIds = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX27;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                // std::string tFileOutput = "./mtk_hex27_cell_shape_ut.exo";
                // tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 27);
                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::RECTANGULAR);
            }
        }

        TEST_CASE("MTK Cell Shape Hex27 Straight", "[MTK],[Cell_Shape],[Shape_Hex27_Straight]")
        {
            if (par_size() <= 1)
            {
                //create a 3D MORIS mesh of hex27's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{ -2.0, -1.0, -1.0 },
                                          { +1.0, -1.0, -1.0 },
                                          { +1.0, +1.0, -1.0 },
                                          { -1.0, +1.0, -1.0 },
                                          { -2.0, -1.0, +1.5 },
                                          { +1.0, -1.0, +1.5 },
                                          { +1.0, +1.0, +1.5 },
                                          { -1.0, +1.0, +1.5 },
                                          { +0.0, -1.0, -1.0 },
                                          { +1.0,  0.0, -1.0 },
                                          {  0.0, +1.0, -1.0 },
                                          { -1.5,  0.0, -1.0 },
                                          { -2.0, -1.0,  0.0 },
                                          { +1.0, -1.0,  0.0 },
                                          { +1.0, +1.0,  0.0 },
                                          { -1.0, +1.0,  0.0 },
                                          {  0.0, -1.0, +1.5 },
                                          { +1.0,  0.0, +1.5 },
                                          {  0.0, +1.0, +1.5 },
                                          { -1.5,  0.0, +1.5 },
                                          {  0.0,  0.0,  0.0 },
                                          {  0.0,  0.0, -1.0 },
                                          {  0.0,  0.0, +1.5 },
                                          { -1.5,  0.0,  0.0 },
                                          { +1.0,  0.0,  0.0 },
                                          {  0.0, -1.0,  0.0 },
                                          {  0.0, +1.0,  0.0 }};

                Matrix<IdMat> tNodeLocalToGlobal = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26}};
                Matrix<IndexMat> tNodeIds = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX27;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                std::string tFileOutput = "./mtk_hex27_cell_shape_ut.exo";
                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 27);
                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::STRAIGHT);
            }
        }

        TEST_CASE("MTK Cell Shape Hex27 General", "[MTK],[Cell_Shape],[Shape_Hex27_General]")
        {
            if (par_size() <= 1)
            {
                //create a 3D MORIS mesh of hex27's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 3; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{ -2.0, -1.0, -1.0 },
                                          { +1.0, -1.0, -1.0 },
                                          { +1.0, +1.0, -1.0 },
                                          { -1.0, +1.0, -1.0 },
                                          { -2.0, -1.0, +1.5 },
                                          { +1.0, -1.0, +1.5 },
                                          { +1.0, +1.0, +1.5 },
                                          { -1.0, +1.0, +1.5 },
                                          { +0.0, -1.0, -1.0 },
                                          { +1.0,  0.0, -1.0 },
                                          {  0.0, +1.0, -1.0 },
                                          { -1.5,  0.0, -1.0 },
                                          { -2.0, -1.0,  0.0 },
                                          { +1.0, -1.0,  0.0 },
                                          { +1.0, +1.0,  0.0 },
                                          { -1.0, +1.0,  0.0 },
                                          {  0.0, -1.0, +1.5 },
                                          { +1.0,  0.0, +1.5 },
                                          {  0.0, +1.0, +1.5 },
                                          { -1.5,  0.0, +1.5 },
                                          {  0.0,  0.0,  0.0 },
                                          {  0.0,  0.0, -1.0 },
                                          {  0.0,  0.0, +1.5 },
                                          { -1.0,  0.0,  0.0 },
                                          { +1.0,  0.0,  0.0 },
                                          {  0.0, -1.0,  0.0 },
                                          {  0.0, +1.0,  0.0 }};

                Matrix<IdMat> tNodeLocalToGlobal = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                Matrix<IndexMat> tNodeIndices = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26}};
                Matrix<IndexMat> tNodeIds = {{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::HEX27;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeLocalToGlobal;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                std::string tFileOutput = "./mtk_hex27_cell_shape_ut.exo";
                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 27);
                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::GENERAL);
            }
        }

        TEST_CASE("MTK Cell Shape Quad4 Rectangular", "[MTK],[Cell_Shape],[Shape_Quad4_Rectangular]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                                          {1.0, 0.0},
                                          {1.0, 1.0},
                                          {0.0, 1.0}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD4;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                //std::string tFileOutput = "./mtk_quad4_cell_shape_ut.exo";
                //tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 4);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::RECTANGULAR);
            }
        }

        TEST_CASE("MTK Cell Shape Quad4 Parallel", "[MTK],[Cell_Shape],[Shape_Quad4_Parallel]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                        {1.0, 0.5},
                        {1.0, 1.5},
                        {0.0, 1.0}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD4;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                //std::string tFileOutput = "./mtk_quad4_cell_shape_ut.exo";
                //tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 4);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::PARALLEL);
            }
        }

        TEST_CASE("MTK Cell Shape Quad4 Parallel2", "[MTK],[Cell_Shape],[Shape_Quad4_Parallel2]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {
                        {1.0, 0.5},
                        {1.0, 1.5},
                        {0.0, 1.0},
                        {0.0, 0.0}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD4;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                //std::string tFileOutput = "./mtk_quad4_cell_shape_ut.exo";
                //tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 4);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::PARALLEL);
            }
        }

        TEST_CASE("MTK Cell Shape Quad4 Straight", "[MTK],[Cell_Shape],[Shape_Quad4_Straight]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                                          {1.0, 0.0},
                                          {1.0, 1.0},
                                          {0.0, 1.5}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD4;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
                //std::string tFileOutput = "./mtk_quad4_cell_shape_ut.exo";
                //tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 4);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::STRAIGHT);
            }
        }

        TEST_CASE("MTK Cell Shape Quad9 Rectangular", "[MTK],[Cell_Shape],[Shape_Quad9_Rectangular]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                                          {1.0, 0.0},
                                          {1.0, 1.5},
                                          {0.0, 1.5},
                                          {0.5, 0.0},
                                          {1.0, 0.5},
                                          {0.5, 1.5},
                                          {0.0, 0.5},
                                          {0.5, 0.5}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD9;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeIds;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
//                std::string tFileOutput = "./mtk_quad9_cell_shape_ut.exo";
//                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 9);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::RECTANGULAR);
            }
        }

        TEST_CASE("MTK Cell Shape Quad9 Straight", "[MTK],[Cell_Shape],[Shape_Quad9_Straight]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                                          {1.5, 0.0},
                                          {1.0, 1.0},
                                          {0.0, 1.0},
                                          {0.5, 0.0},
                                          {1.25, 0.5},
                                          {0.5, 1.0},
                                          {0.0, 0.5},
                                          {0.5, 0.5}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD9;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeIds;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
//                std::string tFileOutput = "./mtk_quad9_cell_shape_ut.exo";
//                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 9);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::STRAIGHT);
            }
        }

        TEST_CASE("MTK Cell Shape Quad9 General", "[MTK],[Cell_Shape],[Shape_Quad9_General]")
        {
            if (par_size() <= 1)
            {
                //create a 2D MORIS mesh of quad4's using MTK database
                //------------------------------------------------------------------------------
                uint tNumDim = 2; // specify number of spatial dimensions

                // Node coordinate matrix
                Matrix<DDRMat> tCoords = {{0.0, 0.0},
                                          {1.5, 0.0},
                                          {1.0, 1.0},
                                          {0.0, 1.0},
                                          {0.5, 0.0},
                                          {1.0, 0.5},
                                          {0.5, 1.0},
                                          {0.0, 0.5},
                                          {0.5, 0.5}};

                Matrix<IndexMat> tNodeIndices = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};
                Matrix<IndexMat> tNodeIds = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

                // specify element connectivity of quad for mesh
                Matrix<IdMat> tElementConnQuad = {{1, 2, 3, 4, 5, 6, 7, 8, 9}};

                // specify the local to global element map for quads
                Matrix<IdMat> tElemLocalToGlobalQuad = {{1}};

                //------------------------------------------------------------------------------
                // create MORIS mesh using MTK database
                MtkMeshData aMeshData;
                aMeshData.CreateAllEdgesAndFaces = true;
                aMeshData.SpatialDim = &tNumDim;
                aMeshData.ElemConn(0) = &tElementConnQuad;
                aMeshData.CellTopology(0) = CellTopology::QUAD9;
                aMeshData.NodeCoords = &tCoords;
                aMeshData.LocaltoGlobalElemMap(0) = &tElemLocalToGlobalQuad;
                aMeshData.LocaltoGlobalNodeMap = &tNodeIds;

                Mesh *tMesh = create_interpolation_mesh(MeshType::STK, aMeshData);
                // Dump to file
//                std::string tFileOutput = "./mtk_quad9_cell_shape_ut.exo";
//                tMesh->create_output_mesh(tFileOutput);

                CHECK(tMesh->get_num_elems() == 1);

                // get the cell
                Cell & tCell = tMesh->get_mtk_cell(0);

                // Some checks on the cell
                CHECK(tCell.get_id() == 1);
                CHECK(tCell.get_index() == 0);
                CHECK(tCell.get_owner() == (moris_id)par_rank());
                CHECK(tCell.get_number_of_vertices() == 9);

                CHECK(tCell.get_cell_info()->compute_cell_shape(&tCell) == CellShape::GENERAL);
            }
        }

    } // namespace mtk
} // namespace moris
