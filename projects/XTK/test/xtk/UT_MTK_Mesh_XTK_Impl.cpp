/*
 * cl_MTK_Mesh_XTK_Impl_UT.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"

#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

namespace xtk
{

//TEST_CASE("Unit test of XTK's Interface with MTK","[XTK_MTK]")
//        {
//    if(par_size() < 2)
//    {
//
//        real tRadius  = 0.55;
//        real tXCenter = 2.0;
//        real tYCenter = 2.0;
//        real tZCenter = 2.0;
//        Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
//        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//        Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);
//
//        // Create Mesh --------------------------------------------------------------------
//        std::string tMeshFileName = "generated:2x2x2";
//        moris::mtk::Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName );
//
//        // Setup XTK Model ----------------------------------------------------------------
//        size_t tModelDimension = 3;
//        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
//        tXTKModel.mVerbose  =  false;
//
//        //Specify decomposition Method and Cut Mesh ---------------------------------------
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//        tXTKModel.decompose(tDecompositionMethods);
//        tXTKModel.unzip_interface();
//        tXTKModel.perform_basis_enrichment();
//
//        // return XTK as an mtk mesh
//        moris::mtk::Mesh* tXTKToMTK = tXTKModel.get_xtk_as_mtk();
//
//        // return XTK as an mtk mesh implemented by STK
//        moris::mtk::Mesh* tStkMesh = tXTKModel.get_output_mesh();
//
//
//        // Verify that the STK and XTK meshes return the same thing
//        // for certain
//        moris::uint tNumElemsXTK = tXTKToMTK->get_num_entities(EntityRank::ELEMENT);
//        moris::uint tNumElemsSTK = tStkMesh->get_num_entities(EntityRank::ELEMENT);
//
//        // TODO: figure out how to handle interface elements
//        //    CHECK(tNumElemsXTK == tNumElemsSTK);
//
//        moris::uint tNumNodesXTK = tXTKToMTK->get_num_entities(EntityRank::NODE);
//        moris::uint tNumNodesSTK = tStkMesh->get_num_entities(EntityRank::NODE);
//        CHECK(tNumNodesSTK == tNumNodesXTK);
//
//        // get cells by index and verify their indices are the same
//        // check element to node connectivity
//        for(moris::uint iEl = 0; iEl<tNumElemsXTK; iEl++)
//        {
//            mtk::Cell  & tCell = tXTKToMTK->get_mtk_cell(iEl);
//            CHECK(tCell.get_index() == (moris_index)iEl);
//
//            // Get stk index using xtk cell id
//            moris::moris_index tXTKCellId = tCell.get_id();
//            moris::moris_index tSTKCellIndex = tStkMesh->get_loc_entity_ind_from_entity_glb_id(tXTKCellId,EntityRank::ELEMENT);
//
//            Matrix<IndexMat> tXTKCellNodesInds   = tCell.get_vertex_inds();
//            Matrix<IndexMat> tSTKCellNodesInds   = tStkMesh->get_entity_connected_to_entity_loc_inds(tSTKCellIndex,EntityRank::ELEMENT,EntityRank::NODE);
//            Matrix<IndexMat> tXTKCellToNodeInds2 = tXTKToMTK->get_entity_connected_to_entity_loc_inds(iEl,EntityRank::ELEMENT,EntityRank::NODE);
//            CHECK(all_true(tXTKCellNodesInds == tSTKCellNodesInds));
//            CHECK(all_true(tXTKCellToNodeInds2 == tXTKCellNodesInds));
//
//            Matrix<IndexMat> tXTKCellNodesIds = tCell.get_vertex_ids();
//            Matrix<IndexMat> tSTKCellNodesIds = tStkMesh->get_entity_connected_to_entity_glob_ids(tXTKCellId,EntityRank::ELEMENT,EntityRank::NODE);
//            CHECK(all_true(tXTKCellNodesIds == tSTKCellNodesIds));
//        }
//
//        // get vertex by index and verify their indices are the same
//        // verify their coordinates too
//        for(moris::uint iN = 0; iN<tNumNodesXTK; iN++)
//        {
//            mtk::Vertex  & tVertex = tXTKToMTK->get_mtk_vertex(iN);
//            CHECK(tVertex.get_index() == (moris_index)iN);
//            CHECK(tVertex.get_id() == tStkMesh->get_glb_entity_id_from_entity_loc_index(iN,EntityRank::NODE));
//
//            moris::Matrix<moris::DDRMat> tXTKCoord = tVertex.get_coords();
//            moris::Matrix<moris::DDRMat> tSTKCoord = tStkMesh->get_node_coordinate((moris_index)iN);
//
//            CHECK(all_true(tXTKCoord == tSTKCoord));
//
//            // Accessing t-matrix
//            moris::mtk::Vertex_Interpolation * tVertexInterpolation = tVertex.get_interpolation(0);
//
//            moris::print(*tVertexInterpolation->get_weights()," weights ");
//            moris::print(tVertexInterpolation->get_indices()," indices ");
//
//        }
//
//
//
//        delete tXTKToMTK;
//        delete tMeshData;
//        delete tStkMesh;
//    }
//        }
}

