/*
 * UT_MTK_Mesh_Comm_Maps.cpp
 *
 *  Created on: Jun 3, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Visualization_STK.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

namespace moris
{

TEST_CASE("MTK Vertex Communication Tables","[VERTEX_COMM_TABLES]")
{
    // create mesh
    const std::string tMeshFile = "generated:1x1x4";

    // add parallel fields to the mesh
    moris::mtk::Visualization_STK tVizTool;
    moris::mtk::MtkFieldsInfo* tFieldsInfo = tVizTool.setup_parallel_cell_fields_for_declaration();

    // Declare some supplementary fields
    moris::mtk::MtkMeshData tMeshData;
    tMeshData.FieldsInfo = tFieldsInfo;

    // Create MORIS mesh using MTK database
    moris::mtk::Mesh* tMesh = mtk::create_mesh( MeshType::STK, tMeshFile, &tMeshData );

    // fill in the parallel fields
    tVizTool.populate_parallel_cell_fields_on_mesh(tMesh);

    // Output the mesh
    std::string tOutputMeshFile = "./mtk_comm_tables.exo";
    tMesh->create_output_mesh(tOutputMeshFile);

    if(par_size() == 2)
    {
        if(par_rank() == 0)
        {
            Matrix<IdMat> tVertexPairsWithProc1Gold = {{5,  12},
                                                       {6,  13},
                                                       {7,  14},
                                                       {8,  15},
                                                       {9,   8},
                                                       {10,  9},
                                                       {11, 10},
                                                       {12, 11},
                                                       {13,  0},
                                                       {14,  1},
                                                       {15,  2},
                                                       {16,  3}};

            Matrix<IdMat> tVertexPairsWithProc2Gold = {{5, 4},
                                                       {6, 5},
                                                       {7, 6},
                                                       {8, 7}};

            Matrix<IdMat> tProcsWithVertexPairsGold = {{1}};

            CHECK(all_true(tVertexPairsWithProc1Gold == tMesh->get_communication_vertex_pairing()(0)));
            CHECK(all_true(tProcsWithVertexPairsGold == tMesh->get_communication_proc_ranks()));

        }

        if(par_rank() == 1)
        {
            Matrix<IdMat> tVertexPairsWithProc0Gold = {{5, 4},
                                                       {6, 5},
                                                       {7, 6},
                                                       {8, 7},
                                                       {9, 12},
                                                       {10, 13},
                                                       {11, 14},
                                                       {12, 15},
                                                       {13, 8},
                                                       {14, 9},
                                                       {15, 10},
                                                       {16, 11}};

            Matrix<IdMat> tProcsWithVertexPairsGold = {{0}};

            CHECK(all_true(tVertexPairsWithProc0Gold == tMesh->get_communication_vertex_pairing()(0)));
            CHECK(all_true(tProcsWithVertexPairsGold == tMesh->get_communication_proc_ranks()));
        }
    }

    if(par_size() == 4)
    {
        if(par_rank() == 0)
        {
            Matrix<IdMat> tVertexPairsWithProc1Gold = {{1, 4},
                                                       {2, 5},
                                                       {3, 6},
                                                       {4, 7},
                                                       {5, 0},
                                                       {6, 1},
                                                       {7, 2},
                                                       {8, 3},
                                                       {9, 12},
                                                       {10, 13},
                                                       {11, 14},
                                                       {12, 15}};


            Matrix<IdMat> tVertexPairsWithProc2Gold = {{5, 4},
                                                       {6, 5},
                                                       {7, 6},
                                                       {8, 7}};

            Matrix<IdMat> tProcsWithVertexPairsGold = {{1},
                                                       {2}};

            CHECK(all_true(tVertexPairsWithProc1Gold == tMesh->get_communication_vertex_pairing()(0)));
            CHECK(all_true(tVertexPairsWithProc2Gold == tMesh->get_communication_vertex_pairing()(1)));
            CHECK(all_true(tProcsWithVertexPairsGold == tMesh->get_communication_proc_ranks()));

        }

        if(par_rank() == 1)
        {
            Matrix<IdMat> tVertexPairsWithProc0Gold = {{1, 0},
                                                       {2, 1},
                                                       {3, 2},
                                                       {4, 3},
                                                       {5, 8},
                                                       {6, 9},
                                                       {7, 10},
                                                       {8, 11},
                                                       {9, 4},
                                                       {10, 5},
                                                       {11, 6},
                                                       {12, 7}};

            Matrix<IdMat> tVertexPairsWithProc2Gold = {{9, 0},
                                                       {10, 1},
                                                       {11, 2},
                                                       {12, 3},
                                                       {13, 12},
                                                       {14, 13},
                                                       {15, 14},
                                                       {16, 15}};

            Matrix<IdMat> tVertexPairsWithProc3Gold = {{9, 8},
                                                       {10, 9},
                                                       {11, 10},
                                                       {12, 11}};

            Matrix<IdMat> tProcsWithVertexPairsGold = {{0},{2},{3}};

            CHECK(all_true(tVertexPairsWithProc0Gold == tMesh->get_communication_vertex_pairing()(0)));
            CHECK(all_true(tVertexPairsWithProc2Gold == tMesh->get_communication_vertex_pairing()(1)));
            CHECK(all_true(tVertexPairsWithProc3Gold == tMesh->get_communication_vertex_pairing()(2)));
            CHECK(all_true(tProcsWithVertexPairsGold == tMesh->get_communication_proc_ranks()));
        }

        if(par_rank() == 2)
        {
            Matrix<IdMat> tVertexPairsWithProc0Gold = {{5, 8},
                                                       {6, 9},
                                                       {7, 10},
                                                       {8, 11}};

            Matrix<IdMat> tVertexPairsWithProc1Gold = {{9, 12},
                                                       {10, 13},
                                                       {11, 14},
                                                       {12, 15},
                                                       {13, 8},
                                                       {14, 9},
                                                       {15, 10},
                                                       {16, 11}};

            Matrix<IdMat> tVertexPairsWithProc3Gold = {{13, 4},
                                                       {14, 5},
                                                       {15, 6},
                                                       {16, 7},
                                                       {17, 0},
                                                       {18, 1},
                                                       {19, 2},
                                                       {20, 3}};

            Matrix<IdMat> tProcsWithVertexPairsGold = {{0},{1},{3}};

            CHECK(all_true(tVertexPairsWithProc0Gold == tMesh->get_communication_vertex_pairing()(0)));
            CHECK(all_true(tVertexPairsWithProc1Gold == tMesh->get_communication_vertex_pairing()(1)));
            CHECK(all_true(tVertexPairsWithProc3Gold == tMesh->get_communication_vertex_pairing()(2)));
            CHECK(all_true(tProcsWithVertexPairsGold == tMesh->get_communication_proc_ranks()));
        }

        if(par_rank() == 3)
        {
            Matrix<IdMat> tVertexPairsWithProc1Gold = {{9, 12},
                                                       {10, 13},
                                                       {11, 14},
                                                       {12, 15}};


            Matrix<IdMat> tVertexPairsWithProc2Gold = {{13, 12},
                                                       {14, 13},
                                                       {15, 14},
                                                       {16, 15},
                                                       {17, 8},
                                                       {18, 9},
                                                       {19, 10},
                                                       {20, 11}};

            Matrix<IdMat> tProcsWithVertexPairsGold = {{1},
                                                       {2}};

            CHECK(all_true(tVertexPairsWithProc1Gold == tMesh->get_communication_vertex_pairing()(0)));
            CHECK(all_true(tVertexPairsWithProc2Gold == tMesh->get_communication_vertex_pairing()(1)));
            CHECK(all_true(tProcsWithVertexPairsGold == tMesh->get_communication_proc_ranks()));
        }
    }
}

TEST_CASE("MTK Cell Communication Tables","[CELL_COMM_TABLES]")
{
    // create mesh
    const std::string tMeshFile = "generated:1x1x4";

    // add parallel fields to the mesh
    moris::mtk::Visualization_STK tVizTool;
    moris::mtk::MtkFieldsInfo* tFieldsInfo = tVizTool.setup_parallel_cell_fields_for_declaration();

    // Declare some supplementary fields
    moris::mtk::MtkMeshData tMeshData;
    tMeshData.FieldsInfo = tFieldsInfo;

    // Create MORIS mesh using MTK database
    moris::mtk::Mesh* tMesh = mtk::create_mesh( MeshType::STK, tMeshFile, &tMeshData );

    // fill in the parallel fields
//    tVizTool.populate_parallel_cell_fields_on_mesh(tMesh);

    // Output the mesh
    std::string tOutputMeshFile = "./mtk_comm_tables.exo";
    tMesh->create_output_mesh(tOutputMeshFile);

    moris_index tCellId = 2;

    if(par_rank() == 0)
    {
        moris_index tCellIndex = tMesh->get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);
        Matrix<IndexMat> tSharingProcs(1,1);
        tMesh->get_processors_whom_share_entity(tCellIndex,EntityRank::ELEMENT,tSharingProcs);

        moris::print_std_initializer_list(tSharingProcs,"sharing on 0");
    }

    barrier("lll");

    if(par_rank() == 1)
    {
        moris_index tCellIndex = tMesh->get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);
        Matrix<IndexMat> tSharingProcs(1,1);
        tMesh->get_processors_whom_share_entity(tCellIndex,EntityRank::ELEMENT,tSharingProcs);
        moris::print_std_initializer_list(tSharingProcs,"sharing on 1");
    }

    barrier("lll");

    if(par_rank() == 2)
    {
        moris_index tCellIndex = tMesh->get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);
        Matrix<IndexMat> tSharingProcs(1,1);
        tMesh->get_processors_whom_share_entity(tCellIndex,EntityRank::ELEMENT,tSharingProcs);
        moris::print_std_initializer_list(tSharingProcs,"sharing on 2");
    }

    barrier("lll");
    if(par_rank() == 3)
    {
//        moris::print_std_initializer_list(tMesh->get_communication_proc_ranks(),"proc ranks on p3");
//        moris::print_std_initializer_list(tMesh->get_communication_vertex_pairing()(0),"vert pairs on p3 ");
//        moris::print_std_initializer_list(tMesh->get_communication_vertex_pairing()(1),"vert pairs on p3 ");
    }

    barrier("lll");

    delete tMesh;

}

}
