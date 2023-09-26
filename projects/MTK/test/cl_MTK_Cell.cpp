/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell.cpp
 *
 */

#include "catch.hpp"
#include "cl_Communication_Tools.hpp"

// base class
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Cell_STK.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Vertex_STK.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
#include "fn_reindex_mat.hpp"

namespace moris
{
namespace mtk
{
TEST_CASE("MTK Cell","[MTK],[MTK_CELL]")
{
    if(par_size()<=1)
    {
        // construct a mesh
        std::string tFilename = "generated:2x2x2";
        Mesh_Core_STK tMesh1( tFilename, NULL );

        // get vertex information attached to element with index 0
        Matrix< IndexMat > tNodeIndices = tMesh1.get_entity_connected_to_entity_loc_inds(0, EntityRank::ELEMENT,EntityRank::NODE);
        Matrix< IdMat >   tNodeIds      = tMesh1.get_nodes_connected_to_element_glob_ids(1);

        // Setup Node Vertices (note: this data structure will be in the STK_Implementation
        moris::Cell<Vertex*> tElementVertices;
        for(size_t i =0; i<tNodeIndices.numel(); i++)
        {
            tElementVertices.push_back( new Vertex_STK(tNodeIds(i),tNodeIndices(i),&tMesh1) );
        }

        // Setup cell associated with element index 0
        std::shared_ptr<Cell_Info> tConn = std::make_shared<Cell_Info_Hex8>();
        Cell_STK tCell(tConn,
                       1,
                       0,
                       tElementVertices,
                       &tMesh1);

        REQUIRE(tCell.get_id() == 1);
        REQUIRE(tCell.get_index() == 0);

        REQUIRE(tCell.get_owner() == (moris_id)par_rank());
        REQUIRE(tCell.get_number_of_vertices() == 8);

        // verify ids and indices
        Matrix< IdMat >  tIdMat = tCell.get_vertex_ids();
        REQUIRE(all_true(tIdMat == tNodeIds));

        Matrix< IndexMat > tIndMat = tCell.get_vertex_inds();
        REQUIRE(all_true(tIndMat == tNodeIndices));

        if(par_rank() == 0)
        {
            Matrix< DDRMat > tGoldVertCoords
            ({{0.0, 0.0, 0.0},
                {1.0, 0.0, 0.0},
                {1.0, 1.0, 0.0},
                {0.0, 1.0, 0.0},
                {0.0, 0.0, 1.0},
                {1.0, 0.0, 1.0},
                {1.0, 1.0, 1.0},
                {0.0, 1.0, 1.0}});

            // Verify coordinates are as expected
            REQUIRE(all_true(tGoldVertCoords == tCell.get_vertex_coords()));
        }

        if(par_rank() == 1)
        {
            Matrix< DDRMat > tGoldVertCoords
           ({{0.0, 0.0, 1.0},
             {1.0, 0.0, 1.0},
             {1.0, 1.0, 1.0},
             {0.0, 1.0, 1.0},
             {0.0, 0.0, 2.0},
             {1.0, 0.0, 2.0},
             {1.0, 1.0, 2.0},
             {0.0, 1.0, 2.0}});

            // Verify coordinates are as expected
            REQUIRE(all_true(tGoldVertCoords == tCell.get_vertex_coords()));
        }

        // Gold Normals
        moris::Cell<moris::Matrix<moris::DDRMat>> tGoldNormals(6);
        tGoldNormals(0) = {{+0.000000000000000e+00},
                           {-1.000000000000000e+00},
                           {+0.000000000000000e+00}};

        tGoldNormals(1) = {{+1.000000000000000e+00},
                           {+0.000000000000000e+00},
                           {+0.000000000000000e+00}};

        tGoldNormals(2) = {{+0.000000000000000e+00},
                           {+1.000000000000000e+00},
                           {+0.000000000000000e+00}};

        tGoldNormals(3) = {{-1.000000000000000e+00},
                           {+0.000000000000000e+00},
                           {+0.000000000000000e+00}};

        tGoldNormals(4) = {{+0.000000000000000e+00},
                           {+0.000000000000000e+00},
                           {-1.000000000000000e+00}};

        tGoldNormals(5) = {{+0.000000000000000e+00},
                           {+0.000000000000000e+00},
                           {+1.000000000000000e+00}};

        // iterate through sides and verify normal
        for(moris::uint i = 0; i < 6 ; i ++)
        {
            Cell_Info_Hex8 tHex8;
            moris::Matrix<moris::IndexMat> tNodeToFaceMap = tHex8.get_node_to_face_map(i);

            moris::Matrix<moris::IndexMat> tNodeIdsOnFace = reindex_mat(tNodeToFaceMap,0,tIdMat);

            moris::Matrix<moris::DDRMat> tOutwardNormal = tCell.compute_outward_side_normal(i);

            CHECK(all_true(tOutwardNormal == tGoldNormals(i) ));
        }

        // ===================================================
        // Dump to file
        // ===================================================
        std::string tFileOutput = "./mtk_hex_cell_ut.exo";
        tMesh1.create_output_mesh(tFileOutput);

        // Delete because vertices were created with new call
        for (auto iT : tElementVertices)
        {
          delete iT;
        }
        tElementVertices.clear();

    }
}

TEST_CASE("MTK Cell Tet","[MTK],[MTK_CELL_TET]")
{
    if(par_size()<=1)
    {
        // construct a mesh
        moris::Matrix< moris::DDRMat >tNodeCoords( {{ 0.0,  0.0, 0.0}, { 1.0,  0.0, 0.0}, { 0.0,  1.0, 0.0},   { 0.0,  0.0, 1.0}});

        moris::Matrix<moris::IndexMat> tNodeIndices = {{0,1,2,3}};
        moris::Matrix<moris::IndexMat> tNodeIds     = {{1,2,3,4}};

        // 0D to 3D connectivity (node to element)
        Matrix< IdMat >     aElemConn = tNodeIds.copy();
        Matrix< IdMat >  aElemLocaltoGlobal = {{1}};

        // spatial dimension
        moris::uint tSpatialDim = 3;

        // Create MORIS mesh using MTK database
        MtkMeshData aMeshData;
        aMeshData.CreateAllEdgesAndFaces  = true;
        aMeshData.SpatialDim              = &tSpatialDim;
        aMeshData.ElemConn(0)             = &aElemConn;
        aMeshData.CellTopology(0)         = CellTopology::TET4;
        aMeshData.NodeCoords              = &tNodeCoords;
        aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

        Mesh* tMesh1  = create_interpolation_mesh( MeshType::STK, aMeshData );

        // Setup Node Vertices (note: this data structure will be in the STK_Implementation
        moris::Cell<Vertex*> tElementVertices;
        for(size_t i =0; i<tNodeIndices.numel(); i++)
        {
            tElementVertices.push_back( new Vertex_STK(tNodeIds(i),tNodeIndices(i),tMesh1) );
        }

        // Setup cell associated with element index 0
        std::shared_ptr<Cell_Info> tTet4Con = std::make_shared<Cell_Info_Tet4>();
        Cell_STK tCell(tTet4Con, 1, 0, tElementVertices, tMesh1);

        REQUIRE(tCell.get_id() == 1);
        REQUIRE(tCell.get_index() == 0);

        REQUIRE(tCell.get_owner() == (moris_id)par_rank());
        REQUIRE(tCell.get_number_of_vertices() == 4);

        // verify ids and indices
        Matrix< IdMat >  tIdMat = tCell.get_vertex_ids();
        REQUIRE(all_true(tIdMat == tNodeIds));

        Matrix< IndexMat > tIndMat = tCell.get_vertex_inds();
        REQUIRE(all_true(tIndMat == tNodeIndices));

        // coordinate checking
        Matrix< DDRMat > tGoldVertCoords = {{0.0, 0.0, 0.0},
                                            {1.0, 0.0, 0.0},
                                            {0.0, 1.0, 0.0},
                                            {0.0, 0.0, 1.0}};

        // Verify coordinates are as expected
        REQUIRE(all_true(tGoldVertCoords == tCell.get_vertex_coords()));

        // Gold Normals
        moris::uint tNumSides = 4;
        moris::Cell<moris::Matrix<moris::DDRMat>> tGoldNormals(tNumSides);
        tGoldNormals(0) = {{+0.000000000000000e+00},
                           {-1.000000000000000e+00},
                           {+0.000000000000000e+00}};

        tGoldNormals(1) = {{+5.773502691896258e-01},
                           {+5.773502691896258e-01},
                           {+5.773502691896258e-01}};

        tGoldNormals(2) = {{-1.000000000000000e+00},
                           {+0.000000000000000e+00},
                           {+0.000000000000000e+00}};

        tGoldNormals(3) = {{+0.000000000000000e+00},
                           {+0.000000000000000e+00},
                           {-1.000000000000000e+00}};

        // Check normals
        // iterate through sides and verify normal
        for(moris::uint i = 0; i < tNumSides ; i ++)
        {
            moris::Matrix<moris::DDRMat> tOutwardNormal = tCell.compute_outward_side_normal(i);

            CHECK(all_true(tOutwardNormal == tGoldNormals(i) ));
        }

        // Dump to file
        std::string tFileOutput = "./mtk_tet_cell_ut.exo";
        tMesh1->create_output_mesh(tFileOutput);

        // Delete because vertices were created with new call
        for (auto iT : tElementVertices)
        {
          delete iT;
        }
        tElementVertices.clear();

        delete tMesh1;

    }
}

}
}

