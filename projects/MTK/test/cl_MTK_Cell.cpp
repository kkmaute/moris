/*
 * cl_MTK_Cell.cpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */


#include "catch.hpp"

// base class
#include"cl_MTK_Cell.hpp"

// implementations to test
#include "cl_MTK_Cell_STK.hpp"
#include "cl_MTK_Mesh_STK_New.hpp"
#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Vertex_STK.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "fn_print.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"
namespace moris
{
namespace mtk
{
TEST_CASE("MTK Cell","[MTK],[MTK_CELL]")
{
    // construct a mesh
    std::string tFilename = "generated:2x2x2";
    Mesh_STK tMesh1( tFilename, NULL );

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
    Cell_STK tCell(CellType::HEX8,
                   1,
                   0,
                   tElementVertices,
                   &tMesh1);

    REQUIRE(tCell.get_id() == 1);
    REQUIRE(tCell.get_index() == 0);

    REQUIRE(tCell.get_owner() == 0);
    REQUIRE(tCell.get_number_of_vertices() == 8);

    // verify ids and indices
    Matrix< IdMat >  tIdMat = tCell.get_vertex_ids();
    REQUIRE(all_true(tIdMat == tNodeIds));

    Matrix< IndexMat > tIndMat = tCell.get_vertex_inds();
    REQUIRE(all_true(tIndMat == tNodeIndices));

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

    // Delete because vertices were created with new call
    for(size_t i =0; i<tNodeIndices.numel(); i++)
    {
        delete tElementVertices(i);
    }

}
}
}
