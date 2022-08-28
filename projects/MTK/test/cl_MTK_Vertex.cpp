/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex.cpp
 *
 */

#include "catch.hpp"

// implementations to test
#include "cl_MTK_Vertex_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris
{
namespace mtk
{
TEST_CASE("MTK Vertex","[MTK],[MTK_VERTEX],[STK_VERTEX]")
{
    if(par_size() <=2)
    {
        // construct a mesh
        std::string tFilename = "generated:2x2x2";
        Mesh_Core_STK tMesh1( tFilename, NULL );

        uint tNodeInd = 3;
        uint tNodeId  = tMesh1.get_glb_entity_id_from_entity_loc_index(tNodeInd, EntityRank::NODE);

        // Setup vertex outside of mesh
        mtk::Vertex_STK tVertex(tNodeId,tNodeInd,&tMesh1);

        // Check index and id of vertex
        REQUIRE(tNodeId == (uint)tVertex.get_id());
        REQUIRE(tNodeInd == (uint)tVertex.get_index());

        // Check coordinates
        if(par_rank() == 0)
        {
            Matrix< DDRMat > tGoldCoords({{0.0, 1.0, 0.0}});
            Matrix< DDRMat > tVertexCoords = tVertex.get_coords();
            Matrix <DDBMat > tEqualCoords = (tGoldCoords == tVertexCoords);

            REQUIRE(tEqualCoords(0,0) == 1);
            REQUIRE(tEqualCoords(0,1) == 1);
            REQUIRE(tEqualCoords(0,2) == 1);
        }

        if(par_rank() == 1)
        {
            Matrix< DDRMat > tGoldCoords({{0.0, 1.0, 2.0}});
            Matrix< DDRMat > tVertexCoords = tVertex.get_coords();
            Matrix <DDBMat > tEqualCoords = (tGoldCoords == tVertexCoords);

            REQUIRE(tEqualCoords(0,0) == 1);
            REQUIRE(tEqualCoords(0,1) == 1);
            REQUIRE(tEqualCoords(0,2) == 1);
        }
    }
}

}
}

