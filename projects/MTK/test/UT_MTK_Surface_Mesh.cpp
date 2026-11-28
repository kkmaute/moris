/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 */

#include "catch.hpp"
#include "cl_MTK_Cell_ISC.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Vertex_Proxy.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Side_Cluster_DataBase.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_Json_Object.hpp"
#include <ostream>

// Exercise incremental face registration before geometry and normal evaluation.
#define private public
#include "cl_MTK_Surface_Mesh.hpp"
#undef private

namespace moris::mtk
{
    class Surface_Test_Mesh : public Integration_Mesh_DataBase_IG
    {
      public:
        uint get_spatial_dim() const override
        {
            return 2;
        }
    };

    TEST_CASE( "MTK surface face registration", "[MTK],[surface_face_registration]" )
    {
        // A single QUAD4 contributes one edge to the contact surface. Repeated
        // registration must preserve its connectivity; a different edge or
        // cluster must be rejected rather than silently losing a contribution.
        Surface_Test_Mesh          tMesh;
        Vector< Side_Set const * > tSideSets;
        Surface_Mesh               tSurface( &tMesh, tSideSets );
        auto                       tInfo = std::make_shared< Cell_Info_Quad4 >();
        Matrix< DDRMat >           tCoordinates;
        tInfo->get_loc_coords_of_cell( tCoordinates );

        Vector< Vertex_Proxy > tVertices( 4 );
        Vector< Vertex* >      tPointers( 4 );
        for ( uint i = 0; i < 4; ++i )
        {
            tVertices( i ).mVertexInd   = i;
            tVertices( i ).mVertexCoord = tCoordinates.get_row( i );
            tPointers( i )              = &tVertices( i );
        }

        Cell_ISC                                  tCell( 42, 17, 0, tInfo, tPointers );
        map< moris_index, Vector< moris_index > > tNeighbors;
        tSurface.mClusterToCellIndices.resize( 2 );
        tSurface.initialize_cell( tNeighbors, &tCell, 0, 0 );

        auto const tCheckOriginalFace = [ & ]() {
            REQUIRE( tSurface.get_number_of_cells() == 1 );
            REQUIRE( tSurface.get_number_of_vertices() == 2 );
            REQUIRE( tSurface.get_global_cell_index( 0 ) == 17 );
            REQUIRE( tSurface.get_local_cell_index( 17 ) == 0 );
            REQUIRE( tSurface.get_cluster_of_cell( 0 ) == 0 );
            REQUIRE( tSurface.mCellSideOrdinals( 0 ) == 0 );
            REQUIRE( tSurface.mClusterToCellIndices( 0 ).size() == 1 );
            REQUIRE( tSurface.mClusterToCellIndices( 1 ).empty() );
            auto const tFaceVertices    = tCell.get_geometric_vertices_on_side_ordinal( 0 );
            auto const tSurfaceVertices = tSurface.get_vertices_of_cell( 0 );
            REQUIRE( tSurfaceVertices.size() == 2 );
            for ( uint i = 0; i < 2; ++i )
            {
                REQUIRE( tSurface.get_global_vertex_index( tSurfaceVertices( i ) ) == tFaceVertices( i )->get_index() );
                REQUIRE( tNeighbors[ tFaceVertices( i )->get_index() ].size() == 1 );
            }
        };

        tCheckOriginalFace();
        tSurface.initialize_cell( tNeighbors, &tCell, 0, 0 );
        tCheckOriginalFace();

        REQUIRE_THROWS( tSurface.initialize_cell( tNeighbors, &tCell, 1, 0 ) );
        tCheckOriginalFace();

        REQUIRE_THROWS( tSurface.initialize_cell( tNeighbors, &tCell, 0, 1 ) );
        tCheckOriginalFace();
    }
}    // namespace moris::mtk
