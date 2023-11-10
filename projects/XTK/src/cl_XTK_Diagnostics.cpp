/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Diagnostics.cpp
 *
 */

#include "cl_XTK_Diagnostics.hpp"

#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_norm.hpp"

using namespace moris;
namespace xtk
{
bool
interpolated_coordinate_check( Cut_Integration_Mesh* aCutMesh )
{
    moris::Cell< std::shared_ptr< Matrix< DDRMat > > >* tCoords = aCutMesh->get_all_vertex_coordinates_loc_inds();

    // keep track of problem vertices
    moris::Cell< moris_index > tProblemIgVertices( 0 );

    // Allocate a basis function weight matrix
    moris::Matrix< moris::DDRMat > tBasisWeights( 1, 8 );

    // tolerance for difference between coordinates
    real tTol = 1e-12;

    // iterate through
    moris::uint tNumIgCellGroups = aCutMesh->get_num_ig_cell_groups();

    // iterate over ig cell groups
    for ( moris::uint iIgCellGroup = 0; iIgCellGroup < tNumIgCellGroups; iIgCellGroup++ )
    {
        // get the cell group
        std::shared_ptr< IG_Cell_Group > tIgCellGroup = aCutMesh->get_ig_cell_group( (moris_index)iIgCellGroup );

        // get the vertex group
        std::shared_ptr< IG_Vertex_Group > tIgVertexGroup = aCutMesh->get_vertex_group( (moris_index)iIgCellGroup );

        // the base cell
        moris::mtk::Cell* tBaseCell = aCutMesh->get_ig_cell_group_parent_cell( (moris_index)iIgCellGroup );

        moris::Matrix< moris::DDRMat > tBaseCellCoords = tBaseCell->get_vertex_coords();

        moris::mtk::Cell_Info const* tCellInfo = tBaseCell->get_cell_info();

        // iterate through the IgCellGroup
        for ( moris::uint iCell = 0; iCell < tIgCellGroup->mIgCellGroup.size(); iCell++ )
        {
            moris::mtk::Cell* tIgCell = tIgCellGroup->mIgCellGroup( iCell );

            moris::Cell< moris::mtk::Vertex* > tVertices = tIgCell->get_vertex_pointers();

            // iterate through vertices attached to the cell
            for ( moris::uint iV = 0; iV < tVertices.size(); iV++ )
            {
                // local coordinate of this vertex wrt the current cell group
                std::shared_ptr< Matrix< DDRMat > > tVertexLocalCoords = tIgVertexGroup->get_vertex_local_coords( tVertices( iV )->get_index() );

                // evaluate the basis function
                tCellInfo->eval_N( *tVertexLocalCoords, tBasisWeights );

                // Evaluate the nodes global coordinate from the basis weights
                moris::Matrix< moris::DDRMat > tInterpNodeCoord = tBasisWeights * tBaseCellCoords;

                // node index
                moris_index tNodeIndex = tVertices( iV )->get_index();

                if ( norm( tInterpNodeCoord - *( ( *tCoords )( tNodeIndex ) ) ) > tTol )
                {
                    std::cout << "Problem Vertex: " << tVertices( iV )->get_id();
                    std::cout << " | Group Index: " << iIgCellGroup;
                    std::cout << " | Interp Coords: ";

                    for ( moris::uint iSp = 0; iSp < tInterpNodeCoord.numel(); iSp++ )
                    {
                        std::cout << std::left << std::scientific << tInterpNodeCoord( iSp ) << " ";
                    }

                    std::cout << " | Coords: ";

                    moris::Matrix< moris::DDRMat > tCoords = tVertices( iV )->get_coords();

                    for ( moris::uint iSp = 0; iSp < tCoords.numel(); iSp++ )
                    {
                        std::cout << std::left << std::scientific << tCoords( iSp ) << " ";
                    }

                    std::cout << " | xsi: ";

                    for ( moris::uint iSp = 0; iSp < tVertexLocalCoords->numel(); iSp++ )
                    {
                        std::cout << std::left << std::scientific << ( *tVertexLocalCoords )( iSp ) << " ";
                    }
                    std::cout << std::endl;
                    tProblemIgVertices.push_back( iV );
                }
            }
        }
    }

    // a problem vertex can appear multiple times if is an issue in more than one ig cell group
    return tProblemIgVertices.size() == 0;
}

bool
verify_interface_vertices(
    xtk::Model*                           aModel,
    moris::Matrix< moris::DDRMat > const& aIsocontourThreshold,
    moris::Matrix< moris::DDRMat > const& aIsocontourTolerance )
{
    // get the cut integration mesh
    Cut_Integration_Mesh* tCutMesh = aModel->get_cut_integration_mesh();

    // get the geometry engine
    moris::ge::Geometry_Engine* tGeomEngine = aModel->get_geom_engine();

    // number of vertices
    moris::uint tNumVertices = tCutMesh->get_num_entities( mtk::EntityRank::NODE, 0 );

    // keep track of problem vertices
    moris::Cell< moris_index > tProblemIgVertices( 0 );

    // iterate through vertices in the integration mesh
    for ( moris::uint iV = 0; iV < tNumVertices; iV++ )
    {
        // iterate through geometries
        for ( moris::uint iGeom = 0; iGeom < tGeomEngine->get_num_geometries(); iGeom++ )
        {
            moris::mtk::Vertex* tVertex = &tCutMesh->get_mtk_vertex( (moris_index)iV );
            if ( tGeomEngine->is_interface_vertex( tVertex->get_index(), iGeom ) )
            {
                moris::real tGeomVal    = tGeomEngine->get_field_value( iGeom, tVertex->get_index(), tVertex->get_coords() );
                moris::real tDifference = std::abs( tGeomVal - aIsocontourThreshold( iGeom ) );
                if ( tDifference > aIsocontourTolerance( iGeom ) )
                {
                    tProblemIgVertices.push_back( tVertex->get_index() );

                    std::cout << "Interface Issue:"
                              << " Vert Index: " << tVertex->get_index() << " | Vert Id: " << tVertex->get_id() << " | Owner: " << tVertex->get_owner() << " | diff: " << tDifference << " | par_rank: " << par_rank() << std::endl;
                }
            }
        }
    }

    return tProblemIgVertices.size() == 0;
}

bool
check_vertices(
    xtk::Model*                                           aModel,
    moris::Cell< moris_index > const&                     aGoldNumVerts,
    std::unordered_map< moris_index, moris_index > const& aGoldVertexMap,
    moris::Matrix< moris::DDRMat > const&                 aGoldVertexCoords,
    moris::real                                           aTolerance )
{
    xtk::Cut_Integration_Mesh* tCutIgMesh = aModel->get_cut_integration_mesh();

    if ( (moris_index)tCutIgMesh->get_num_nodes() != aGoldNumVerts( moris::par_rank() ) )
    {
        std::cout << "Number of vertex mismatch on proc" << par_rank() << std::endl;
        return false;
    }
    // iterate through verts
    for ( moris::uint iV = 0; iV < tCutIgMesh->get_num_nodes(); iV++ )
    {
        moris::mtk::Vertex const& tVertex = tCutIgMesh->get_mtk_vertex( (moris_index)iV );
        auto                      tIter   = aGoldVertexMap.find( tVertex.get_id() );

        MORIS_ERROR( tIter != aGoldVertexMap.end(), "Vertex not in map" );

        moris::Matrix< moris::DDRMat > tVertexCoords = tVertex.get_coords();

        if ( moris::norm( tVertexCoords - aGoldVertexCoords.get_row( tIter->second ) ) > aTolerance )
        {
            std::cout << "tVertex.get_id() = " << tVertex.get_id() << " | gold index = " << tIter->second << std::endl;
            moris::print( aGoldVertexCoords.get_row( tIter->second ), "Gold Coordinates" );
            moris::print( tVertexCoords, "Coordinates" );
            std::cout << "Vertex coordinate failure on " << par_rank() << std::endl;
            return false;
        }
    }

    return true;
}

bool
check_cells(
    xtk::Model*                                                         aModel,
    moris::Cell< moris::moris_index > const&                            aGoldNumCells,
    std::unordered_map< moris::moris_index, moris::moris_index > const& aGoldCellMap,
    moris::Matrix< moris::IndexMat > const&                             aGoldCellConn )
{
    xtk::Cut_Integration_Mesh* tCutIgMesh = aModel->get_cut_integration_mesh();

    if ( (moris_index)tCutIgMesh->get_num_elems() != aGoldNumCells( moris::par_rank() ) )
    {
        std::cout << "Number of cell mismatch on proc" << par_rank() << std::endl;
        return false;
    }
    // iterate through verts
    for ( moris::uint iCell = 0; iCell < tCutIgMesh->get_num_elems(); iCell++ )
    {
        moris::mtk::Cell const& tCell = tCutIgMesh->get_mtk_cell( (moris_index)iCell );

        auto tIter = aGoldCellMap.find( tCell.get_id() );

        moris_index const tGoldCellIndex = tIter->second;

        MORIS_ERROR( tIter != aGoldCellMap.end(), "Cell not in map" );

        moris::Matrix< moris::IndexMat > tCellVertIds = tCell.get_vertex_ids();
        for ( moris::uint iV = 0; iV < tCellVertIds.numel(); iV++ )
        {
            if ( tCellVertIds( iV ) != aGoldCellConn( tGoldCellIndex, iV ) )
            {
                std::cout << "Cell connectivity failure on " << par_rank() << std::endl;
                return false;
            }
        }
    }
    return true;
}
}// namespace xtk
