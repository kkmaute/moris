//
// Created by frank on 12/4/23.
//

#include <unordered_map>
#include "cl_MTK_Json_Debug_Output.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "iostream"
#include "cl_Json_Object.hpp"

namespace moris
{
    namespace mtk
    {
        void Json_Debug_Output::write_to_json( std::string const &aFileName )
        {
            Json tMeshObj;
            mAllIGVertices.clear();
            mAllIPVertices.clear();
            mAllIGCells.clear();
            mAllIPCells.clear();

            Json tSets;

            // get all block sets
            auto tBlockSets = mMesh->get_block_sets();
            tSets.push_back( { "block_sets", serialize_block_sets( tBlockSets ) } );

            // side sets
            auto tSideSets = mMesh->get_side_sets();
            tSets.push_back( { "side_sets", serialize_side_sets( tSideSets ) } );

            // double side sets
            auto tDoubleSideSets = mMesh->get_double_side_sets();
            tSets.push_back( { "double_side_sets", serialize_double_side_sets( tDoubleSideSets ) } );

            // store vertices
            tMeshObj.push_back( { "sets", tSets } );
            tMeshObj.push_back( { "vertices_ig", serialize_all_vertices( mAllIGVertices, mIGVertexDisplacements ) } );
            tMeshObj.push_back( { "vertices_ip", serialize_all_vertices( mAllIPVertices, mIPVertexDisplacements ) } );
            tMeshObj.push_back( { "cells_ig", serialize_all_cells( mAllIGCells ) } );
            tMeshObj.push_back( { "cells_ip", serialize_all_cells( mAllIPCells ) } );

            write_json( aFileName, tMeshObj );
        }

        Json Json_Debug_Output::serialize_all_vertices(
                const std::unordered_set< Vertex const * >                  &aVertices,
                const std::unordered_map< moris::moris_index, Vector< moris::real > > &aVertexDisplacements )
        {
            Json tVertexObj;

            Vector< moris_id >    tIds;
            Vector< moris_index > tIndices;
            Vector< real >        tX;
            Vector< real >        tY;
            Vector< real >        tZ;
            for ( auto const *tVertex : aVertices )
            {
                // Add values to corresponding arrays
                tIds.push_back( tVertex->get_id() );
                tIndices.push_back( tVertex->get_index() );

                Vector< real > tDisplacement;
                if ( aVertexDisplacements.find( tVertex->get_index() ) != aVertexDisplacements.end() )
                {
                    tDisplacement = aVertexDisplacements.at( tVertex->get_index() );
                    MORIS_ASSERT( tDisplacement.size() == mMesh->get_spatial_dim(), "Displacement vector size does not match mesh dimension" );
                }
                else
                {
                    tDisplacement = { 0.0, 0.0, 0.0 };
                }

                tX.push_back( tVertex->get_coords()( 0 ) + tDisplacement( 0 ) );
                tY.push_back( tVertex->get_coords()( 1 ) + tDisplacement( 1 ) );
                if ( mMesh->get_spatial_dim() == 3 )
                {
                    tZ.push_back( tVertex->get_coords()( 2 ) + tDisplacement( 2 ) );
                }
            }

            tVertexObj.put_child( "id", to_json( tIds ) );
            tVertexObj.put_child( "index", to_json( tIndices ) );
            tVertexObj.put_child( "x", to_json( tX ) );
            tVertexObj.put_child( "y", to_json( tY ) );

            if ( mMesh->get_spatial_dim() == 3 )
            {
                tVertexObj.put_child( "z", to_json( tZ ) );
            }
            return tVertexObj;
        }

        Json Json_Debug_Output::serialize_all_cells( const std::unordered_set< Cell const * > &aCells )
        {
            Json tAllCells;

            for ( auto *tCell : aCells )
            {
                Json tCellObj;
                // Assuming get_id and get_index are member functions of Cell
                tCellObj.put( "id", tCell->get_id() );
                tCellObj.put( "index", tCell->get_index() );

                Json tFacets;
                for ( unsigned iSide = 0; iSide < tCell->get_number_of_facets(); ++iSide )
                {
                    Json tSideObj;
                    tSideObj.put( "ordinal", iSide );

                    Json tVertexIds;
                    auto tSideVertices = tCell->get_vertices_on_side_ordinal( iSide );
                    for ( auto *tVertex : tSideVertices )
                    {
                        Json tVertexId;
                        tVertexId.put( "", tVertex->get_id() );
                        tVertexIds.push_back( { "", tVertexId } );
                    }
                    tSideObj.push_back( { "vertex_ids", tVertexIds } );
                    tFacets.push_back( { "", tSideObj } );
                }
                tCellObj.push_back( { "facets", tFacets } );
                tAllCells.push_back( { "", tCellObj } );
            }
            return tAllCells;
        }

        Json Json_Debug_Output::serialize_double_side_sets( Vector< moris::mtk::Double_Side_Set * > &aDoubleSideSets )
        {
            Json tDoubleSideSetsObj;
            for ( auto tSideSet : aDoubleSideSets )
            {
                Json tSideSetObj;
                tSideSetObj.put( "name", tSideSet->get_set_name() );
                tSideSetObj.put( "index", tSideSet->get_set_index() );

                Json tClusterObj;
                for ( auto tCluster : tSideSet->get_clusters_on_set() )
                {
                    auto tDblCluster = dynamic_cast< Double_Side_Cluster const * >( tCluster );
                    Json tDblClusterObj;
                    tDblClusterObj.put_child( "follower_cluster", serialize_side_cluster( &tDblCluster->get_follower_side_cluster() ) );
                    tDblClusterObj.put_child( "leader_cluster", serialize_side_cluster( &tDblCluster->get_leader_side_cluster() ) );
                    tClusterObj.push_back( { "", tDblClusterObj } );
                }
                tSideSetObj.push_back( { "clusters", tClusterObj } );
                tDoubleSideSetsObj.push_back( { tSideSet->get_set_name(), tSideSetObj } );
            }
            return tDoubleSideSetsObj;
        }

        Json Json_Debug_Output::serialize_block_sets( Vector< moris::mtk::Block_Set * > &tBlockSets )
        {
            Json tBlockSetsObj;
            for ( auto tBlockSet : tBlockSets )
            {
                Json tBlockSetObj;
                tBlockSetObj.put( "name", tBlockSet->get_set_name() );
                tBlockSetObj.put( "index", tBlockSet->get_set_index() );

                auto tClusters = tBlockSet->get_clusters_on_set();
                tBlockSetObj.push_back( { "clusters", serialize_cell_clusters( tClusters ) } );

                tBlockSetsObj.push_back( { tBlockSet->get_set_name(), tBlockSetObj } );
            }
            return tBlockSetsObj;
        }

        Json Json_Debug_Output::serialize_cell_clusters( Vector< Cluster const * > &tClusters )
        {
            Json tClusterArr;
            for ( auto tCluster : tClusters )
            {
                tClusterArr.push_back( { "", serialize_cell_cluster( tCluster ) } );
            }
            return tClusterArr;
        }

        Json Json_Debug_Output::serialize_cell_cluster( Cluster const *aCluster )
        {
            Json tClusterObj;

            auto tDataBase = dynamic_cast< Cell_Cluster const * >( aCluster );

            auto tPrimaryCells = tDataBase->get_primary_cells_in_cluster();
            tClusterObj.put_child( "primary_ig_cells", serialize_ig_cells( tPrimaryCells ) );

            auto tVoidCells = tDataBase->get_void_cells_in_cluster();
            tClusterObj.put_child( "void_ig_cells", serialize_ig_cells( tVoidCells ) );

            const moris::mtk::Cell *tIPCell = &tDataBase->get_interpolation_cell();
            tClusterObj.put_child( "ip_cell", serialize_ip_cell( tIPCell ) );

            return tClusterObj;
        }

        Json Json_Debug_Output::serialize_ig_cells( Vector< moris::mtk::Cell const * > &aCells )
        {
            Json tPrimaryCells;
            for ( auto tCell : aCells )
            {
                for ( auto tVertex : tCell->get_vertex_pointers() )
                {
                    mAllIGVertices.insert( tVertex );
                }
                mAllIGCells.insert( tCell );
                Json tCellId;
                tCellId.put( "", tCell->get_id() );
                tPrimaryCells.push_back( { "", tCellId } );
            }
            return tPrimaryCells;
        }

        Json Json_Debug_Output::serialize_ip_cell( moris::mtk::Cell const *&aCell )
        {
            Json tPrimaryCells;
            for ( auto tVertex : aCell->get_vertex_pointers() )
            {
                mAllIPVertices.insert( tVertex );
            }
            mAllIPCells.insert( aCell );
            tPrimaryCells.put( "", aCell->get_id() );
            return tPrimaryCells;
        }

        Json Json_Debug_Output::serialize_side_sets( Vector< moris::mtk::Side_Set * > &aSideSets )
        {
            Json tSideSetsObj;
            for ( auto tSideSet : aSideSets )
            {
                Json tObj;
                auto tSideClusters = tSideSet->get_clusters_on_set();
                tObj.put_child( "clusters", serialize_side_clusters( tSideClusters ) );
                tObj.put( "index", tSideSet->get_set_index() );
                tObj.put( "name", tSideSet->get_set_name() );
                tSideSetsObj.push_back( { tSideSet->get_set_name(), tObj } );
            }
            return tSideSetsObj;
        }

        Json Json_Debug_Output::serialize_side_clusters( Vector< Cluster const * > &aSideClusters )
        {
            Json tSideClustersList;
            for ( auto tCluster : aSideClusters )
            {
                tSideClustersList.push_back( { "", serialize_side_cluster( tCluster ) } );
            }
            return tSideClustersList;
        }

        Json Json_Debug_Output::serialize_side_cluster( Cluster const *tCluster )
        {
            Json tClusterObj;
            auto tDataBase = dynamic_cast< Side_Cluster const * >( tCluster );
            for ( auto tVertex : tDataBase->get_vertices_in_cluster() )
            {
                mAllIGVertices.insert( tVertex );
            }

            auto tCells = tDataBase->get_cells_in_side_cluster();
            tClusterObj.put_child( "primary_ig_cells", serialize_ig_cells( tCells ) );

            const moris::mtk::Cell *tIPCell = &tDataBase->get_interpolation_cell();
            tClusterObj.put_child( "ip_cell", serialize_ip_cell( tIPCell ) );


            return tClusterObj;
        }
    }    // namespace mtk
}    // namespace moris