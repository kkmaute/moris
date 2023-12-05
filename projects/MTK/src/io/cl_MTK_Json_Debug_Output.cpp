//
// Created by frank on 12/4/23.
//

#include <unordered_map>
#include "cl_MTK_Json_Debug_Output.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Side_Cluster_DataBase.hpp"
#include "iostream"

using json = nlohmann::json;


namespace moris
{
    namespace mtk
    {
        void Json_Debug_Output::write_to_json( std::string const &aFileName )
        {
            json tMeshObj;
            mAllIGVertices.clear();
            mAllIPVertices.clear();
            mAllIGCells.clear();
            mAllIPCells.clear();

            // get all block sets
            auto tBlockSets                    = mMesh->get_block_sets();
            tMeshObj[ "sets" ][ "block_sets" ] = serialize_block_sets( tBlockSets );

            // side sets
            auto tSideSets                    = mMesh->get_side_sets();
            tMeshObj[ "sets" ][ "side_sets" ] = serialize_side_sets( tSideSets );

            // double side sets
            auto tDoubleSideSets                     = mMesh->get_double_side_sets();
            tMeshObj[ "sets" ][ "double_side_sets" ] = serialize_double_side_sets( tDoubleSideSets );

            // store vertices
            tMeshObj[ "vertices_ig" ] = serialize_all_vertices( mAllIGVertices );
            tMeshObj[ "vertices_ip" ] = serialize_all_vertices( mAllIPVertices );
            tMeshObj[ "cells_ig" ]    = serialize_all_cells( mAllIGCells );
            tMeshObj[ "cells_ip" ]    = serialize_all_cells( mAllIPCells );

            std::ofstream tFile( aFileName );
            //            std::cout << to_string( tMeshObj );
            tFile << to_string( tMeshObj );
            tFile.close();
        }

        json Json_Debug_Output::serialize_all_vertices( const std::unordered_set< Vertex const * > &aVertices )
        {
            json tVertexObj;

            json::array_t tIds;
            json::array_t tIndices;
            json::array_t tX;
            json::array_t tY;
            json::array_t tZ;

            for ( auto tVertexPair : aVertices )
            {
                auto tDataBase = dynamic_cast< Vertex_DataBase const * >( tVertexPair );

                tIds.emplace_back( tDataBase->get_id() );
                tIndices.emplace_back( tDataBase->get_index() );
                tX.emplace_back( tDataBase->get_coords()( 0 ) );
                tY.emplace_back( tDataBase->get_coords()( 1 ) );
                if ( mMesh->get_spatial_dim() == 3 )
                {
                    tZ.emplace_back( tDataBase->get_coords()( 2 ) );
                }
            }
            tVertexObj[ "id" ]    = tIds;
            tVertexObj[ "index" ] = tIndices;
            tVertexObj[ "x" ]     = tX;
            tVertexObj[ "y" ]     = tY;
            if ( mMesh->get_spatial_dim() == 3 )
            {
                tVertexObj[ "z" ] = tZ;
            }
            return tVertexObj;
        }

        json::array_t Json_Debug_Output::serialize_all_cells( const std::unordered_set< Cell const * > &aCells )
        {
            json::array_t tAllCells;

            for ( auto tCellPair : aCells )
            {
                auto tDataBase = dynamic_cast< Cell_DataBase const * >( tCellPair );
                if ( tDataBase == nullptr )
                {
                    std::cout << "Cell with id " << tCellPair->get_id() << " is not a cell database" << std::endl;
                    continue;
                }
                json     tCellObj;
                moris_id tId         = tDataBase->get_id();
                tCellObj[ "id" ]     = tId;
                tCellObj[ "index" ]  = tDataBase->get_index();
                tCellObj[ "facets" ] = json::array();

                for ( uint iSide = 0; iSide < tDataBase->get_number_of_facets(); ++iSide )
                {
                    json tSideObj;
                    tSideObj[ "ordinal" ]    = iSide;
                    tSideObj[ "vertex_ids" ] = json::array();
                    //                    tSideObj[ "vertex_indices" ] = json::array();
                    auto tSideVertices = tDataBase->get_vertices_on_side_ordinal( iSide );
                    for ( auto tVertex : tSideVertices )
                    {
                        tSideObj[ "vertex_ids" ].push_back( tVertex->get_id() );
                        //                        tSideObj[ "vertex_indices" ].push_back( tVertex->get_index() );
                    }
                    tCellObj[ "facets" ].push_back( tSideObj );
                }
                tAllCells.push_back( tCellObj );
            }
            return tAllCells;
        }

        json Json_Debug_Output::serialize_double_side_sets( moris::Cell< moris::mtk::Double_Side_Set * > &aDoubleSideSets )
        {
            json tDoubleSideSetsObj;
            for ( auto tSideSet : aDoubleSideSets )
            {
                json tSideSetObj;
                tSideSetObj[ "name" ]     = tSideSet->get_set_name();
                tSideSetObj[ "index" ]    = tSideSet->get_set_index();
                tSideSetObj[ "clusters" ] = json::array();
                for ( auto tCluster : tSideSet->get_clusters_on_set() )
                {
                    json tClusterObj;
                    tClusterObj[ "follower_cluster" ] = serialize_side_cluster( &tCluster->get_follower_side_cluster() );
                    tClusterObj[ "leader_cluster" ]   = serialize_side_cluster( &tCluster->get_leader_side_cluster() );
                }

                tDoubleSideSetsObj[ tSideSet->get_set_name() ] = tSideSetObj;
            }
            return tDoubleSideSetsObj;
        }

        json Json_Debug_Output::serialize_block_sets( moris::Cell< moris::mtk::Block_Set * > &tBlockSets )
        {
            json tBlockSetsObj;
            for ( auto tBlockSet : tBlockSets )
            {
                json tBlockSetObj;
                tBlockSetObj[ "name" ]  = tBlockSet->get_set_name();
                tBlockSetObj[ "index" ] = tBlockSet->get_set_index();

                auto tClusters             = tBlockSet->get_clusters_on_set();
                tBlockSetObj[ "clusters" ] = serialize_cell_clusters( tClusters );

                tBlockSetsObj[ tBlockSet->get_set_name() ] = tBlockSetObj;
            }
            return tBlockSetsObj;
        }

        json::array_t Json_Debug_Output::serialize_cell_clusters( moris::Cell< Cluster const * > &tClusters )
        {
            json::array_t tClusterArr;
            for ( auto tCluster : tClusters )
            {
                tClusterArr.push_back( serialize_cell_cluster( tCluster ) );
            }
            return tClusterArr;
        }

        json Json_Debug_Output::serialize_cell_cluster( Cluster const *aCluster )
        {
            json tClusterObj;
            auto tDataBase                    = dynamic_cast< Cell_Cluster_DataBase const                    *>( aCluster );
            auto tPrimaryCells                = tDataBase->get_primary_cells_in_cluster();
            tClusterObj[ "primary_ig_cells" ] = serialize_cells( tPrimaryCells );
            auto tVoidCells                   = tDataBase->get_void_cells_in_cluster();
            tClusterObj[ "void_ig_cells" ]    = serialize_cells( tVoidCells );
            moris::Cell< const moris::mtk::Cell * > tIPCells{ &tDataBase->get_interpolation_cell() };
            tClusterObj[ "ip_cell" ] = serialize_cells( tIPCells, true )[ 0 ];
            return tClusterObj;
        }

        json::array_t Json_Debug_Output::serialize_cells( moris::Cell< moris::mtk::Cell const * > &aCells, bool aIsIPCell )
        {
            json::array_t tPrimaryCells;
            for ( auto tCell : aCells )
            {
                for ( auto tVertex : tCell->get_vertex_pointers() )
                {
                    if ( aIsIPCell )
                    {
                        mAllIPVertices.insert( tVertex );
                    }
                    else
                    {
                        mAllIGVertices.insert( tVertex );
                    }
                }
                if ( aIsIPCell )
                {
                    mAllIPCells.insert( tCell );
                }
                else
                {
                    mAllIGCells.insert( tCell );
                }
                tPrimaryCells.push_back( tCell->get_id() );
            }
            return tPrimaryCells;
        }

        json Json_Debug_Output::serialize_side_sets( moris::Cell< moris::mtk::Side_Set * > &aSideSets )
        {
            json tSideSetsObj;
            for ( auto tSideSet : aSideSets )
            {
                json tObj;
                auto tSideClusters = tSideSet->get_clusters_on_set();
                tObj[ "clusters" ] = serialize_side_clusters( tSideClusters );
                tObj[ "index" ]    = tSideSet->get_set_index();
                tObj[ "name" ]     = tSideSet->get_set_name();


                tSideSetsObj[ tSideSet->get_set_name() ] = tObj;
            }
            return tSideSetsObj;
        }

        json::array_t Json_Debug_Output::serialize_side_clusters( moris::Cell< Cluster const * > &aSideClusters )
        {
            json::array_t tSideClustersList;
            for ( auto tCluster : aSideClusters )
            {
                tSideClustersList.push_back( serialize_side_cluster( tCluster ) );
            }
            return tSideClustersList;
        }

        json Json_Debug_Output::serialize_side_cluster( Cluster const *tCluster )
        {
            json tClusterObj;
            auto tDataBase = dynamic_cast< Side_Cluster_DataBase const * >( tCluster );
            for ( auto tVertex : tDataBase->get_vertices_in_cluster() )
            {
                mAllIGVertices.insert( tVertex );
                auto tCells                       = tDataBase->get_cells_in_side_cluster();
                tClusterObj[ "primary_ig_cells" ] = serialize_cells( tCells );
                moris::Cell< const moris::mtk::Cell * > tIPCells{ &tDataBase->get_interpolation_cell() };
                tClusterObj[ "ip_cell" ] = serialize_cells( tIPCells, true )[ 0 ];
            }
            return tClusterObj;
        }

    }    // namespace mtk
}    // namespace moris