/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mesh_DataBase_IG.cpp
 * 
 */

#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Side_Cluster_DataBase.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_Tracer.hpp"

#include "cl_MTK_Mesh_DataBase_IP.hpp"

namespace moris::mtk
{

    //-----------------------------------------------------------------------------
    Integration_Mesh_DataBase_IG::Integration_Mesh_DataBase_IG( Integration_Mesh_DataBase* aIGDataBase, mtk::Integration_Mesh& aIGMesh, Interpolation_Mesh_Analysis& aIPMesh )                                                                                                                                                                        // mIGDataBase( aIGDataBase ), mIGMesh( aIGMesh ), mIPMesh( aIPMesh )
    {
    }

    //-----------------------------------------------------------------------------

    Integration_Mesh_DataBase_IG::~Integration_Mesh_DataBase_IG()
    {
        delete mCellClusterVertexCoords;
        delete mSecondaryClusterVertexCoords;
    }

    //-----------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_spatial_dim() const
    {
        return mVertexCoordinates.n_rows();
    }

    //-----------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_num_entities(
        enum EntityRank   aEntityRank,
        const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                return mVertices.size();
                break;
            }
            case EntityRank::ELEMENT:
            {
                return mCells.size();
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    //-----------------------------------------------------------------------------

    MeshType
    Integration_Mesh_DataBase_IG::get_mesh_type() const
    {
        return MeshType::MTK;
    }

    //-----------------------------------------------------------------------------

    Matrix< IdMat >
    Integration_Mesh_DataBase_IG::get_communication_table() const
    {
        MORIS_ERROR( 0, "get_communication_table not implemented for Integration_Mesh_DataBase_IG" );
        return { {} };
        // return mIGMesh.get_communication_table();
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Integration_Mesh_DataBase_IG::get_entity_connected_to_entity_loc_inds(
        moris_index       aEntityIndex,
        enum EntityRank   aInputEntityRank,
        enum EntityRank   aOutputEntityRank,
        const moris_index aDiscretizationIndex ) const
    {
        MORIS_ERROR( 0, "get_entity_connected_to_entity_loc_inds not implemented for Integration_Mesh_DataBase_IG" );
        return { {} };
        // return mIGMesh.get_entity_connected_to_entity_loc_inds( aEntityIndex, aInputEntityRank, aOutputEntityRank, aDiscretizationIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Integration_Mesh_DataBase_IG::get_node_coordinate( moris_index aNodeIndex ) const
    {
        return mVertexCoordinates.get_column( aNodeIndex );
    }

    // ----------------------------------------------------------------------------

    moris_id
    Integration_Mesh_DataBase_IG::get_glb_entity_id_from_entity_loc_index(
            moris_index aEntityIndex,
            EntityRank  aEntityRank,
            moris_index aDiscretizationIndex ) const
    {
        switch ( aEntityRank )
        {
            case ( EntityRank::NODE ):
                return mVertexIdList( aEntityIndex );
            case ( EntityRank::ELEMENT ):
                return mCellIdList( aEntityIndex );
            default:
                MORIS_ERROR( false, "Integration_Mesh_DataBase_IG::get_glb_entity_id_from_entity_loc_index() does not support given entity type" );
                return 0;
        }
    }

    // ----------------------------------------------------------------------------

    moris_index
    Integration_Mesh_DataBase_IG::get_loc_entity_ind_from_entity_glb_id(
            moris_id    aEntityId,
            EntityRank  aEntityRank,
            moris_index aDiscretizationIndex ) const
    {
        switch ( aEntityRank )
        {
            case ( EntityRank::NODE ):
                return mVertexGlobalIdToLocalIndex.find( aEntityId )->second;
            default:
                MORIS_ERROR( false, "Integration_Mesh_DataBase_IG::get_loc_entity_ind_from_entity_glb_id() does not support given entity type" );
                return 0;
        }
    }

    // ----------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_node_owner( moris_index aNodeIndex ) const
    {
        MORIS_ERROR( 0, "get_node_owner not implemented for Integration_Mesh_DataBase_IG" );
        return 0;
        // return mIGMesh.get_node_owner( aNodeIndex );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Integration_Mesh_DataBase_IG::get_element_indices_in_block_set( uint aSetIndex )
    {
        MORIS_ERROR( 0, "get_element_indices_in_block_set not implemented for Integration_Mesh_DataBase_IG" );
        return { {} };
        // return mIGMesh.get_element_indices_in_block_set( aSetIndex );
    }

    // ----------------------------------------------------------------------------

    enum CellTopology
    Integration_Mesh_DataBase_IG::get_blockset_topology( const std::string& aSetName )
    {
        // chceck if we have the same toplogy as the old mesh
        // MORIS_ASSERT( mCellTopologyToNameMap[aSetName] == mIGMesh.get_blockset_topology( aSetName ), "No the Same Cell Topo" );

        return mCellTopologyToNameMap[aSetName];
    }

    // ----------------------------------------------------------------------------
    enum CellShape
    Integration_Mesh_DataBase_IG::get_IG_blockset_shape( const std::string& aSetName )
    {
        // get the clusters in the set
        moris::Cell< mtk::Cluster const* > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        CellShape tCellShape = CellShape::EMPTY;

        // if the set isn't empty exist
        if ( tSetClusters.size() > 0 )
        {
            // get the cells in the first cluster
            moris::Cell< moris::mtk::Cell const* > tClusterCells = tSetClusters( 0 )->get_primary_cells_in_cluster();

            // compute the cell shape based on the first cell
            if( tClusterCells.size() > 0 )
            {
                tCellShape = tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
            }
            else // in case there are void clusters, look at the void cells
            {
                tClusterCells = tSetClusters( 0 )->get_void_cells_in_cluster();
                tCellShape =tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
            }
        }

// within debug, check all cells to make sure that they are the same Cell Shape
#ifdef MORIS_HAVE_DEBUG

        // skip check for sets only used for visualization purposes
        if( !std::strstr( aSetName.c_str(), "Vis" ) )
        {
            // checking all clusters in set
            for ( uint iCluster = 0; iCluster < tSetClusters.size(); iCluster++ )
            {
                // get cell of cells in the cluster
                moris::Cell< moris::mtk::Cell const* > tClusterCellsCheck = tSetClusters( iCluster )->get_primary_cells_in_cluster();

                // looping through the cells in the cluster
                for ( uint iCheckCell = 0; iCheckCell < tClusterCellsCheck.size(); iCheckCell++ )
                {
                    MORIS_ASSERT( tClusterCellsCheck( iCheckCell )->get_cell_info()->compute_cell_shape( tClusterCellsCheck( iCheckCell ) ) == tCellShape,
                        "Integration_Mesh_DataBase_IG::get_IG_blockset_shape - cell shape is not consistent in the block" );
                }
            }
        }
#endif

        // return the shape of the IG cells
        return tCellShape;
    }

    // ----------------------------------------------------------------------------

    enum CellShape
    Integration_Mesh_DataBase_IG::get_IP_blockset_shape( const std::string& aSetName )
    {
        // MORIS_ERROR( 0, "get_IP_blockset_shape not implemented for Integration_Mesh_DataBase_IG" );
        // return mIGMesh.get_IP_blockset_shape( aSetName );

        // get the clusters in the set
        moris::Cell< mtk::Cluster const* > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        CellShape tCellShape = CellShape::EMPTY;

        // if the set isn't empty exist
        if ( tSetClusters.size() > 0 )
        {
            // get the cells in the first cluster
            mtk::Cell const& tClusterCell = tSetClusters( 0 )->get_interpolation_cell();

            // compute the cell shape of the first cell
            tCellShape = tClusterCell.get_cell_info()->compute_cell_shape( &tClusterCell );
        }

        // within debug, checking all cells to make sure that they are the same Cell Shape
        // if cells exist
        // looping through the clusters
        for ( uint iCluster = 1; iCluster < tSetClusters.size(); iCluster++ )
        {
            MORIS_ASSERT( tSetClusters( iCluster )->get_interpolation_cell().get_cell_info()->compute_cell_shape( &tSetClusters( iCluster )->get_interpolation_cell() ) == tCellShape,
                "Enriched_Integration_Mesh::get_IP_blockset_shape - cell shape is not consistent in the block" );
        }

        return tCellShape;
    }

    // ----------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_element_owner( moris_index aElementIndex ) const
    {
        MORIS_ERROR( 0, "get_element_owner not implemented for Integration_Mesh_DataBase_IG" );
        return 0;
        // return mIGMesh.get_element_owner( aElementIndex );
    }

    // ----------------------------------------------------------------------------
    std::unordered_map< moris_id, moris_index >
    Integration_Mesh_DataBase_IG::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        return mVertexGlobalIdToLocalIndex;
    }

    // ----------------------------------------------------------------------------
    Vertex&
    Integration_Mesh_DataBase_IG::get_mtk_vertex( moris_index aVertexIndex )
    {
        // check vertex indices bound
        MORIS_ASSERT( aVertexIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );

        // return the vertex
        return mVertices( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vertex const&
    Integration_Mesh_DataBase_IG::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        // check vertex indices bound
        MORIS_ASSERT( aVertexIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );

        // return the vertex
        return mVertices( aVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< IndexMat >
    Integration_Mesh_DataBase_IG::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
    {
        MORIS_ERROR( 0, "get_elements_connected_to_element_and_face_ind_loc_inds Not implemented for the mesh class" );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell_Cluster const&
    Integration_Mesh_DataBase_IG::get_cell_cluster( Cell const& aInterpCell ) const
    {
        MORIS_ASSERT( aInterpCell.get_index() < (moris_index)mCellClusters.size(), "Interpolation cell index is out of bounds." );
        return mCellClusters( aInterpCell.get_index() );
    }

    // ----------------------------------------------------------------------------

    Cell_Cluster const&
    Integration_Mesh_DataBase_IG::get_cell_cluster( moris_index aInterpCellIndex ) const
    {
        MORIS_ASSERT( aInterpCellIndex < (moris_index)mCellClusters.size(), "Interpolation Cell index out of bounds" );
        return mCellClusters( aInterpCellIndex );
    }


    // ----------------------------------------------------------------------------

    moris::Cell< std::string >
    Integration_Mesh_DataBase_IG::get_block_set_names() const
    {
        MORIS_ERROR( 0, "get_block_set_names not implemented for Integration_Mesh_DataBase_IG" );
        return { "" };
        // return mIGMesh.get_block_set_names();
    }

    // ----------------------------------------------------------------------------

    moris::Cell< Cluster const* >
    Integration_Mesh_DataBase_IG::get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const
    {
        MORIS_ERROR( 0, "get_cell_clusters_in_set not implemented for Integration_Mesh_DataBase_IG" );
        return { nullptr };
        // return mIGMesh.get_cell_clusters_in_set( aBlockSetOrdinal );
    }

    // ----------------------------------------------------------------------------

    moris::Cell< Cluster const* >
    Integration_Mesh_DataBase_IG::get_side_set_cluster( moris_index aSideSetOrdinal ) const
    {
        MORIS_ERROR( 0, "get_side_set_cluster not implemented for Integration_Mesh_DataBase_IG" );
        return { nullptr };
        // return mIGMesh.get_side_set_cluster( aSideSetOrdinal );
    }

    // ----------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_num_side_sets() const
    {
        MORIS_ERROR( 0, "get_num_side_sets not implemented for Integration_Mesh_DataBase_IG" );
        return 0;
        // return mIGMesh.get_num_side_sets();
    }


    // ----------------------------------------------------------------------------

    std::string
    Integration_Mesh_DataBase_IG::get_side_set_label( moris_index aSideSetOrdinal ) const
    {
        MORIS_ERROR( 0, "get_side_set_label not implemented for Integration_Mesh_DataBase_IG" );
        return "";
        // return mIGMesh.get_side_set_label( aSideSetOrdinal );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Integration_Mesh_DataBase_IG::get_side_set_index( std::string aSideSetLabel ) const
    {
        MORIS_ERROR( 0, "get_side_set_index not implemented for Integration_Mesh_DataBase_IG" );
        return 0;
        // return mIGMesh.get_side_set_index( aSideSetLabel );
    }

    // ----------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_num_double_sided_sets() const
    {
        MORIS_ERROR( 0, "get_num_double_sided_sets not implemented for Integration_Mesh_DataBase_IG" );
        return 0;
        // return mIGMesh.get_num_double_sided_sets();
    }

    // ----------------------------------------------------------------------------

    std::string
    Integration_Mesh_DataBase_IG::get_double_sided_set_label( moris_index aSideSetOrdinal ) const
    {
        MORIS_ERROR( 0, "get_double_sided_set_label not implemented for Integration_Mesh_DataBase_IG" );
        return "";
        // return mIGMesh.get_double_sided_set_label( aSideSetOrdinal );
    }

    // ----------------------------------------------------------------------------

    moris::Cell< Cluster const* >
    Integration_Mesh_DataBase_IG::get_double_side_set_cluster( moris_index aSideSetOrdinal ) const
    {
        MORIS_ERROR( 0, "get_double_side_set_cluster not implemented for Integration_Mesh_DataBase_IG" );
        return { nullptr };
        // return mIGMesh.get_double_side_set_cluster( aSideSetOrdinal );
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Cell&
    Integration_Mesh_DataBase_IG::get_mtk_cell( moris_index aCellIndex )
    {
        MORIS_ASSERT( aCellIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
        return mCells( aCellIndex );
    }

    // ----------------------------------------------------------------------------

    mtk::Cell const&
    Integration_Mesh_DataBase_IG::get_mtk_cell( moris_index aCellIndex ) const
    {
        MORIS_ASSERT( aCellIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
        return mCells( aCellIndex );
    }

    // ----------------------------------------------------------------------------

    moris_id
    Integration_Mesh_DataBase_IG::get_entity_id( enum EntityRank aEntityRank, moris_index aEntityIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );
                return mVertexIdList( aEntityIndex );
                break;
            }
            case EntityRank::ELEMENT:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
                return mCellIdList( aEntityIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id 
    Integration_Mesh_DataBase_IG::get_entity_owner( enum EntityRank aEntityRank, moris_index aEntityIndex ) const
    {
        switch ( aEntityRank )
        {
            case EntityRank::NODE:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mVertices.size(), "index of the vertex specified exceeds the bounds" );
                return mVertexOwnerList( aEntityIndex );
                break;
            }
            case EntityRank::ELEMENT:
            {
                MORIS_ASSERT( aEntityIndex < (moris_index)mCells.size(), "index of the vertex specified exceeds the bounds" );
                return mCellOwnerList( aEntityIndex );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Cell_Info >
    Integration_Mesh_DataBase_IG::get_cell_info_sp( moris_index aEntityIndex ) const
    {
        return mCellInfoList( aEntityIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::real*
    Integration_Mesh_DataBase_IG::get_vertex_coords_ptr( moris_index aVertexIndex )
    {
        moris::real* tCoordsPointer = const_cast< real* >( mVertexCoordinates.colptr( aVertexIndex ) );
        return tCoordsPointer;
    }

    //--------------------------------------------------------------------------------------------------------------


    Vertex**
    Integration_Mesh_DataBase_IG::get_cell_vertices( moris_index aCellIndex )
    {
        return mCellToVertices.memptr() + mCellToVertexOffSet( aCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Cell* 
    Integration_Mesh_DataBase_IG::get_ip_cell_in_cluster( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)mCellClusters.size(), "index of the vertex specified exceeds the bounds" );
                return &mIPMesh->get_mtk_cell( aClusterIndex );
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return &mIPMesh->get_mtk_cell( mSideClusterToIPCell( aClusterIndex ) );
                }
                else
                {
                    return &mIPMesh->get_mtk_cell( mGhostLeaderFollowerIPCellList( aClusterIndex - mSideClusters.size() ) );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Cell* const*
    Integration_Mesh_DataBase_IG::get_ig_cells_in_cluster( enum ClusterType aClusterType, Primary_Void aPrimaryOrVoid, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                switch ( aPrimaryOrVoid )
                {
                    case mtk::Primary_Void::PRIMARY:
                    {
                        return mCellClusterToPrimaryIGCell.memptr() + mCellClusterToPrimaryIGCellOffSet( aClusterIndex );
                        break;
                    }
                    case mtk::Primary_Void::VOID:
                    {
                        return mCellClusterToVoidIGCell.memptr() + mCellClusterToVoidIGCellOffset( aClusterIndex );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
                    }
                        return nullptr;
                }
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return mSideClusterToPrimaryIGCell.memptr() + mSideClusterToPrimaryIGCellOffset( aClusterIndex );
                }
                else
                {
                    return mGhostLeaderFollowerIGCellList.memptr() + aClusterIndex - mSideClusters.size();
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint 
    Integration_Mesh_DataBase_IG::get_num_cells_in_cluster( enum ClusterType aClusterType, Primary_Void aPrimaryOrVoid, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                switch ( aPrimaryOrVoid )
                {
                    case mtk::Primary_Void::PRIMARY:
                    {
                        return mCellClusterToPrimaryIGCellOffSet( aClusterIndex + 1 ) - mCellClusterToPrimaryIGCellOffSet( aClusterIndex );
                        break;
                    }
                    case mtk::Primary_Void::VOID:
                    {
                        return mCellClusterToVoidIGCellOffset( aClusterIndex + 1 ) - mCellClusterToVoidIGCellOffset( aClusterIndex );
                        break;
                    }
                    default:
                    {
                        MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
                    }

                        return 0;
                }
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return mSideClusterToPrimaryIGCellOffset( aClusterIndex + 1 ) - mSideClusterToPrimaryIGCellOffset( aClusterIndex );
                }
                else
                {
                    return 1;
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index*
    Integration_Mesh_DataBase_IG::get_side_ordinals_in_cluster( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        if ( aClusterIndex < (moris_index)mSideClusters.size() )
        {
            return const_cast< moris_index* >( mSideClusterToPrimaryIGCellSideOrd.memptr() + mSideClusterToPrimaryIGCellOffset( aClusterIndex ) );
        }
        else
        {
            return const_cast< moris_index* >( mGhostLeaderFollowerOrd.memptr() + aClusterIndex - mSideClusters.size() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_DataBase_IG::cluster_is_trivial( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                return mCellClusterIsTrivial( aClusterIndex );
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return mSideClusterIsTrivial( aClusterIndex );
                }
                else
                {
                    return mGhostLeaderFollowerIsTrivial( aClusterIndex - mSideClusters.size() );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return false;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vertex* const*
    Integration_Mesh_DataBase_IG::get_vertices_in_cluster( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                return mCellClusterToVeretx.memptr() + mCellClusterToVertexOffset( aClusterIndex );
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return mSideClusterToVeretx.memptr() + mSideClusterToVertexOffSet( aClusterIndex );
                }
                else
                {
                    return mGhostLeaderFollowerToVertex.memptr() + mGhostLeaderFollowerVertexOffSet( aClusterIndex - mSideClusters.size() );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_num_vertices_in_cluster( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                return mCellClusterToVertexOffset( aClusterIndex + 1 ) - mCellClusterToVertexOffset( aClusterIndex );
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return mSideClusterToVertexOffSet( aClusterIndex + 1 ) - mSideClusterToVertexOffSet( aClusterIndex );
                }
                else
                {
                    return mGhostLeaderFollowerVertexOffSet( aClusterIndex - mSideClusters.size() + 1 ) - mGhostLeaderFollowerVertexOffSet( aClusterIndex - mSideClusters.size() );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >*
    Integration_Mesh_DataBase_IG::get_local_coord_matrix_ptr( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                return mCellClusterVertexCoords;
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    return mCellClusterVertexCoords;
                }
                else
                {
                    return mSecondaryClusterVertexCoords;
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Integration_Mesh_DataBase_IG::get_row_number_local_coords_matrix( enum ClusterType aClusterType, moris_index aClusterIndex ) const
    {
        switch ( aClusterType )
        {
            case ClusterType::CELL:
            {
                return mCellClusterIndexToRowNumber.at( aClusterIndex );
                break;
            }
            case ClusterType::SIDE:
            {
                if ( aClusterIndex < (moris_index)mSideClusters.size() )
                {
                    if ( aClusterIndex < (moris_index)mNumSideClusters )
                    {       
                        return mCellClusterIndexToRowNumber.at( mSideClusterToIPCell( aClusterIndex ) );
                    }    
                    else
                    {
                         return mSideClusterIndexToRowNumber.at( aClusterIndex );
                    }
                       
                }
                else
                {
                    return mSecondaryClusterIndexToRowNumber.at( aClusterIndex - mSideClusters.size() );
                    // return mSecondaryClusterIndexToRowNumber[aClusterIndex - mSideClusters.size()];
                }
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }

                return 0;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Cell_Cluster const*
    Integration_Mesh_DataBase_IG::get_associated_cell_cluster( moris_index aClusterIndex ) const
    {
        if ( aClusterIndex < (moris_index)mSideClusters.size() )
        {
            return &mCellClusters( mSideClusterToIPCell( aClusterIndex ) );
        }
        else
        {
            return &mCellClusters( mGhostLeaderFollowerIPCellList( aClusterIndex ) );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_DataBase_IG::is_secondary_cluster( moris_index aClusterIndex ) const
    {
        if ( aClusterIndex < (moris_index)mNumSideClusters )
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::Memory_Map
    Integration_Mesh_DataBase_IG::get_memory_usage()
    {
        moris::Memory_Map tMemoryMap;

        // tMemoryMap.mMemoryMapData.reserve( 4 );

        // // pointers and references
        // tMemoryMap.mMemoryMapData["misc"] = sizeof( mIGDataBase ) * 3;

        // // objects,  vertices, cells, ...
        // tMemoryMap.mMemoryMapData["objects"] = moris::internal_capacity( mVertices )
        //                                        + moris::internal_capacity( mCells )
        //                                        + moris::internal_capacity( mCellClusters )
        //                                        + moris::internal_capacity( mSideClusters )
        //                                        + moris::internal_capacity( mDblSideClusters )
        //                                        + moris::internal_capacity( mGhostLeader )
        //                                        + moris::internal_capacity( mGhostFollower )
        //                                        + moris::internal_capacity( mGhostDblSidedSet );

        // // connectivity pointers
        // tMemoryMap.mMemoryMapData["connectivity"] = mCellToVertices.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mCellClusterToPrimaryIGCell.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mCellClusterToVoidIGCell.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mCellClusterToVeretx.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mSideClusterToPrimaryIGCell.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mSideClusterToVoidIGCell.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mSideClusterToVeretx.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mGhostLeaderFollowerIGCellList.capacity() * ( 1 + sizeof( void* ) )
        //                                             + mGhostLeaderFollowerToVertex.capacity() * ( 1 + sizeof( void* ) );

        // // listed memebr data
        // tMemoryMap.mMemoryMapData["Indexed cells"] = mVertexIdList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mCellIdList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mVertexOwnerList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mVertexIdList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mCellOwnerList.capacity() * ( 1 + sizeof( moris_id ) )
        //                                              + mGhostLeaderFollowerIPCellList.capacity() * ( 1 + sizeof( moris_index ) )
        //                                              + mCellClusterIsTrivial.capacity() * ( 1 + sizeof( bool ) )
        //                                              + mCellClusterIsTrivial.capacity() * ( 1 + sizeof( bool ) );

        // tMemoryMap.mMemoryMapData["Sets"] = moris::internal_capacity_ptr( mListOfAllSets )
        //                                     + mListOfBlocks.capacity() * ( 1 + sizeof( void* ) )
        //                                     + mListOfSideSets.capacity() * ( 1 + sizeof( void* ) )
        //                                     + mListOfDoubleSideSets.capacity() * ( 1 + sizeof( void* ) );

        // // vertex global to local map size
        // size_t tVertexGlbToLocalMapSize = 0;

        // uint  n = mVertexGlobalIdToLocalIndex.bucket_count();
        // float m = mVertexGlobalIdToLocalIndex.max_load_factor();
        // if ( m > 1.0 )
        // {
        //     tVertexGlbToLocalMapSize = (size_t)n * m;
        // }
        // else
        // {
        //     tVertexGlbToLocalMapSize = (size_t)n;
        // }

        // // cell toplogy map data
        // size_t tDofMapSize = ( sizeof( enum CellTopology ) + sizeof( std::string ) + sizeof( std::_Rb_tree_node_base ) ) * mCellTopologyToNameMap.size();//+ sizeof(std::_Rb_tree);

        // // add up the map data
        // tMemoryMap.mMemoryMapData["Map"] = tVertexGlbToLocalMapSize + tDofMapSize;

        // // IG data base memory map
        // // moris::Memory_Map tIGDataBaseMM;

        // // tIGDataBaseMM = mIGDataBase->get_memory_usage();

        // // tMemoryMap.mMemoryMapData["database"] = tIGDataBaseMM.sum();

        // // // print the memory
        // // tMemoryMap.par_print( "IG-Mesh" );

        return tMemoryMap;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Integration_Mesh_DataBase_IG::free_memory()
    {
    }


}// namespace moris::mtk
