/**
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Enriched_Integration_Mesh.cpp
 *
 */

#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Side_Cluster.hpp"
#include "cl_XTK_Cell_Cluster_Group.hpp"
#include "cl_XTK_Side_Cluster_Group.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Set_Communicator.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "fn_isempty.hpp"
#include "cl_TOL_Memory_Map.hpp"
#include <memory>
#include "cl_Logger.hpp"
#include "fn_XTK_match_normal_to_side_ordinal.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris::xtk
{
    //------------------------------------------------------------------------------

    Enriched_Integration_Mesh::Enriched_Integration_Mesh( Model *aXTKModel )
            : mModel( aXTKModel )
            , mCutIgMesh( mModel->get_cut_integration_mesh() )
            , mCellClusters( 0, nullptr )
            , mFields( 0 )
            , mGlobalSetFieldLabelToIndex( 3 )
            , mSetWiseFieldLabelToIndex( 3 )
            , mCellInfo( nullptr )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Enriched Integration Mesh", "Construction", mModel->mVerboseLevel, 0 );

        // standard assumption: Enr. IG Mesh is associated with first B-spline mesh
        mBsplineMeshIndices = { { 0 } };

        this->setup_cell_clusters();

        this->setup_blockset_with_cell_clusters();

        this->setup_side_set_clusters();

        this->setup_double_side_set_clusters();

        this->setup_interface_side_sets();

        // this->setup_interface_vertex_sets();

        this->setup_color_to_set();

        this->collect_all_sets();

        // get the Cell info for trivial integration clusters
        if ( this->get_spatial_dim() == 2 )
        {
            mCellInfo = new moris::mtk::Cell_Info_Quad4();
        }
        else if ( this->get_spatial_dim() == 3 )
        {
            mCellInfo = new moris::mtk::Cell_Info_Hex8();
        }
    }

    //------------------------------------------------------------------------------

    Enriched_Integration_Mesh::Enriched_Integration_Mesh(
            Model                   *aXTKModel,
            const Matrix< IndexMat > aBsplineMeshIndices )
            : mModel( aXTKModel )
            , mCutIgMesh( mModel->get_cut_integration_mesh() )
            , mCellClusters( 0, nullptr )
            , mFields( 0 )
            , mGlobalSetFieldLabelToIndex( 3 )
            , mSetWiseFieldLabelToIndex( 3 )
            , mCellInfo( nullptr )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Enriched Integration Mesh", "SPG Based Construction", mModel->mVerboseLevel, 0 );

        // copy list of mesh indices
        mBsplineMeshIndices = aBsplineMeshIndices;

        this->setup_cell_clusters_new();

        this->setup_blockset_with_cell_clusters();

        this->setup_side_set_clusters();

        this->setup_double_side_set_clusters();

        this->setup_interface_side_sets();

        this->setup_color_to_set();

        this->collect_all_sets();

        // get the Cell info for trivial integration clusters
        if ( this->get_spatial_dim() == 2 )
        {
            mCellInfo = new moris::mtk::Cell_Info_Quad4();
        }
        else if ( this->get_spatial_dim() == 3 )
        {
            mCellInfo = new moris::mtk::Cell_Info_Hex8();
        }
    }

    //------------------------------------------------------------------------------

    Enriched_Integration_Mesh::~Enriched_Integration_Mesh()
    {
        delete mCellInfo;

        for ( auto p : mListOfBlocks )
        {
            delete p;
        }
        mListOfBlocks.clear();

        for ( auto p : mListOfSideSets )
        {
            delete p;
        }
        mListOfSideSets.clear();

        for ( auto p : mListOfDoubleSideSets )
        {
            delete p;
        }
        mListOfDoubleSideSets.clear();

        // for ( auto p : mDoubleSideClusters )
        // {
        //     delete p;
        // }
        // mDoubleSideClusters.clear();

        mDoubleSideSingleSideClusters.clear();
    }

    //------------------------------------------------------------------------------

    mtk::MeshType
    Enriched_Integration_Mesh::get_mesh_type() const
    {
        return mtk::MeshType::XTK;
    }

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_mesh_index_map()
    {
        // setup conversion table from MSI discretization mesh index to position of B-spline mesh in list of XTK meshes
        if ( !mLocalMeshIndexMapIsSet )
        {
            this->setup_mesh_index_map();
            mLocalMeshIndexMapIsSet = true;
        }

        // find mesh list index in map
        for ( uint iMeshIndex = 0; iMeshIndex < (uint)mBsplineMeshIndices.max() + 1; iMeshIndex++ )
        {
            mMeshIndexToLocMeshIndex[ iMeshIndex ] = iMeshIndex;
        }
    }

    // ----------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::get_local_mesh_index_xtk( moris_index const & aDiscretizationMeshIndex ) const
    {
        auto tIter = mMeshIndexToLocMeshIndex.find( aDiscretizationMeshIndex );

        MORIS_ASSERT( tIter != mMeshIndexToLocMeshIndex.end(), "Enriched_Integration_Mesh::get_local_mesh_index_xtk() - Mesh index not in map" );

        return tIter->second;
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_enriched_mesh_indices() const
    {
        return mBsplineMeshIndices;
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_interpolation_types() const
    {
        return mBsplineMeshIndices.numel();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_cell_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
    {
        // get the list index for the discretization mesh index
        moris_index tBsplineMeshListIndex = this->get_local_mesh_index_xtk( aDiscretizationMeshIndex );

        // get the number of cluster groups for the current B-spline mesh index
        return mCellClusterGroups( tBsplineMeshListIndex ).size();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
    {
        // get the list index for the discretization mesh index
        moris_index tBsplineMeshListIndex = this->get_local_mesh_index_xtk( aDiscretizationMeshIndex );

        // get the number of cluster groups for the current B-spline mesh index
        return mDblSideClusterGroups( tBsplineMeshListIndex ).size();
    }

    // ----------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_dbl_side_single_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
    {
        // get the list index for the discretization mesh index
        moris_index tBsplineMeshListIndex = this->get_local_mesh_index_xtk( aDiscretizationMeshIndex );

        // get the number of cluster groups for the current B-spline mesh index
        return mDblSideClusterGroups( tBsplineMeshListIndex ).size();
    }

    // ----------------------------------------------------------------------------

    Vector< std::shared_ptr< xtk::Cell_Cluster_Group > > const &
    Enriched_Integration_Mesh::get_cell_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
    {
        // get the list index for the discretization mesh index
        moris_index tBsplineMeshListIndex = this->get_local_mesh_index_xtk( aDiscretizationMeshIndex );

        // get the number of cluster groups for the current B-spline mesh index
        return mCellClusterGroups( tBsplineMeshListIndex );
    }

    // ----------------------------------------------------------------------------

    Vector< std::shared_ptr< xtk::Side_Cluster_Group > > const &
    Enriched_Integration_Mesh::get_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
    {
        // get the list index for the discretization mesh index
        moris_index tBsplineMeshListIndex = this->get_local_mesh_index_xtk( aDiscretizationMeshIndex );

        // get the number of cluster groups for the current B-spline mesh index
        return mDblSideClusterGroups( tBsplineMeshListIndex );
    }

    // ----------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_spatial_dim() const
    {
        return mCutIgMesh->get_spatial_dim();
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_entities( mtk::EntityRank aEntityRank, const moris_index aIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mCutIgMesh->get_num_entities( mtk::EntityRank::NODE, 0 );
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                return mCutIgMesh->get_num_entities( mtk::EntityRank::ELEMENT, 0 );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Only support get num entities for nodes and elements currently" );
            }
                return 0;
        }
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_owned_cells() const
    {
        uint tNumEntities = this->get_num_entities( mtk::EntityRank::ELEMENT );

        uint tNumOwnedEntities = 0;

        moris_id tParRank = moris::par_rank();

        // iterate and find out how many I own
        for ( uint i = 0; i < tNumEntities; i++ )
        {
            mtk::Cell const &tCell = this->get_mtk_cell( (moris_index)i );
            if ( tCell.get_owner() == tParRank )
            {
                tNumOwnedEntities++;
            }
        }

        return tNumOwnedEntities;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index       aEntityIndex,
            mtk::EntityRank   aInputEntityRank,
            mtk::EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        MORIS_ERROR( aInputEntityRank == mtk::EntityRank::ELEMENT && aOutputEntityRank == mtk::EntityRank::NODE,
                "Only support element to node connectivity" );

        return this->get_mtk_cell( aEntityIndex ).get_vertex_inds();
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
    {
        MORIS_ERROR( 0,
                "XTK ENRICHED MESH ERROR: get_elements_connected_to_element_and_face_ind_loc_inds no implemented" );

        return Matrix< IndexMat >( 0, 0 );
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Vertex const * >
    Enriched_Integration_Mesh::get_all_vertices() const
    {
        uint                          tNumNodes = this->get_num_entities( mtk::EntityRank::NODE );
        Vector< mtk::Vertex const * > tVertices( tNumNodes );

        for ( uint i = 0; i < tNumNodes; i++ )
        {
            tVertices( i ) = &mCutIgMesh->get_mtk_vertex( i );
        }
        return tVertices;
    }

    //------------------------------------------------------------------------------

    moris_id
    Enriched_Integration_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index       aEntityIndex,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        return mCutIgMesh->get_glb_entity_id_from_entity_loc_index( aEntityIndex, aEntityRank );
    }

    //------------------------------------------------------------------------------

    std::unordered_map< moris_id, moris_index >
    Enriched_Integration_Mesh::get_vertex_glb_id_to_loc_vertex_ind_map() const
    {
        return mCutIgMesh->get_vertex_glb_id_to_loc_vertex_ind_map();
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id          aEntityId,
            mtk::EntityRank   aEntityRank,
            const moris_index aIndex ) const
    {
        return mCutIgMesh->get_loc_entity_ind_from_entity_glb_id( aEntityId, aEntityRank );
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Integration_Mesh::get_entity_connected_to_entity_glob_ids(
            moris_id          aEntityId,
            mtk::EntityRank   aInputEntityRank,
            mtk::EntityRank   aOutputEntityRank,
            const moris_index aIndex ) const
    {
        moris_index tEntityIndex = get_loc_entity_ind_from_entity_glb_id( aEntityId, aInputEntityRank );

        Matrix< IndexMat > tEntityToEntityLoc = this->get_entity_connected_to_entity_loc_inds( tEntityIndex, aInputEntityRank, aOutputEntityRank );

        return convert_indices_to_ids( tEntityToEntityLoc, aOutputEntityRank );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    Enriched_Integration_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
    {
        mtk::Vertex const &tVertex = get_mtk_vertex( aNodeIndex );
        return tVertex.get_coords();
    }

    //------------------------------------------------------------------------------

    mtk::Vertex &
    Enriched_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex )
    {
        return mCutIgMesh->get_mtk_vertex( aVertexIndex );
    }

    //------------------------------------------------------------------------------

    mtk::Vertex const &
    Enriched_Integration_Mesh::get_mtk_vertex( moris_index aVertexIndex ) const
    {
        return mCutIgMesh->get_mtk_vertex( aVertexIndex );
    }

    //------------------------------------------------------------------------------

    mtk::Cell &
    Enriched_Integration_Mesh::get_writable_mtk_cell( moris_index aElementIndex )
    {
        return mCutIgMesh->get_mtk_cell( aElementIndex );
    }

    //------------------------------------------------------------------------------

    mtk::Cell &
    Enriched_Integration_Mesh::get_mtk_cell( moris_index aElementIndex )
    {
        return mCutIgMesh->get_mtk_cell( aElementIndex );
    }

    //------------------------------------------------------------------------------

    mtk::Cell const &
    Enriched_Integration_Mesh::get_mtk_cell( moris_index aElementIndex ) const
    {
        return mCutIgMesh->get_mtk_cell( aElementIndex );
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Integration_Mesh::get_communication_table() const
    {
        return mCutIgMesh->get_communication_table();
    }

    //------------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Integration_Mesh::get_set_names( mtk::EntityRank aSetEntityRank ) const
    {
        switch ( aSetEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                return mVertexSetNames;
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                MORIS_ASSERT( this->get_facet_rank() == mtk::EntityRank::EDGE,
                        "Enriched_Integration_Mesh::get_set_names() - side sets are defined on edges in 2d" );
                return mSideSetLabels;
                break;
            }
            case mtk::EntityRank::FACE:
            {
                MORIS_ASSERT( this->get_facet_rank() == mtk::EntityRank::FACE,
                        "Enriched_Integration_Mesh::get_set_names() - side sets are defined on faces in 3d" );
                return mSideSetLabels;
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                return mBlockSetNames;
                break;
            }
            default:
            {
                MORIS_ERROR( false,
                        "Enriched_Integration_Mesh::get_set_names() - "
                        "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
            }
                return Vector< std::string >( 0 );
                break;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_set_entity_loc_inds(
            mtk::EntityRank aSetEntityRank,
            std::string     aSetName ) const
    {
        switch ( aSetEntityRank )
        {
            case mtk::EntityRank::NODE:
            {
                // get the vertex set index
                auto tSetIndex = mVertexSetLabelToOrd.find( aSetName );

                Vector< moris::mtk::Vertex * > tVerticesInSet = mVerticesInVertexSet( tSetIndex->second );
                Matrix< IndexMat >             tVerticesInSetMat( 1, tVerticesInSet.size() );
                for ( uint i = 0; i < tVerticesInSet.size(); i++ )
                {
                    tVerticesInSetMat( i ) = tVerticesInSet( i )->get_index();
                }

                return tVerticesInSetMat;
                break;
            }
            case mtk::EntityRank::EDGE:
            {
                MORIS_ASSERT( this->get_facet_rank() == mtk::EntityRank::EDGE, "side sets are defined on edges in 2d" );
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
            case mtk::EntityRank::FACE:
            {
                MORIS_ASSERT( this->get_facet_rank() == mtk::EntityRank::FACE, "side sets are defined on faces in 3d" );
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
            case mtk::EntityRank::ELEMENT:
            {
                return this->get_block_entity_loc_inds( aSetName );
                break;
            }
            default:
            {
                MORIS_ERROR( 0, "Currently only supporting block, node and side sets in XTK enriched integration meshes" );
                return Matrix< IndexMat >( 0, 0 );
                break;
            }
        }
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_element_indices_in_block_set( uint aSetIndex )
    {
        std::string tSetName = this->get_set_names( mtk::EntityRank::ELEMENT )( aSetIndex );
        return this->get_block_entity_loc_inds( tSetName );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::get_sideset_elems_loc_inds_and_ords(
            const std::string  &aSetName,
            Matrix< IndexMat > &aElemIndices,
            Matrix< IndexMat > &aSidesetOrdinals ) const
    {
        // get the index
        moris_index tSideSetIndex = this->get_side_set_index( aSetName );

        // get the cell clusters
        Vector< mtk::Cluster const * > tSideClusters = this->get_side_set_cluster( tSideSetIndex );

        // iterate through side clusters and count number of sides in set
        uint tNumSides = 0;
        for ( auto iCluster : tSideClusters )
        {
            tNumSides = tNumSides + iCluster->get_num_primary_cells();
        }

        // size outputs
        aElemIndices.resize( 1, tNumSides );
        aSidesetOrdinals.resize( 1, tNumSides );

        // reset count
        tNumSides = 0;
        for ( auto iCluster : tSideClusters )
        {
            Matrix< IndexMat > tSideOrdinals = iCluster->get_cell_side_ordinals();
            Matrix< IndexMat > tCellIndices  = iCluster->get_primary_cell_indices_in_cluster();

            aElemIndices( { 0, 0 }, { tNumSides, tNumSides + tCellIndices.numel() - 1 } )     = tCellIndices.matrix_data();
            aSidesetOrdinals( { 0, 0 }, { tNumSides, tNumSides + tCellIndices.numel() - 1 } ) = tSideOrdinals.matrix_data();

            tNumSides = tNumSides + tSideOrdinals.numel();
        }
    }

    //------------------------------------------------------------------------------

    moris_id
    Enriched_Integration_Mesh::get_max_entity_id( mtk::EntityRank aEntityRank, const moris_index aIndex ) const
    {
        MORIS_ASSERT( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT, "Only Elements or Nodes have max id" );

        uint tNumEntities = this->get_num_entities( aEntityRank );

        moris_id tLocalMaxId = 0;

        for ( uint i = 0; i < tNumEntities; i++ )
        {
            moris_id tId = this->get_glb_entity_id_from_entity_loc_index( i, aEntityRank );

            if ( tId > tLocalMaxId )
            {
                tLocalMaxId = tId;
            }
        }

        moris_id tGlobalMaxId = moris::max_all( tLocalMaxId );
        return tGlobalMaxId;
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_node_owner( moris_index aNodeIndex ) const
    {
        return this->get_mtk_vertex( aNodeIndex ).get_owner();
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_element_owner( moris_index aElementIndex ) const
    {
        return this->get_mtk_cell( aElementIndex ).get_owner();
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_block_entity_loc_inds( std::string aSetName ) const
    {
        // get index
        moris_index tBlockIndex = this->get_block_set_index( aSetName );

        // get clusters in block
        Vector< mtk::Cluster const * > tCellClusters = this->get_cell_clusters_in_set( tBlockIndex );

        // iterate through clusters and count all primary integration cells
        uint tCount = 0;
        for ( auto it : tCellClusters )
        {
            tCount = tCount + it->get_num_primary_cells();
        }

        // allocate output
        Matrix< IndexMat > tCellIndices( 1, tCount );

        // reset count to use it for something else
        tCount = 0;

        // iterate through clusters and collect all primary integration cell indices
        for ( auto it : tCellClusters )
        {
            // get primary cell clusters
            Matrix< IndexMat > tPrimaryCellIndices = it->get_primary_cell_indices_in_cluster();

            tCellIndices( { 0, 0 }, { tCount, tCount + tPrimaryCellIndices.numel() - 1 } ) = tPrimaryCellIndices.matrix_data();

            tCount = tCount + tPrimaryCellIndices.numel();
        }

        return tCellIndices;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_dbl_sided_interface_sets(
            Vector< moris_index > aLeaderBulkPhaseIndices,
            Vector< moris_index > aFollowerBulkPhaseIndices )
    {
        // get the number sets to be created
        uint tNumDblSideSets = aLeaderBulkPhaseIndices.size();

        // check that the input is valid
        MORIS_ERROR( tNumDblSideSets == aFollowerBulkPhaseIndices.size(),
                "Enriched_Integration_Mesh::create_dbl_sided_interface_set() - "
                "List of Leader and Follower bulk phases are of different length." );

        // initialize list of double sided side set names
        Vector< std::string > tDblSideSetNames( tNumDblSideSets );

        // get the names for all dbl side sets to be created
        for ( uint iDblSideSet = 0; iDblSideSet < tNumDblSideSets; iDblSideSet++ )
        {
            // get the leader and follower bulk phase indices
            moris_index tLeaderBulkPhaseIndex   = aLeaderBulkPhaseIndices( iDblSideSet );
            moris_index tFollowerBulkPhaseIndex = aFollowerBulkPhaseIndices( iDblSideSet );

            // get the name of this side set's name
            tDblSideSetNames( iDblSideSet ) = this->get_dbl_interface_side_set_name( tLeaderBulkPhaseIndex, tFollowerBulkPhaseIndex );
        }

        // register the double sided side set names and get the sets' ordinals
        Vector< moris_index > tDblSideSetOrds = this->register_double_side_set_names( tDblSideSetNames );

        //
        for ( uint iDblSideSet = 0; iDblSideSet < tNumDblSideSets; iDblSideSet++ )
        {
            // get the leader and follower bulk phase indices
            moris_index tLeaderBulkPhaseIndex   = aLeaderBulkPhaseIndices( iDblSideSet );
            moris_index tFollowerBulkPhaseIndex = aFollowerBulkPhaseIndices( iDblSideSet );

            // get the set ordinal
            moris_index tDblSideSetOrd = tDblSideSetOrds( iDblSideSet );

            // get the side set colors
            Matrix< IndexMat > tLeaderColor   = { { tLeaderBulkPhaseIndex } };
            Matrix< IndexMat > tFollowerColor = { { tFollowerBulkPhaseIndex } };
            this->set_double_side_set_colors( tDblSideSetOrd, tLeaderColor, tFollowerColor );

            // get the other interface index
            moris_index tOtherInterfaceIndex = this->get_dbl_side_set_index( tFollowerBulkPhaseIndex, tLeaderBulkPhaseIndex );

            // appropriately resize member data
            uint tNumClusters = mDoubleSideSetsLeaderIndex( tOtherInterfaceIndex ).size();
            mDoubleSideSetsLeaderIndex( tDblSideSetOrd ).resize( tNumClusters );
            mDoubleSideSetsFollowerIndex( tDblSideSetOrd ).resize( tNumClusters );

            // create double sided side clusters
            for ( uint iCluster = 0; iCluster < tNumClusters; iCluster++ )
            {
                // leader is follower, follower is leader
                mDoubleSideSetsLeaderIndex( tDblSideSetOrd )( iCluster )   = mDoubleSideSetsFollowerIndex( tOtherInterfaceIndex )( iCluster );
                mDoubleSideSetsFollowerIndex( tDblSideSetOrd )( iCluster ) = mDoubleSideSetsLeaderIndex( tOtherInterfaceIndex )( iCluster );

                // get leader and follower clusters
                Side_Cluster *tLeaderSideCluster   = mDoubleSideSingleSideClusters( mDoubleSideSetsLeaderIndex( tDblSideSetOrd )( iCluster ) ).get();
                Side_Cluster *tFollowerSideCluster = mDoubleSideSingleSideClusters( mDoubleSideSetsFollowerIndex( tDblSideSetOrd )( iCluster ) ).get();

                // create double side cluster
                std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster =
                        std::make_shared< mtk::Double_Side_Cluster >(
                                tLeaderSideCluster,
                                tFollowerSideCluster,
                                tLeaderSideCluster->get_vertices_in_cluster() );

                // store dbl sided cluster on list of all clusters (globally and per set)
                mDoubleSideClusters.push_back( tDblSideCluster );
                mDoubleSideSets( tDblSideSetOrd ).push_back( tDblSideCluster );
            }

        }    // end for: loop over double sided side sets to be created

        // finalize the set coloring
        this->setup_color_to_set();

        // commit the created double sided side sets to the mesh
        for ( uint iDblSideSet = 0; iDblSideSet < tNumDblSideSets; iDblSideSet++ )
        {
            moris_index tDblSideSetOrd = tDblSideSetOrds( iDblSideSet );
            this->commit_double_side_set( tDblSideSetOrd );
        }

        // finalize the list of double sided side sets
        this->collect_all_sets();

        // communicate the double sided side sets after new ones have been added
        this->communicate_sets_of_type( mtk::SetType::DOUBLE_SIDED_SIDESET );

    }    // end function: Enriched_Integration_Mesh::create_dbl_sided_interface_sets()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::deactivate_empty_sets()
    {
        this->deactivate_empty_block_sets();
        this->deactivate_empty_side_sets();
        this->setup_color_to_set();
        this->collect_all_sets();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::deactivate_empty_side_sets()
    {
        // copy old data
        std::unordered_map< std::string, moris_index >           tOldSetMap      = mSideSideSetLabelToOrd;
        Vector< std::string >                                    tOldSetNames    = mSideSetLabels;
        Vector< Vector< std::shared_ptr< xtk::Side_Cluster > > > tOldSetClusters = mSideSets;

        // clear member data
        mSideSideSetLabelToOrd.clear();
        mSideSetLabels.clear();
        mSideSets.clear();

        for ( auto iB : mListOfSideSets )
        {
            delete iB;
        }
        mListOfSideSets.clear();

        // current index
        moris_index tSetIndex = 0;

        // build list of side set indices
        Vector< moris_index > tSideSetIndexList;
        tSideSetIndexList.reserve( tOldSetClusters.size() );

        for ( uint i = 0; i < tOldSetClusters.size(); i++ )
        {
            uint tMySize  = tOldSetClusters( i ).size();
            uint tAllSize = sum_all( tMySize );

            if ( tAllSize > 0 )
            {
                mSideSetLabels.push_back( tOldSetNames( i ) );
                mSideSets.push_back( tOldSetClusters( i ) );

                MORIS_ASSERT( mSideSideSetLabelToOrd.find( tOldSetNames( i ) ) == mSideSideSetLabelToOrd.end(), "Duplicate block set in mesh" );
                mSideSideSetLabelToOrd[ tOldSetNames( i ) ] = tSetIndex;
                tSideSetIndexList.push_back( tSetIndex );
                tSetIndex++;
            }
        }

        // commit and communicate the side sets
        this->commit_side_set( tSideSetIndexList );
        this->communicate_sets_of_type( mtk::SetType::SIDESET );

        this->setup_color_to_set();
        this->collect_all_sets();
    }

    void
    Enriched_Integration_Mesh::add_mesh_field_real_scalar_data_loc_inds(
            std::string const      &aFieldName,
            mtk::EntityRank const  &aFieldEntityRank,
            Matrix< DDRMat > const &aFieldData )
    {

        MORIS_ASSERT( aFieldEntityRank == mtk::EntityRank::ELEMENT, "Only tested for nodal and element scalar field" );

        moris_index tFieldIndex = this->create_field( aFieldName, mtk::EntityRank::ELEMENT, 0 );

        this->add_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, aFieldData );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_cell_id_fields()
    {
        // Fields constructed here
        Vector< std::string > tCellFields = { "cell_id" };

        moris_index tFieldIndex = this->create_field( tCellFields( 0 ), mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );

        Matrix< DDRMat > tCellIdField( 1, this->get_num_elems() );

        for ( uint iCell = 0; iCell < this->get_num_elems(); iCell++ )
        {
            tCellIdField( iCell ) = (moris::real)this->get_mtk_cell( iCell ).get_id();
        }

        this->add_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, tCellIdField );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::deactivate_empty_block_sets()
    {
        std::unordered_map< std::string, moris_index > tOldSetMap      = mBlockSetLabelToOrd;
        Vector< std::string >                          tOldSetNames    = mBlockSetNames;
        Vector< mtk::CellTopology >                    tOldSetTopo     = mBlockSetTopology;
        Vector< Vector< xtk::Cell_Cluster const * > >  tOldSetClusters = mPrimaryBlockSetClusters;
        Vector< Matrix< IndexMat > >                   tOldColors      = mBlockSetColors;

        // clear member data
        mBlockSetLabelToOrd.clear();
        mBlockSetNames.clear();
        mBlockSetTopology.clear();
        mPrimaryBlockSetClusters.clear();
        mBlockSetColors.clear();

        for ( auto iB : mListOfBlocks )
        {
            delete iB;
        }
        mListOfBlocks.clear();

        // current index
        moris_index tSetIndex = 0;

        for ( uint i = 0; i < tOldSetClusters.size(); i++ )
        {
            uint tMySize  = tOldSetClusters( i ).size();
            uint tAllSize = sum_all( tMySize );

            if ( tAllSize > 0 )
            {
                mBlockSetNames.push_back( tOldSetNames( i ) );
                mPrimaryBlockSetClusters.push_back( tOldSetClusters( i ) );
                mBlockSetTopology.push_back( tOldSetTopo( i ) );
                mBlockSetColors.push_back( tOldColors( i ) );

                MORIS_ASSERT( mBlockSetLabelToOrd.find( tOldSetNames( i ) ) == mBlockSetLabelToOrd.end(), "Duplicate block set in mesh" );
                mBlockSetLabelToOrd[ tOldSetNames( i ) ] = tSetIndex;
                this->commit_block_set( tSetIndex );
                tSetIndex++;
            }
        }

        this->setup_color_to_set();
        this->collect_all_sets();

        // communicate committed block sets
        this->communicate_sets_of_type( mtk::SetType::BULK );
    }

    //------------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Integration_Mesh::create_basis_support_fields( Matrix< DDRMat > const &aProbeSpheres )
    {

        Vector< std::string > tFieldNames;

// DEBUG because basis coordinates is only defined on debug mode
#ifdef MORIS_HAVE_DEBUG
        MORIS_ASSERT( aProbeSpheres.n_cols() == 4, "Probe sphere should be r, xc, yc, zc" );
        moris_index tNumSpheres = aProbeSpheres.n_rows();

        // background mesh data
        moris::mtk::Interpolation_Mesh &tMeshData = mModel->get_background_mesh();

        // get the enriched interpolation mesh
        Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( 0 );

        // base string of field
        std::string tBaseStr = "weights";

        // determine which basis functions we are visualizing
        Vector< Vector< moris_index > >                          tActiveBasis( tEnrInterpMesh->get_num_interpolation_types() );
        Vector< std::unordered_map< moris_index, moris_index > > tEnrCoeffActiveIndexFieldIndex( tEnrInterpMesh->get_num_interpolation_types() );

        moris_index tFieldIndex = 0;

        for ( uint iBT = 0; iBT < tEnrInterpMesh->mMeshIndices.numel(); iBT++ )
        {
            moris_index tMeshIndex = iBT;

            // iterate through background basis functions
            for ( uint iBackBasisIndex = 0; iBackBasisIndex < tEnrInterpMesh->get_num_background_coefficients( tMeshIndex ); iBackBasisIndex++ )
            {
                // get the basis coordinate of the background basis function
                Matrix< DDRMat > tBasisCoords = tMeshData.get_basis_coords( tMeshIndex, (moris_index)iBackBasisIndex );

                // iterate through circles, see if the basis is active
                for ( moris_index iSp = 0; iSp < tNumSpheres; iSp++ )
                {
                    // initialize  the level set value
                    moris::real tLSVal = 0;

                    // determine the value of level set at the basis coordinates based on the dimension
                    if ( this->get_spatial_dim() == 3 )
                    {
                        tLSVal = sqrt( pow( tBasisCoords( 0 ) - aProbeSpheres( iSp, 1 ), 2 ) + pow( tBasisCoords( 1 ) - aProbeSpheres( iSp, 2 ), 2 ) + pow( tBasisCoords( 2 ) - aProbeSpheres( iSp, 3 ), 2 ) ) - aProbeSpheres( iSp, 0 );
                    }
                    else
                    {
                        tLSVal = sqrt( pow( tBasisCoords( 0 ) - aProbeSpheres( iSp, 1 ), 2 ) + pow( tBasisCoords( 1 ) - aProbeSpheres( iSp, 2 ), 2 ) ) - aProbeSpheres( iSp, 0 );
                    }

                    if ( tLSVal < 0.0 )
                    {
                        // iterate through enriched interpolation coeffs
                        Matrix< IndexMat > const &tEnrCoeffs = tEnrInterpMesh->get_enriched_coefficients_at_background_coefficient( tMeshIndex, (moris_index)iBackBasisIndex );

                        for ( uint iEnrBasisOrd = 0; iEnrBasisOrd < tEnrCoeffs.numel(); iEnrBasisOrd++ )
                        {
                            const moris_index tEnrIndex = tEnrCoeffs( iEnrBasisOrd );

                            tActiveBasis( tMeshIndex ).push_back( tEnrIndex );
                            tEnrCoeffActiveIndexFieldIndex( tMeshIndex )[ tEnrIndex ] = tFieldIndex;
                            tFieldIndex++;
                        }
                    }
                }
            }
        }

        // field names for output
        Vector< std::string > tOutputFieldNames;

        // field information for internal use
        tFieldNames.resize( tFieldIndex );
        Vector< moris_index >      tFieldIndices( tFieldIndex );
        Vector< Matrix< DDRMat > > tFieldData( tFieldIndex, Matrix< DDRMat >( 1, this->get_num_nodes(), -10 ) );

        // iterate through interpolation types and for each basis declare the field in mesh
        for ( uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
        {

            moris_index tMeshIndex     = iBT;
            std::string tInterpTypeStr = "_mi_" + std::to_string( tMeshIndex );

            // iterate through basis functions
            for ( uint iB = 0; iB < tActiveBasis( tMeshIndex ).size(); iB++ )
            {
                MORIS_ASSERT( tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).find( tActiveBasis( tMeshIndex )( iB ) ) != tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).end(), "Not in map" );
                tFieldIndex = tEnrCoeffActiveIndexFieldIndex( tMeshIndex ).find( tActiveBasis( tMeshIndex )( iB ) )->second;

                tFieldNames( tFieldIndex ) = tBaseStr + tInterpTypeStr + "_ind_" + std::to_string( tActiveBasis( tMeshIndex )( iB ) );

                // declare the field in this mesh
                tFieldIndices( tFieldIndex ) = this->create_field( tFieldNames( tFieldIndex ), mtk::EntityRank::NODE, MORIS_INDEX_MAX );
            }
        }

        // populate field data
        for ( uint iCl = 0; iCl < this->mCellClusters.size(); iCl++ )
        {
            xtk::Cell_Cluster *tCluster = mCellClusters( iCl ).get();

            // get the interpolation cell
            moris::mtk::Cell const &tInterpCell = tCluster->get_interpolation_cell();

            // get the vertices attached to this cell
            Vector< moris::mtk::Vertex * > tVertices = tInterpCell.get_vertex_pointers();

            // allocate place to put coefficients interpolating into these vertices
            Vector< Vector< moris_index > > tCoeffsIPIntoCluster( tEnrInterpMesh->get_num_interpolation_types() );

            // collect coefficients of this interpolation cell
            for ( uint iV = 0; iV < tVertices.size(); iV++ )
            {
                for ( uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
                {
                    mtk::Vertex_Interpolation *tVertexIp = tVertices( iV )->get_interpolation( tEnrInterpMesh->get_interpolation_index( (moris_index)iBT ) );

                    // get indices of coefficients
                    Matrix< IndexMat > tCoeffInds = tVertexIp->get_indices();

                    for ( uint iC = 0; iC < tCoeffInds.numel(); iC++ )
                    {
                        tCoeffsIPIntoCluster( iBT ).push_back( tCoeffInds( iC ) );
                    }
                }
            }

            // get primary cells from the cluster
            Vector< moris::mtk::Cell const * > const &tPrimaryCells = tCluster->get_primary_cells_in_cluster();

            // iterate through primary cells
            for ( uint iCell = 0; iCell < tPrimaryCells.size(); iCell++ )
            {
                // get vertices attached to primary cells
                Vector< moris::mtk::Vertex * > tVertices = tPrimaryCells( iCell )->get_vertex_pointers();

                // iterate through vertices and mark them as in support of all coefficients in tCoeffsIPIntoCluster
                for ( uint iV = 0; iV < tVertices.size(); iV++ )
                {
                    for ( uint iBT = 0; iBT < tEnrInterpMesh->get_num_interpolation_types(); iBT++ )
                    {
                        for ( uint iC = 0; iC < tCoeffsIPIntoCluster( iBT ).size(); iC++ )
                        {
                            auto tFieldIndIter = tEnrCoeffActiveIndexFieldIndex( iBT ).find( tCoeffsIPIntoCluster( iBT )( iC ) );
                            if ( tFieldIndIter != tEnrCoeffActiveIndexFieldIndex( iBT ).end() )
                            {
                                tFieldData( tFieldIndIter->second )( tVertices( iV )->get_index() ) = 1;
                            }
                        }
                    }
                }
            }
        }

        // add field data to mesh
        // iterate through interpolation
        for ( uint iField = 0; iField < tFieldIndices.size(); iField++ )
        {
            this->add_field_data( tFieldIndices( iField ), mtk::EntityRank::NODE, tFieldData( iField ), MORIS_INDEX_MAX );
        }

#endif

        return tFieldNames;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::write_mesh( moris::ParameterList *aParamList )
    {
        if ( aParamList->get< bool >( "deactivate_empty_sets" ) )
        {
            this->deactivate_empty_sets();
        }

        // get path to output XTK files to
        std::string tOutputPath = aParamList->get< std::string >( "output_path" );
        std::string tOutputFile = aParamList->get< std::string >( "output_file" );
        std::string tOutputBase = tOutputFile.substr( 0, tOutputFile.find( "." ) );
        std::string tOutputExt  = tOutputFile.substr( tOutputFile.find( "." ), tOutputFile.length() );

        // make sure the output mesh file name has correct file ending
        MORIS_ASSERT( tOutputExt == ".exo" || tOutputExt == ".e", "Invalid file extension, needs to be .exo or .e" );

        // Write mesh
        moris::mtk::Writer_Exodus tExodusWriter( this );

        // if user requests to keep XTK output for all iterations, add iteration count to output file name
        if ( aParamList->get< bool >( "keep_all_opt_iters" ) )
        {
            // get optimization iteration ( function returns zero if no optimization )
            uint tOptIter = gLogger.get_opt_iteration();

            tExodusWriter.write_mesh(
                    "",
                    tOutputPath + tOutputBase + "." + std::to_string( tOptIter ) + tOutputExt,
                    "",
                    tOutputPath + "xtk_temp." + std::to_string( tOptIter ) + tOutputExt );
        }
        // otherwise, proceed as usual and overwrite xtk_temp.exo each iteration
        else
        {
            tExodusWriter.write_mesh( "", tOutputPath + tOutputFile, "", tOutputPath + "xtk_temp.exo" );
        }

        //----------------------------------------------------------------
        // create data for XTK data if requested

        if ( aParamList->get< bool >( "write_enrichment_fields" ) )
        {
            std::string           tProbeSpheresStr = aParamList->get< std::string >( "write_enrichment_fields_probe_spheres" );
            Vector< std::string > tNodeFields;

            if ( !tProbeSpheresStr.empty() )
            {
                Matrix< DDRMat > tProbeSpheres = string_to_mat< DDRMat >( tProbeSpheresStr );

                // set up the nodal fields for basis support
                this->create_basis_support_fields( tProbeSpheres );
            }

            // Vector<std::string> tEnrichmentFieldNames =  mModel->get_basis_enrichment().get_cell_enrichment_field_names();
            // tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames, this);

            // place an element field in the mesh
            this->create_bg_cell_id_field();

            // subphase neighbor field
            this->create_subphase_fields();
        }

        //----------------------------------------------------------------
        // write nodal fields

        Vector< std::string > tNodeFields = this->get_field_names( mtk::EntityRank::NODE, MORIS_INDEX_MAX );
        tExodusWriter.set_nodal_fields( tNodeFields );

        for ( uint iF = 0; iF < tNodeFields.size(); iF++ )
        {
            moris::moris_index tFieldIndex = this->get_field_index( tNodeFields( iF ), mtk::EntityRank::NODE );
            tExodusWriter.write_nodal_field( tNodeFields( iF ), this->get_field_data( tFieldIndex, mtk::EntityRank::NODE, MORIS_INDEX_MAX ) );
        }

        //----------------------------------------------------------------
        // write elemental fields

        // create elemental field containing the cell IDs
        this->create_cell_id_fields();

        // get a list of names of all the fields to be written
        Vector< std::string > tCellFields = this->get_field_names( mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );

        // set the field names in the exodus file
        tExodusWriter.set_elemental_fields( tCellFields );

        // get a list of the labels of all the blocks to be written to exodus
        Vector< std::string > tBlockNames = this->get_set_names( mtk::EntityRank::ELEMENT );

        // write every field
        for ( uint iField = 0; iField < tCellFields.size(); iField++ )
        {
            // get the index of the of the current field in the stored list of fields
            moris::moris_index      tFieldIndex = this->get_field_index( tCellFields( iField ), mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );
            Matrix< DDRMat > const &tFieldData  = this->get_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );

            // write field on every block
            for ( uint iBlock = 0; iBlock < this->get_num_blocks(); iBlock++ )
            {
                // get the label and index of the current block in the internal mesh and all its IG elements
                std::string        tBlockName   = tBlockNames( iBlock );
                moris_index        tBlockIndex  = this->get_block_set_index( tBlockName );
                Matrix< IndexMat > tCellIndices = this->get_element_indices_in_block_set( tBlockIndex );

                // initialize and collect field data local to the current block
                Matrix< DDRMat > tBlockFieldData( 1, tCellIndices.numel(), -10.0 );
                for ( uint iCell = 0; iCell < tCellIndices.numel(); iCell++ )
                {
                    tBlockFieldData( iCell ) = tFieldData( tCellIndices( iCell ) );
                }

                // write data to mesh (unless it is empty)
                if ( tBlockFieldData.numel() > 0 )
                {
                    tExodusWriter.write_elemental_field( tBlockName, tCellFields( iField ), tBlockFieldData );
                }
            }
        }

        //----------------------------------------------------------------
        // write side set fields

        // collect names for all side set fields
        Vector< std::string >                tAllSideSetFields;
        std::map< std::string, moris_index > tAllSideSetFieldsMap;
        for ( uint iSideSet = 0; iSideSet < this->get_num_side_sets(); iSideSet++ )
        {
            // get a list of names of all the fields stored for the current side set
            Vector< std::string > tFieldLabels = this->get_field_names( this->get_facet_rank(), iSideSet );

            // go through fields and add them to the map if not in there already
            for ( std::string iFieldName : tFieldLabels )
            {
                if ( tAllSideSetFieldsMap.find( iFieldName ) == tAllSideSetFieldsMap.end() )
                {
                    tAllSideSetFieldsMap[ iFieldName ] = tAllSideSetFields.size();
                    tAllSideSetFields.push_back( iFieldName );
                }
            }
        }

        // add field labels to exodus mesh
        tExodusWriter.set_side_set_fields( tAllSideSetFields );

        // get a list of the labels of all the blocks to be written to exodus
        Vector< std::string > tSideSetNames = this->get_set_names( this->get_facet_rank() );

        // write field for every side set
        for ( uint iSideSet = 0; iSideSet < this->get_num_side_sets(); iSideSet++ )
        {
            // get a list of names of all the fields to be written to the current side set
            Vector< std::string > tFieldLabels = this->get_field_names( this->get_facet_rank(), iSideSet );

            // get the label of the current set
            std::string tSideSetName = tSideSetNames( iSideSet );

            // write every field
            for ( uint iField = 0; iField < tFieldLabels.size(); iField++ )
            {
                // get the name of the current field
                std::string tFieldName = tFieldLabels( iField );

                // get the index of the of the current field in the stored list of fields
                moris::moris_index      tFieldIndex = this->get_field_index( tFieldName, this->get_facet_rank(), iSideSet );
                Matrix< DDRMat > const &tFieldData  = this->get_field_data( tFieldIndex, this->get_facet_rank(), iSideSet );

                // write data to mesh (unless it is empty)
                if ( tFieldData.numel() > 0 )
                {
                    // make sure the data written is of correct size
                    MORIS_ASSERT( this->get_num_elements_in_side_set( tSideSetName ) == tFieldData.numel(),
                            "XTK::Enriched_Integration_Mesh::write_mesh() - "
                            "Trying to write side set field data of incorrect size to xtk mesh output." );

                    // output to mesh writer
                    tExodusWriter.write_side_set_field( tSideSetName, tFieldName, tFieldData );
                }
            }

        }    // end for: write field for every side set

        //----------------------------------------------------------------
        // finalize exodus writer

        tExodusWriter.set_time( 0.0 );
        tExodusWriter.close_file();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_bg_cell_id_field()
    {
        // Fields constructed here
        Vector< std::string > tCellFields = { "bg_cell_id" };

        // set field index
        moris_index tFieldIndex = this->create_field( tCellFields( 0 ), mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );

        Matrix< DDRMat > tCellIdField( 1, this->get_num_elems() );

        //
        for ( uint iEnrIpCell = 0; iEnrIpCell < mSubphaseIndexToClusterIndex.numel(); iEnrIpCell++ )
        {
            // get the subphase index
            moris_index tSubphaseClusterIndex = mSubphaseIndexToClusterIndex( iEnrIpCell );

            // Cell Cluster
            moris_index                               tBaseIpCellId     = mCellClusters( tSubphaseClusterIndex )->get_xtk_interpolation_cell()->get_base_cell()->get_id();
            Vector< moris::mtk::Cell const * > const &tIgCellsInCluster = mCellClusters( tSubphaseClusterIndex )->get_primary_cells_in_cluster();

            for ( uint iCell = 0; iCell < tIgCellsInCluster.size(); iCell++ )
            {
                tCellIdField( tIgCellsInCluster( iCell )->get_index() ) = std::floor( (moris::real)tBaseIpCellId );
            }
        }

        this->add_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, tCellIdField );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_subphase_fields()
    {
        // this->compute_subphase_centroids( "sp_centroid.csv" );

        // this->write_bg_cell_to_subphase_index( "bg_cell_to_sp_index.csv" );

        this->write_subphase_neighborhood( "subphase_neighborhood.csv" );

        // Fields constructed here
        Vector< std::string > tCellFields = { "sp_index", "bulk_phase" };

        Vector< moris::moris_index > tFieldIndices( tCellFields.size() );

        Matrix< DDRMat > tCellToSubphase( 1, this->get_num_elems() );

        Matrix< DDRMat > tCellToBulkPhase( 1, this->get_num_elems() );

        for ( uint i = 0; i < mCellClusters.size(); i++ )
        {
            // Cell Cluster
            std::shared_ptr< xtk::Cell_Cluster > tCluster = mCellClusters( i );

            moris_index tSubphaseIndex  = tCluster->get_xtk_interpolation_cell()->get_subphase_index();
            moris_index tBulkPhaseIndex = tCluster->get_xtk_interpolation_cell()->get_bulkphase_index();

            // get the cells in cluster
            Vector< moris::mtk::Cell const * > const &tIgCellsInCluster = tCluster->get_primary_cells_in_cluster();

            for ( uint iCell = 0; iCell < tIgCellsInCluster.size(); iCell++ )
            {
                tCellToSubphase( tIgCellsInCluster( iCell )->get_index() )  = std::floor( (moris::real)tSubphaseIndex );
                tCellToBulkPhase( tIgCellsInCluster( iCell )->get_index() ) = std::floor( (moris::real)tBulkPhaseIndex );
            }
        }

        for ( uint iF = 0; iF < tCellFields.size(); iF++ )
        {
            tFieldIndices( iF ) = this->create_field( tCellFields( iF ), mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );
        }

        this->add_field_data( tFieldIndices( 0 ), mtk::EntityRank::ELEMENT, tCellToSubphase, MORIS_INDEX_MAX );
        this->add_field_data( tFieldIndices( 1 ), mtk::EntityRank::ELEMENT, tCellToBulkPhase, MORIS_INDEX_MAX );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::write_subphase_neighborhood( std::string aFile )
    {
        Vector< std::shared_ptr< Vector< moris_index > > > const &tSubphaseToSubphase = mCutIgMesh->get_subphase_neighborhood()->mSubphaseToSubPhase;

        std::ostringstream tStringStream;
        for ( uint iC = 0; iC < tSubphaseToSubphase.size(); iC++ )
        {

            for ( uint iN = 0; iN < tSubphaseToSubphase( iC )->size(); iN++ )
            {
                tStringStream << ( *tSubphaseToSubphase( iC ) )( iN );
                if ( iN != tSubphaseToSubphase( iC )->size() - 1 )
                {
                    tStringStream << ",";
                }
            }

            if ( tSubphaseToSubphase( iC )->size() == 0 )
            {
                tStringStream << "NaN";
            }
            tStringStream << std::endl;
        }

        if ( aFile.empty() == false )
        {
            std::ofstream tOutputFile( aFile );
            tOutputFile << tStringStream.str() << std::endl;
            tOutputFile.close();
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_union_block( Vector< std::string > const &aBlocks,
            std::string                                                         aNewBlock,
            Matrix< IndexMat > const                                           &aNewBlockColor )
    {
        MORIS_ERROR( aBlocks.size() >= 2, "Union needs to happen between two blocks or more" );

        mtk::CellTopology tCellTopo = mtk::CellTopology::UNDEFINED;

        uint                  tCount = 0;
        Vector< moris_index > tBlockIndices( aBlocks.size() );
        for ( uint i = 0; i < aBlocks.size(); i++ )
        {
            tBlockIndices( i ) = this->get_block_set_index( aBlocks( i ) );
            tCount             = tCount + mPrimaryBlockSetClusters( tBlockIndices( i ) ).size();

            moris::mtk::Set *tSet = this->get_set_by_index( this->get_set_index_by_name( aBlocks( i ) ) );
            if ( i == 0 )
            {
                tCellTopo = tSet->get_cell_topology();
            }
            else
            {
                MORIS_ERROR( tCellTopo == tSet->get_cell_topology(), "Invalid merge detected, verify that all blocks have the same cell topology" );
            }
        }

        Vector< moris_index > tBlockSetIndex = this->register_block_set_names_with_cell_topo( { aNewBlock }, tCellTopo );
        mPrimaryBlockSetClusters( tBlockSetIndex( 0 ) ).reserve( tCount );

        for ( uint i = 0; i < aBlocks.size(); i++ )
        {
            mPrimaryBlockSetClusters( tBlockSetIndex( 0 ) ).append( mPrimaryBlockSetClusters( tBlockIndices( i ) ) );
        }

        // Check compatibility of union

        this->commit_block_set( tBlockSetIndex( 0 ) );
        this->set_block_set_colors( tBlockSetIndex( 0 ), aNewBlockColor );
        this->setup_color_to_set();
        this->collect_all_sets( false );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_union_side_set(
            Vector< std::string > const &aSideSets,
            std::string                  aNewSideSet,
            Matrix< IndexMat > const    &aNewSideSetColor )
    {
        MORIS_ERROR( aSideSets.size() >= 2, "Union needs to happen between two side sets or more" );

        uint tCount = 0;

        Vector< moris_index > tSideSetIndices( aSideSets.size() );

        for ( uint i = 0; i < aSideSets.size(); i++ )
        {
            tSideSetIndices( i ) = this->get_side_set_index( aSideSets( i ) );
            tCount               = tCount + mSideSets( tSideSetIndices( i ) ).size();
        }

        Vector< moris_index > tSideSetIndex = this->register_side_set_names( { aNewSideSet } );
        mSideSets( tSideSetIndex( 0 ) ).reserve( tCount );

        for ( uint i = 0; i < aSideSets.size(); i++ )
        {
            mSideSets( tSideSetIndex( 0 ) ).append( mSideSets( tSideSetIndices( i ) ) );
        }

        this->commit_side_set( tSideSetIndex( 0 ) );

        this->set_side_set_colors( tSideSetIndex( 0 ), aNewSideSetColor );

        this->setup_color_to_set();
        this->collect_all_sets( false );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::deactivate_all_blocks_except_selected( Vector< std::string > const &aBlockSetsToKeep )
    {
        std::unordered_map< std::string, moris_index > tBlocksToKeepMap;
        for ( uint i = 0; i < aBlockSetsToKeep.size(); i++ )
        {
            tBlocksToKeepMap[ aBlockSetsToKeep( i ) ] = 1;
        }

        std::unordered_map< std::string, moris_index > tOldSetMap      = mBlockSetLabelToOrd;
        Vector< std::string >                          tOldSetNames    = mBlockSetNames;
        Vector< mtk::CellTopology >                    tOldSetTopo     = mBlockSetTopology;
        Vector< Vector< xtk::Cell_Cluster const * > >  tOldSetClusters = mPrimaryBlockSetClusters;
        Vector< Matrix< IndexMat > >                   tOldColors      = mBlockSetColors;

        // clear member data
        mBlockSetLabelToOrd.clear();
        mBlockSetNames.clear();
        mBlockSetTopology.clear();
        mPrimaryBlockSetClusters.clear();
        mBlockSetColors.clear();

        for ( auto iB : mListOfBlocks )
        {
            delete iB;
        }
        mListOfBlocks.clear();

        // current index
        moris_index tSetIndex = 0;

        for ( uint i = 0; i < tOldSetClusters.size(); i++ )
        {
            if ( tBlocksToKeepMap.find( tOldSetNames( i ) ) != tBlocksToKeepMap.end() )
            {
                mBlockSetNames.push_back( tOldSetNames( i ) );
                mPrimaryBlockSetClusters.push_back( tOldSetClusters( i ) );
                mBlockSetTopology.push_back( tOldSetTopo( i ) );
                mBlockSetColors.push_back( tOldColors( i ) );

                MORIS_ASSERT( mBlockSetLabelToOrd.find( tOldSetNames( i ) ) == mBlockSetLabelToOrd.end(), "Duplicate block set in mesh" );
                mBlockSetLabelToOrd[ tOldSetNames( i ) ] = tSetIndex;
                this->commit_block_set( tSetIndex );
                tSetIndex++;
            }
        }

        // communicate committed block sets
        this->communicate_sets_of_type( mtk::SetType::BULK );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::deactivate_all_side_sets_except_selected( Vector< std::string > const &aSideSetsToKeep )
    {
        std::unordered_map< std::string, moris_index > tSideSetsToKeepMap;

        for ( uint i = 0; i < aSideSetsToKeep.size(); i++ )
        {
            tSideSetsToKeepMap[ aSideSetsToKeep( i ) ] = 1;
        }

        // copy old data
        std::unordered_map< std::string, moris_index >           tOldSetMap      = mSideSideSetLabelToOrd;
        Vector< std::string >                                    tOldSetNames    = mSideSetLabels;
        Vector< Vector< std::shared_ptr< xtk::Side_Cluster > > > tOldSetClusters = mSideSets;

        // clear member data
        mSideSideSetLabelToOrd.clear();
        mSideSetLabels.clear();
        mSideSets.clear();

        for ( auto iB : mListOfSideSets )
        {
            delete iB;
        }
        mListOfSideSets.clear();

        // build list of side set indices
        Vector< moris_index > tSideSetIndexList;
        tSideSetIndexList.reserve( tOldSetClusters.size() );

        // current index
        moris_index tSetIndex = 0;

        for ( uint i = 0; i < tOldSetClusters.size(); i++ )
        {
            if ( tSideSetsToKeepMap.find( tOldSetNames( i ) ) != tSideSetsToKeepMap.end() )
            {
                mSideSetLabels.push_back( tOldSetNames( i ) );
                mSideSets.push_back( tOldSetClusters( i ) );

                MORIS_ASSERT( mSideSideSetLabelToOrd.find( tOldSetNames( i ) ) == mSideSideSetLabelToOrd.end(),
                        "Duplicate block set in mesh" );

                mSideSideSetLabelToOrd[ tOldSetNames( i ) ] = tSetIndex;

                tSideSetIndexList.push_back( tSetIndex );

                tSetIndex++;
            }
        }

        // set side set indices
        this->commit_side_set( tSideSetIndexList );
        this->communicate_sets_of_type( mtk::SetType::SIDESET );
    }

    //------------------------------------------------------------------------------

    moris::Memory_Map
    Enriched_Integration_Mesh::get_memory_usage()
    {
        // memory map of ig mesh
        moris::Memory_Map tMM;
        tMM.mMemoryMapData[ "mVertexSetNames" ]              = moris::internal_capacity( mVertexSetNames );
        tMM.mMemoryMapData[ "mVerticesInVertexSet" ]         = moris::internal_capacity( mVerticesInVertexSet );
        tMM.mMemoryMapData[ "mVertexSetColors" ]             = moris::internal_capacity( mVertexSetColors );
        tMM.mMemoryMapData[ "mBlockSetNames" ]               = moris::internal_capacity( mBlockSetNames );
        tMM.mMemoryMapData[ "mBlockSetTopology" ]            = mBlockSetTopology.capacity();
        tMM.mMemoryMapData[ "mBlockSetNames" ]               = moris::internal_capacity( mBlockSetNames );
        tMM.mMemoryMapData[ "mPrimaryBlockSetClusters" ]     = moris::internal_capacity( mPrimaryBlockSetClusters );
        tMM.mMemoryMapData[ "mBlockSetColors" ]              = moris::internal_capacity( mBlockSetColors );
        tMM.mMemoryMapData[ "mColorsBlockSets" ]             = moris::internal_capacity( mColorsBlockSets );
        tMM.mMemoryMapData[ "mSideSetLabels" ]               = moris::internal_capacity( mSideSetLabels );
        tMM.mMemoryMapData[ "mSideSets" ]                    = moris::internal_capacity_nested_ptr( mSideSets );
        tMM.mMemoryMapData[ "mSideSetColors" ]               = moris::internal_capacity( mSideSetColors );
        tMM.mMemoryMapData[ "mColorsSideSets" ]              = moris::internal_capacity( mColorsSideSets );
        tMM.mMemoryMapData[ "mDoubleSideSetLabels" ]         = moris::internal_capacity( mDoubleSideSetLabels );
        tMM.mMemoryMapData[ "mDoubleSideSets" ]              = moris::internal_capacity( mDoubleSideSets );
        tMM.mMemoryMapData[ "mDoubleSideSetsLeaderIndex" ]   = moris::internal_capacity( mDoubleSideSetsLeaderIndex );
        tMM.mMemoryMapData[ "mDoubleSideSetsFollowerIndex" ] = moris::internal_capacity( mDoubleSideSetsFollowerIndex );
        // FIXME: Implement capacities down through MTK children
        //     tMM.mMemoryMapData["mDoubleSideClusters"] = moris::internal_capacity(mDoubleSideClusters);
        tMM.mMemoryMapData[ "mDoubleSideSingleSideClusters" ] = moris::internal_capacity_ptr( mDoubleSideSingleSideClusters );
        tMM.mMemoryMapData[ "mBulkPhaseToDblSideIndex" ]      = mBulkPhaseToDblSideIndex.capacity();
        tMM.mMemoryMapData[ "mLeaderDoubleSideSetColor" ]     = moris::internal_capacity( mLeaderDoubleSideSetColor );
        tMM.mMemoryMapData[ "mFollowerDoubleSideSetColor" ]   = moris::internal_capacity( mFollowerDoubleSideSetColor );
        tMM.mMemoryMapData[ "mColorLeaderDoubleSideSet" ]     = moris::internal_capacity( mColorLeaderDoubleSideSet );
        tMM.mMemoryMapData[ "mColorFollowerDoubleSideSet" ]   = moris::internal_capacity( mColorFollowerDoubleSideSet );
        return tMM;
    }

    //------------------------------------------------------------------------------

    mtk::CellTopology
    Enriched_Integration_Mesh::get_blockset_topology( const std::string &aSetName )
    {
        moris_index tIndex = this->get_block_set_index( aSetName );
        return mBlockSetTopology( tIndex );
    }

    //------------------------------------------------------------------------------

    mtk::CellShape
    Enriched_Integration_Mesh::get_IG_blockset_shape( const std::string &aSetName )
    {
        // get the clusters in the set
        Vector< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        mtk::CellShape tCellShape = mtk::CellShape::EMPTY;

        // get the number of clusters on the set
        uint tNumClustersInSet = tSetClusters.size();

        // if the set isn't empty, it exist
        if ( tNumClustersInSet > 0 )
        {
            // initialize variables for the while-loop below
            bool tIsVoidCluster = true;

            for ( uint iTestCluster = 0; iTestCluster < tNumClustersInSet; iTestCluster++ )
            {
                // get the cells in the current cluster
                Vector< moris::mtk::Cell const * > tClusterCells = tSetClusters( iTestCluster )->get_primary_cells_in_cluster();

                // check if it is a void cluster
                tIsVoidCluster = ( tClusterCells.size() == 0 );

                // if it is not a void cluster, use cluster for reference for cell shape
                if ( !tIsVoidCluster )
                {
                    tCellShape = tClusterCells( 0 )->get_cell_info()->compute_cell_shape( tClusterCells( 0 ) );
                    break;
                }
            }

            // if it made it here without finding a non-empty cluster, something's weird
            MORIS_ERROR( tCellShape != mtk::CellShape::EMPTY,
                    "Enriched_Integration_Mesh::get_IG_blockset_shape() - No non-void clusters found on Set: %s",
                    aSetName.c_str() );
        }

        // within debug, checking all cells to make sure that they are the same Cell Shape
        // if cells exist

        // looping through the clusters
        for ( uint iCluster = 0; iCluster < tSetClusters.size(); iCluster++ )
        {
            // get cell of cells in the cluster
            Vector< moris::mtk::Cell const * > tClusterCellsCheck = tSetClusters( iCluster )->get_primary_cells_in_cluster();

            // looping through the cells in the cluster
            for ( uint iCheckCell = 0; iCheckCell < tClusterCellsCheck.size(); iCheckCell++ )
            {
                MORIS_ASSERT( tClusterCellsCheck( iCheckCell )->get_cell_info()->compute_cell_shape( tClusterCellsCheck( iCheckCell ) ) == tCellShape,
                        "Enriched_Integration_Mesh::get_IG_blockset_shape() - cell shape is not consistent in the block" );
            }
        }

        return tCellShape;
    }

    //------------------------------------------------------------------------------

    mtk::CellShape
    Enriched_Integration_Mesh::get_IP_blockset_shape( const std::string &aSetName )
    {
        // get the clusters in the set
        Vector< mtk::Cluster const * > tSetClusters = this->get_set_by_name( aSetName )->get_clusters_on_set();

        // init cell shape
        mtk::CellShape tCellShape = mtk::CellShape::EMPTY;

        // if the set isn't empty exist
        if ( tSetClusters.size() > 0 )
        {
            // get the cells in the first cluster
            mtk::Cell const &tClusterCell = tSetClusters( 0 )->get_interpolation_cell();

            // compute the cell shape of the first cell
            tCellShape = tClusterCell.get_cell_info()->compute_cell_shape( &tClusterCell );
        }

        // within debug, checking all cells to make sure that they are the same Cell Shape
        // if cells exist
        // looping through the clusters
        for ( uint iCluster = 1; iCluster < tSetClusters.size(); iCluster++ )
        {
            MORIS_ASSERT( tSetClusters( iCluster )->get_interpolation_cell().get_cell_info()->compute_cell_shape( &tSetClusters( iCluster )->get_interpolation_cell() )
                                  == tCellShape,
                    "Enriched_Integration_Mesh::get_IP_blockset_shape - cell shape is not consistent in the block" );
        }

        return tCellShape;
    }

    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Enriched_Integration_Mesh::convert_indices_to_ids(
            Matrix< IndexMat > const &aIndices,
            mtk::EntityRank           aEntityRank ) const
    {
        uint            tNRow = aIndices.n_rows();
        uint            tNCol = aIndices.n_cols();
        Matrix< IdMat > tIds( tNRow, tNCol );
        for ( uint i = 0; i < tNRow; i++ )
        {
            for ( uint j = 0; j < tNCol; j++ )
            {
                tIds( i, j ) = this->get_glb_entity_id_from_entity_loc_index( aIndices( i, j ), aEntityRank );
            }
        }
        return tIds;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::convert_ids_to_indices(
            Matrix< IdMat > const &aIds,
            mtk::EntityRank        aEntityRank ) const
    {
        uint            tNRow = aIds.n_rows();
        uint            tNCol = aIds.n_cols();
        Matrix< IdMat > tIndices( tNRow, tNCol );
        for ( uint i = 0; i < tNRow; i++ )
        {
            for ( uint j = 0; j < tNCol; j++ )
            {
                tIndices( i, j ) = this->get_loc_entity_ind_from_entity_glb_id( aIds( i, j ), aEntityRank );
            }
        }

        return tIndices;
    }

    //------------------------------------------------------------------------------

    Vector< moris::mtk::Cell const * >
    Enriched_Integration_Mesh::get_mtk_cells_loc_inds( Matrix< IndexMat > const &aCellIndices )
    {
        uint                               tNumCells = aCellIndices.numel();
        Vector< moris::mtk::Cell const * > tCells( tNumCells );

        for ( uint i = 0; i < tNumCells; i++ )
        {
            tCells( i ) = &this->get_mtk_cell( aCellIndices( i ) );
        }

        return tCells;
    }
    //------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Enriched_Integration_Mesh::get_mtk_vertices_loc_inds( Matrix< IndexMat > const &aVertexIndices )
    {
        uint                                 tNumVerts = aVertexIndices.numel();
        Vector< moris::mtk::Vertex const * > tVertices( tNumVerts );

        for ( uint i = 0; i < tNumVerts; i++ )
        {
            tVertices( i ) = &this->get_mtk_vertex( (moris_index)aVertexIndices( i ) );
        }

        return tVertices;
    }

    //------------------------------------------------------------------------------

    xtk::Cell_Cluster const &
    Enriched_Integration_Mesh::get_xtk_cell_cluster( mtk::Cell const &aInterpCell ) const
    {
        return get_xtk_cell_cluster( aInterpCell.get_index() );
    }

    //------------------------------------------------------------------------------

    xtk::Cell_Cluster const &
    Enriched_Integration_Mesh::get_xtk_cell_cluster( moris_index aInterpCellIndex ) const
    {
        MORIS_ASSERT( aInterpCellIndex < (moris_index)mCellClusters.size(), "Interpolation Cell index out of bounds" );
        return *mCellClusters( aInterpCellIndex );
    }

    //------------------------------------------------------------------------------

    mtk::Cell_Cluster const &
    Enriched_Integration_Mesh::get_cell_cluster( mtk::Cell const &aInterpCell ) const
    {
        return get_cell_cluster( aInterpCell.get_index() );
    }
    //------------------------------------------------------------------------------

    mtk::Cell_Cluster const &
    Enriched_Integration_Mesh::get_cell_cluster( moris_index aInterpCellIndex ) const
    {
        MORIS_ASSERT( aInterpCellIndex < (moris_index)mCellClusters.size(), "Interpolation Cell index out of bounds" );
        return *mCellClusters( aInterpCellIndex );
    }
    //------------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Integration_Mesh::get_block_set_names() const
    {
        return mBlockSetNames;
    }

    //------------------------------------------------------------------------------

    std::string
    Enriched_Integration_Mesh::get_block_set_label( moris_index aBlockSetOrdinal ) const
    {
        MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mSideSetLabels.size(), "Block set ordinal out of bounds" );
        return mBlockSetNames( aBlockSetOrdinal );
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::get_block_set_index( std::string aBlockSetLabel ) const
    {
        auto tIter = mBlockSetLabelToOrd.find( aBlockSetLabel );

        MORIS_ERROR( tIter != mBlockSetLabelToOrd.end(), "block set set label not found" );

        return tIter->second;
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Cluster const * >
    Enriched_Integration_Mesh::get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const
    {
        MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mBlockSetNames.size(), "Requested block set ordinal out of bounds." );

        Vector< xtk::Cell_Cluster const * > const &tXTKClustersInSet = mPrimaryBlockSetClusters( aBlockSetOrdinal );

        Vector< mtk::Cluster const * > tClusterInSet( tXTKClustersInSet.size() );

        for ( uint i = 0; i < tXTKClustersInSet.size(); i++ )
        {
            tClusterInSet( i ) = tXTKClustersInSet( i );
        }

        return tClusterInSet;
    }
    //------------------------------------------------------------------------------

    Vector< xtk::Cell_Cluster const * > const &
    Enriched_Integration_Mesh::get_xtk_cell_clusters_in_block_set( moris_index aBlockSetOrdinal ) const
    {
        MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mBlockSetNames.size(), "Requested block set ordinal out of bounds." );

        return mPrimaryBlockSetClusters( aBlockSetOrdinal );
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_block_set_colors( moris_index aBlockSetOrdinal ) const
    {
        MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mBlockSetColors.size(), "Block set ordinal out of bounds" );
        return mBlockSetColors( aBlockSetOrdinal );
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Cluster const * >
    Enriched_Integration_Mesh::get_side_set_cluster( moris_index aSideSetOrdinal ) const
    {
        MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSets.size(), "Side set ordinal out of bounds" );

        uint tNumSideClustersInSet = mSideSets( aSideSetOrdinal ).size();

        Vector< mtk::Cluster const * > tSideClustersInSet( tNumSideClustersInSet );

        for ( uint i = 0; i < tNumSideClustersInSet; i++ )
        {
            tSideClustersInSet( i ) = mSideSets( aSideSetOrdinal )( i ).get();
        }

        return tSideClustersInSet;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_side_set_colors( moris_index aSideSetOrdinal ) const
    {
        MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSetColors.size(), "Side set ordinal out of bounds" );
        return mSideSetColors( aSideSetOrdinal );
    }
    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_side_sets() const
    {
        return mSideSets.size();
    }
    //------------------------------------------------------------------------------

    std::string
    Enriched_Integration_Mesh::get_side_set_label( moris_index aSideSetOrdinal ) const
    {
        MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSetLabels.size(), "Side set ordinal out of bounds" );
        return mSideSetLabels( aSideSetOrdinal );
    }

    //------------------------------------------------------------------------------
    moris_index
    Enriched_Integration_Mesh::get_side_set_index( std::string aSideSetLabel ) const
    {
        auto tIter = mSideSideSetLabelToOrd.find( aSideSetLabel );

        MORIS_ERROR( tIter != mSideSideSetLabelToOrd.end(), "side side set label not found" );

        return tIter->second;
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_double_sided_sets() const
    {
        return mDoubleSideSetLabels.size();
    }

    //------------------------------------------------------------------------------

    std::string
    Enriched_Integration_Mesh::get_double_sided_set_label( moris_index aSideSetOrdinal ) const
    {
        MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(), "Double side set ordinal out of bounds" );
        return mDoubleSideSetLabels( aSideSetOrdinal );
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::get_double_sided_set_index( std::string aDoubleSideSetLabel ) const
    {
        auto tIter = mDoubleSideSetLabelToOrd.find( aDoubleSideSetLabel );

        MORIS_ERROR( tIter != mDoubleSideSetLabelToOrd.end(), "double side set label not found" );

        return tIter->second;
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Cluster const * >
    Enriched_Integration_Mesh::get_double_side_set_cluster( moris_index aSideSetOrdinal ) const
    {
        MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(),
                "Enriched_Integration_Mesh::get_double_side_set_cluster() - Double side set ordinal out of bounds" );

        uint tNumDblSideClustersInSet = mDoubleSideSets( aSideSetOrdinal ).size();

        Vector< mtk::Cluster const * > tDblSideClusters( tNumDblSideClustersInSet );

        for ( uint i = 0; i < tNumDblSideClustersInSet; i++ )
        {
            tDblSideClusters( i ) = mDoubleSideSets( aSideSetOrdinal )( i ).get();
        }

        return tDblSideClusters;
    }

    Vector< std::string >
    Enriched_Integration_Mesh::get_field_names(
            mtk::EntityRank          aEntityRank,
            const moris::moris_index aSetOrdinal )
    {
        // initialize output
        Vector< std::string > tOutputFieldNames;

        // get the field storage index corresponding to the entity type (nodal, elemental, facts)
        moris_index tRankFieldIndex = this->get_entity_rank_field_index( aEntityRank );

        // return fields existing globally if requested
        if ( aSetOrdinal == MORIS_INDEX_MAX )
        {
            // collect the field labels
            for ( auto const &iter : mGlobalSetFieldLabelToIndex( tRankFieldIndex ) )
            {
                tOutputFieldNames.push_back( iter.first );
            }
        }

        // return fields constructed set-wise
        else
        {
            if ( mSetWiseFieldLabelToIndex( tRankFieldIndex ).size() == 0 )
            {
                return tOutputFieldNames;
            }

            // collect the field labels
            for ( auto const &iter : mSetWiseFieldLabelToIndex( tRankFieldIndex )( aSetOrdinal ) )
            {
                tOutputFieldNames.push_back( iter.first );
            }
        }

        // return the
        return tOutputFieldNames;
    }

    //------------------------------------------------------------------------------

    moris::moris_index
    Enriched_Integration_Mesh::create_field(
            std::string        aLabel,
            mtk::EntityRank    aEntityRank,
            moris::moris_index aSetOrdinal )
    {
        // make sure there are no redundant field labels/names
        MORIS_ASSERT( !field_exists( aLabel, aEntityRank, aSetOrdinal ), "Enriched_Integration_Mesh::create_field() - Field already created." );

        // get the field index for the requested entity type
        moris_index tEntityFieldIndex = this->get_entity_rank_field_index( aEntityRank );

        // initialize field index
        moris::moris_index tFieldIndex;

        // create a global field
        if ( aSetOrdinal == MORIS_INDEX_MAX )
        {
            // get first free index is index for the new field
            tFieldIndex = mFields.size();

            // store field label in map
            mGlobalSetFieldLabelToIndex( tEntityFieldIndex )[ aLabel ] = tFieldIndex;

            // create and save field
            mFields.push_back( Field( aLabel, aSetOrdinal ) );
        }

        // create set-local field
        else
        {
            // if the array of side set fields has not been initialized yet, do so here
            if ( mSideSetFields.size() < this->get_num_side_sets() )
            {
                mSideSetFields.resize( this->get_num_side_sets() );
            }

            // get first free index is index for the new field
            tFieldIndex = mSideSetFields( aSetOrdinal ).size();

            // resize the map of set-wise fields if necessary // TODO: this should be done exactly once, when it is known how many sets/blocks there are
            uint tMinSize = (uint)aSetOrdinal + 1;
            if ( mSetWiseFieldLabelToIndex( tEntityFieldIndex ).size() < tMinSize )
            {
                mSetWiseFieldLabelToIndex( tEntityFieldIndex ).resize( tMinSize );
            }

            // store field name in map
            mSetWiseFieldLabelToIndex( tEntityFieldIndex )( aSetOrdinal )[ aLabel ] = tFieldIndex;

            // store set-wise field
            mSideSetFields( aSetOrdinal ).push_back( Field( aLabel, aSetOrdinal ) );
        }

        // return index of the field on the set or globally
        return tFieldIndex;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Enriched_Integration_Mesh::get_double_side_set_colors( moris_index aSideSetOrdinal ) const
    {
        MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(), "Double side set ordinal out of bounds" );
        return mLeaderDoubleSideSetColor( aSideSetOrdinal );
    }
    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_sidesets_num_faces( Vector< moris_index > aSideSetIndex ) const
    {
        uint tNumSideSetFaces = 0;

        for ( moris_index Ik = 0; Ik < (moris_index)aSideSetIndex.size(); ++Ik )
        {
            MORIS_ASSERT( aSideSetIndex( Ik ) < (moris_index)mSideSets.size(), "Side set index out of bounds" );

            // add up the sideset number of faces
            tNumSideSetFaces = tNumSideSetFaces + mSideSets( aSideSetIndex( Ik ) ).size();
        }

        return tNumSideSetFaces;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print() const
    {
        this->print_general();
        this->print_block_sets();
        this->print_side_sets();
        this->print_double_side_sets();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print_general() const
    {
        uint tNumCells = this->get_num_owned_cells();

        uint tNumGlobalCells = sum_all( tNumCells );
        if ( par_rank() == 0 )
        {
            std::cout << "Num Cells: " << std::setw( 8 ) << tNumGlobalCells << std::endl;
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print_cell_clusters( uint aVerbosityLevel ) const
    {
        std::cout << "\nCell Clusters:" << std::endl;
        for ( uint i = 0; i < mCellClusters.size(); i++ )
        {
            xtk::Cell_Cluster *tCluster    = mCellClusters( i ).get();
            std::string        tTrivialStr = "f";
            if ( tCluster->is_trivial() )
            {
                tTrivialStr = "t";
            }

            Interpolation_Cell_Unzipped const *tBaseInterpCell = tCluster->get_xtk_interpolation_cell();

            std::cout << "    Cluster Index: " << std::setw( 9 ) << i
                      << " | Interp Cell Id: " << std::setw( 9 ) << tCluster->get_interpolation_cell_id()
                      << " | Base Interp Cell Id: " << std::setw( 9 ) << tBaseInterpCell->get_base_cell()->get_id()
                      << " | Trivial: " << tTrivialStr
                      << " | Num Primary: " << std::setw( 9 ) << tCluster->get_num_primary_cells()
                      << " | Num Void: " << tCluster->get_num_void_cells()
                      << std::endl;

            moris::print( tCluster->get_vertices_local_coordinates_wrt_interp_cell(), "Local Coords" );

            if ( aVerbosityLevel > 0 )
            {
                Vector< moris::mtk::Cell const * > const &tPrimaryCells = tCluster->get_primary_cells_in_cluster();
                std::cout << "        Primary Integration Cell Ids: ";
                for ( uint i = 0; i < tCluster->get_num_primary_cells(); i++ )
                {
                    std::cout << std::setw( 9 ) << tPrimaryCells( i )->get_id();
                }
                std::cout << "\n"
                          << std::endl;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print_block_sets( uint aVerbosityLevel ) const
    {
        std::cout << "\nBlock Sets:" << std::endl;
        std::cout << "    Num Block Sets: " << this->get_num_blocks();

        for ( uint iBS = 0; iBS < this->get_num_blocks(); iBS++ )
        {
            std::cout << "\n    Block Name: " << std::setw( 20 ) << mBlockSetNames( iBS ) <<                   //
                    " | Block Set Ord: " << std::setw( 9 ) << iBS <<                                           //
                    " | Num Cell Clusters: " << std::setw( 9 ) << mPrimaryBlockSetClusters( iBS ).size() <<    //
                    " | Bulk Phase: " << std::setw( 9 ) << mBlockSetColors( iBS )( 0 );

            if ( aVerbosityLevel > 0 )
            {
                Vector< xtk::Cell_Cluster const * > tClusters = this->mPrimaryBlockSetClusters( iBS );
                std::cout << "\n            Cluster in set\n";
                for ( uint i = 0; i < tClusters.size(); i++ )
                {
                    xtk::Cell_Cluster const *tCluster    = tClusters( i );
                    std::string              tTrivialStr = "f";
                    if ( tCluster->is_trivial() )
                    {
                        tTrivialStr = "t";
                    }

                    Interpolation_Cell_Unzipped const *tBaseInterpCell = tCluster->get_xtk_interpolation_cell();

                    std::cout << "            Cluster Index: " << std::setw( 9 ) << i
                              << " | Interp Cell Id: " << std::setw( 9 ) << tCluster->get_interpolation_cell_id()
                              << " | Base Interp Cell Id: " << std::setw( 9 ) << tBaseInterpCell->get_base_cell()->get_id()
                              << " | Trivial: " << tTrivialStr
                              << " | Num Primary: " << std::setw( 9 ) << tCluster->get_num_primary_cells()
                              << " | Num Void: " << tCluster->get_num_void_cells()
                              << std::endl;
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print_side_sets( uint aVerbosityLevel ) const
    {
        std::cout << "\nSide Sets:" << std::endl;
        std::cout << "    Num Side Sets: " << this->get_num_side_sets() << std::endl;

        for ( uint iSS = 0; iSS < this->get_num_side_sets(); iSS++ )
        {
            std::cout << "    Side Set Name: " << std::setw( 20 ) << mSideSetLabels( iSS ) <<         //
                    " | Side Set Ord: " << std::setw( 9 ) << iSS <<                                   //
                    " | Num Cell Clusters: " << std::setw( 9 ) << this->mSideSets( iSS ).size() <<    //
                    " | Bulk Phase: " << std::setw( 9 ) << mSideSetColors( iSS )( 0 ) << std::endl;

            if ( aVerbosityLevel > 0 )
            {
                for ( uint iCL = 0; iCL < this->mSideSets( iSS ).size(); iCL++ )
                {
                    const std::shared_ptr< xtk::Side_Cluster > tCluster = this->mSideSets( iSS )( iCL );

                    Vector< moris::mtk::Cell const * > const &tPrimaryCells = tCluster->get_primary_cells_in_cluster();

                    for ( uint i = 0; i < tCluster->get_num_primary_cells(); i++ )
                    {
                        std::cout << "Integration Cell Id = " << std::setw( 9 ) << tPrimaryCells( i )->get_id() << std::endl;
                        moris::print( tCluster->get_cell_local_coords_on_side_wrt_interp_cell( i ), "Local Coords on Side" );
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print_double_side_sets( uint aVerbosityLevel ) const
    {
        std::cout << "\nDouble Side Sets:" << std::endl;
        std::cout << "    Num Side Sets: " << this->get_num_double_side_set() << std::endl;

        for ( uint iSS = 0; iSS < this->get_num_double_side_set(); iSS++ )
        {
            std::cout << "    Dbl Side Set Name: " << std::setw( 20 ) << mDoubleSideSetLabels( iSS ) <<       //
                    " | Dbl Side Set Ord: " << std::setw( 9 ) << iSS <<                                       //
                    " | Num Cell Clusters: " << std::setw( 9 ) << this->mDoubleSideSets( iSS ).size() <<      //
                    " | Leader Bulk Phase: " << std::setw( 9 ) << mLeaderDoubleSideSetColor( iSS )( 0 ) <<    //
                    " | Follower Bulk Phase: " << std::setw( 9 ) << mFollowerDoubleSideSetColor( iSS )( 0 );

            if ( aVerbosityLevel > 0 )
            {
                for ( uint i = 0; i < mDoubleSideSets( iSS ).size(); i++ )
                {
                    std::cout << "\n      Leader Interpolation Cell: " << std::setw( 9 ) <<    //
                            mDoubleSideSets( iSS )( i )->get_interpolation_cell( mtk::Leader_Follower::LEADER ).get_id();
                    std::cout << " | Follower Interpolation Cell: " << std::setw( 9 ) <<       //
                            mDoubleSideSets( iSS )( i )->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER ).get_id();
                }
            }

            std::cout << std::endl;
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::print_double_side_clusters( uint aVerbosityLevel ) const
    {
        std::cout << "\nDouble Side Clusters:" << std::endl;
        std::cout << "    Num Double Side Clusters: " << mDoubleSideClusters.size() << std::endl;

        for ( uint i = 0; i < mDoubleSideClusters.size(); i++ )
        {
            std::cout << mDoubleSideClusters( i ) << std::endl;
        }
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::create_side_set_from_dbl_side_set(
            moris_index const &aDblSideSetIndex,
            std::string const &aSideSetName,
            bool               aCollectSets )
    {
        Vector< moris_index > tSideSetIndex = this->register_side_set_names( { aSideSetName } );

        Vector< std::shared_ptr< mtk::Double_Side_Cluster > > &tDblSideClusters = mDoubleSideSets( aDblSideSetIndex );

        uint tCount = 0;
        for ( uint i = 0; i < tDblSideClusters.size(); i++ )
        {
            // get the index
            moris_index tLeaderIndex   = mDoubleSideSetsLeaderIndex( aDblSideSetIndex )( i );
            moris_index tFollowerIndex = mDoubleSideSetsFollowerIndex( aDblSideSetIndex )( i );

            mSideSets( tSideSetIndex( 0 ) ).push_back( mDoubleSideSingleSideClusters( tLeaderIndex ) );
            mSideSets( tSideSetIndex( 0 ) ).push_back( mDoubleSideSingleSideClusters( tFollowerIndex ) );
            tCount++;
        }

        this->commit_side_set( tSideSetIndex( 0 ) );
        this->communicate_sets_of_type( mtk::SetType::SIDESET );

        this->set_side_set_colors( tSideSetIndex( 0 ), this->get_double_side_set_colors( aDblSideSetIndex ) );

        if ( aCollectSets )
        {
            this->setup_color_to_set();
            this->collect_all_sets( false );
        }

        return tSideSetIndex( 0 );
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::create_block_set_from_cells_of_side_set(
            moris_index const       &aSideSetIndex,
            std::string const       &aBlockSetName,
            mtk::CellTopology const &aCellTopo,
            bool                     aCreateOnlyForVis )
    {
        Vector< std::shared_ptr< xtk::Side_Cluster > > &tSideClusters = mSideSets( aSideSetIndex );

        Vector< moris_index > tBlockSetIndex = this->register_block_set_names_with_cell_topo( { aBlockSetName }, aCellTopo );

        std::unordered_map< moris_index, moris_index > tIpCellInSet;

        for ( uint i = 0; i < tSideClusters.size(); i++ )
        {
            // cast to xtk side cluster
            xtk::Side_Cluster *tSideCluster = tSideClusters( i ).get();

            if ( tIpCellInSet.find( tSideCluster->mIntegrationCells( 0 )->get_id() ) == tIpCellInSet.end() )
            {
                // create a new cell cluster
                std::shared_ptr< xtk::Cell_Cluster > tCellCluster = std::make_shared< Cell_Cluster >( aCreateOnlyForVis );

                // get the ith enriched interpolation cell
                tCellCluster->mInterpolationCell = tSideCluster->mInterpolationCell;

                // mark as trivial
                tCellCluster->mTrivial = true;

                // add interp cell as integration cell
                tCellCluster->mPrimaryIntegrationCells.append( tSideCluster->mIntegrationCells );

                mCellClusters.push_back( tCellCluster );

                mPrimaryBlockSetClusters( tBlockSetIndex( 0 ) ).push_back( tCellCluster.get() );

                tIpCellInSet[ tSideCluster->mIntegrationCells( 0 )->get_id() ] = i;
            }
        }

        this->commit_block_set( tBlockSetIndex( 0 ) );
        this->set_block_set_colors( tBlockSetIndex( 0 ), this->get_side_set_colors( aSideSetIndex ) );
        this->setup_color_to_set();
        this->collect_all_sets( false );

        // communicate block sets after committing a new one
        this->communicate_sets_of_type( mtk::SetType::BULK );

        return tBlockSetIndex( 0 );
    }

    //------------------------------------------------------------------------------

    std::string
    Enriched_Integration_Mesh::get_interface_side_set_name(
            moris_index aGeomIndex,
            moris_index aBulkPhaseIndex0,
            moris_index aBulkPhaseIndex1 )
    {
        MORIS_ASSERT( aGeomIndex < (moris_index)mModel->get_geom_engine()->get_number_of_geometries(), "Geometry index out of bounds" );
        MORIS_ASSERT( aBulkPhaseIndex0 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 0 out of bounds" );
        MORIS_ASSERT( aBulkPhaseIndex1 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 1 out of bounds" );

        return "iside_b0_" + std::to_string( aBulkPhaseIndex0 ) + "_b1_" + std::to_string( aBulkPhaseIndex1 );
    }

    //------------------------------------------------------------------------------

    std::string
    Enriched_Integration_Mesh::get_dbl_interface_side_set_name(
            moris_index aBulkPhaseIndex0,
            moris_index aBulkPhaseIndex1 )
    {
        MORIS_ASSERT( aBulkPhaseIndex0 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 0 out of bounds" );
        MORIS_ASSERT( aBulkPhaseIndex1 < (moris_index)mModel->get_geom_engine()->get_num_bulk_phase(), "Bulk phase index 1 out of bounds" );

        return "dbl_iside_p0_" + std::to_string( aBulkPhaseIndex0 ) + "_p1_" + std::to_string( aBulkPhaseIndex1 );
    }

    //------------------------------------------------------------------------------

    moris::moris_index
    Enriched_Integration_Mesh::get_field_index(
            const std::string        aLabel,
            const mtk::EntityRank    aEntityRank,
            const moris::moris_index aSetOrdinal )
    {
        // first check that the field exists
        MORIS_ASSERT( field_exists( aLabel, aEntityRank, aSetOrdinal ),
                "Enriched_Integration_Mesh::get_field_index() - Field does not exist in mesh" );

        // get the index in list
        moris_index tIndex = get_entity_rank_field_index( aEntityRank );

        // t
        if ( aSetOrdinal == MORIS_INDEX_MAX )
        {
            auto tIter = mGlobalSetFieldLabelToIndex( tIndex ).find( aLabel );
            return tIter->second;
        }

        // fields local to certain sets
        else
        {
            auto tIter = mSetWiseFieldLabelToIndex( tIndex )( aSetOrdinal ).find( aLabel );
            return tIter->second;
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::add_field_data(
            moris::moris_index      aFieldIndex,
            mtk::EntityRank         aEntityRank,
            Matrix< DDRMat > const &aFieldData,
            moris::moris_index      aSetOrdinal )
    {
        // if field is faceted
        if ( aEntityRank == this->get_facet_rank() )
        {
            mSideSetFields( aSetOrdinal )( aFieldIndex ).mFieldData = aFieldData.copy();
        }

        // if field is global, elemental, or nodal
        else
        {
            mFields( aFieldIndex ).mFieldData = aFieldData.copy();
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > const &
    Enriched_Integration_Mesh::get_field_data(
            moris::moris_index aFieldIndex,
            mtk::EntityRank    aEntityRank,
            moris::moris_index aSetOrdinal ) const
    {
        // if field is faceted
        if ( aEntityRank == this->get_facet_rank() )
        {
            return mSideSetFields( aSetOrdinal )( aFieldIndex ).mFieldData;
        }

        // if field is global, elemental, or nodal
        else
        {
            return mFields( aFieldIndex ).mFieldData;
        }
    }

    //------------------------------------------------------------------------------

    moris_id
    Enriched_Integration_Mesh::allocate_entity_ids(
            moris::size_t   aNumReqs,
            mtk::EntityRank aEntityRank )
    {
        MORIS_ASSERT( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
                "Enriched_Integration_Mesh::allocate_entity_ids() - Only Elements or Nodes have ids" );

        moris_id tGlobalMax = this->get_max_entity_id( aEntityRank );

        int tProcRank = par_rank();
        int tProcSize = par_size();

        Vector< moris_id > aGatheredInfo;
        Vector< moris_id > tFirstId( 1 );
        Vector< moris_id > tNumIdsRequested( 1 );

        tNumIdsRequested( 0 ) = (moris_id)aNumReqs;

        moris::gather( tNumIdsRequested, aGatheredInfo );

        Vector< moris_id > tProcFirstID( tProcSize );

        if ( tProcRank == 0 )
        {
            // Loop over entities print the number of entities requested by each processor
            for ( int iProc = 0; iProc < tProcSize; ++iProc )
            {
                // Give each processor their desired amount of IDs
                tProcFirstID( iProc ) = tGlobalMax;

                // Increment the first available node ID
                tGlobalMax = tGlobalMax + aGatheredInfo( iProc );
            }
        }

        moris::scatter( tProcFirstID, tFirstId );

        return tFirstId( 0 );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::commit_double_side_set( moris_index const &aDoubleSideSetIndex )
    {

        MORIS_ASSERT( mListOfDoubleSideSets.size() == (uint)aDoubleSideSetIndex,
                "Enriched_Integration_Mesh::commit_double_side_set() - "
                "Committing double side set failed. aDoubleSideSetIndex needs to be equivalent to the size of the list of double side sets" );

        mListOfDoubleSideSets.resize( mListOfDoubleSideSets.size() + 1, nullptr );

        mListOfDoubleSideSets( aDoubleSideSetIndex ) =
                new moris::mtk::Double_Side_Set( mDoubleSideSetLabels( aDoubleSideSetIndex ),
                        this->get_double_side_set_cluster( aDoubleSideSetIndex ),
                        this->get_double_side_set_colors( aDoubleSideSetIndex ),
                        this->get_spatial_dim() );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::commit_double_side_set( const Vector< moris_index > &aSideSetIndexList )
    {
        mListOfDoubleSideSets.resize( mListOfDoubleSideSets.size() + aSideSetIndexList.size(), nullptr );

        for ( uint iDblSS = 0; iDblSS < aSideSetIndexList.size(); ++iDblSS )
        {
            const moris_index tSideSetIndex = aSideSetIndexList( iDblSS );

            MORIS_ASSERT( (uint)tSideSetIndex < mListOfDoubleSideSets.size(),
                    "Enriched_Integration_Mesh::commit_double_side_set() - "
                    "Committing double side set failed. aDoubleSideSetIndex needs to be equivalent to the size of the list of double side sets" );

            mListOfDoubleSideSets( tSideSetIndex ) =
                    new moris::mtk::Double_Side_Set(
                            mDoubleSideSetLabels( tSideSetIndex ),
                            this->get_double_side_set_cluster( tSideSetIndex ),
                            this->get_double_side_set_colors( tSideSetIndex ),
                            this->get_spatial_dim() );
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::commit_side_set( moris_index const &aSideSetIndex )
    {
        MORIS_ASSERT( mListOfSideSets.size() == (uint)aSideSetIndex,
                "Enriched_Integration_Mesh::commit_side_set() - "
                "Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of single side sets" );

        mListOfSideSets.resize( mListOfSideSets.size() + 1, nullptr );

        mListOfSideSets( aSideSetIndex ) =
                new moris::mtk::Side_Set(
                        mSideSetLabels( aSideSetIndex ),
                        this->get_side_set_cluster( aSideSetIndex ),
                        this->get_side_set_colors( aSideSetIndex ),
                        this->get_spatial_dim() );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::commit_side_set( const Vector< moris_index > &aSideSetIndexList )
    {
        mListOfSideSets.resize( mListOfSideSets.size() + aSideSetIndexList.size(), nullptr );

        for ( uint iI = 0; iI < aSideSetIndexList.size(); ++iI )
        {
            const moris_index tSideSetIndex = aSideSetIndexList( iI );

            MORIS_ASSERT( (uint)tSideSetIndex < mListOfSideSets.size(),
                    "Enriched_Integration_Mesh::commit_side_set() - "
                    "Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of single side sets" );

            mListOfSideSets( tSideSetIndex ) =
                    new moris::mtk::Side_Set(
                            mSideSetLabels( tSideSetIndex ),
                            this->get_side_set_cluster( tSideSetIndex ),
                            this->get_side_set_colors( tSideSetIndex ),
                            this->get_spatial_dim() );
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::commit_block_set( moris_index const &aBlockSetIndex )
    {
        MORIS_ASSERT(
                mListOfBlocks.size() == (uint)aBlockSetIndex,
                "Enriched_Integration_Mesh::commit_block_set() -  "
                "Committing side set failed. aSideSetIndex needs to be equivalent to the size of the list of double side sets" );

        mListOfBlocks.resize( mListOfBlocks.size() + 1, nullptr );

        mListOfBlocks( aBlockSetIndex ) =
                new moris::mtk::Block_Set(
                        mBlockSetNames( aBlockSetIndex ),
                        this->get_cell_clusters_in_set( aBlockSetIndex ),
                        this->get_block_set_colors( aBlockSetIndex ),
                        this->get_spatial_dim() );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::communicate_sets_of_type( const mtk::SetType aSetType )
    {
        switch ( aSetType )
        {
            case mtk::SetType::BULK:
            {
                // communicate block sets
                mtk::Set_Communicator tSetCommunicator( mListOfBlocks );
                break;
            }
            case mtk::SetType::SIDESET:
            {
                // communicate side sets
                mtk::Set_Communicator tSetCommunicator( mListOfSideSets );
                break;
            }
            case mtk::SetType::DOUBLE_SIDED_SIDESET:
            {
                // communicate double sided side sets
                mtk::Set_Communicator tSetCommunicator( mListOfDoubleSideSets );
                break;
            }
            default:
            {
                MORIS_ERROR( false,
                        "xtk::Enriched_Integration_Mesh::communicate_sets_of_type() - Unknown set type to communicate." );
            }

        }    // end switch: set type

    }        // end function: Enriched_Integration_Mesh::communicate_sets_of_type()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_cell_clusters()
    {
        // trace/log this function
        Tracer tTracer( "XTK", "Enriched Integration Mesh", "setup_cell_clusters", mModel->mVerboseLevel, 1 );

        // get pointer to enr. IP mesh
        Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( 0 );

        // Number of interpolation cells
        uint tNumInterpCells = tEnrInterpMesh->get_num_entities( mtk::EntityRank::ELEMENT );

        // allocate cell cluster member data
        mCellClusters.resize( tNumInterpCells, nullptr );

        // Allocate subphase index to cluster index
        mSubphaseIndexToClusterIndex.resize( 1, tNumInterpCells );

        // reference the enriched cells
        Vector< Interpolation_Cell_Unzipped * > const &tEnrichedInterpCells = tEnrInterpMesh->get_enriched_interpolation_cells();

        // iterate through interpolation cells to create cell clusters
        for ( uint i = 0; i < tNumInterpCells; i++ )
        {
            // index
            moris_index tInterpCellIndex = tEnrichedInterpCells( i )->get_index();    // TODO: where does the enr. IP cell's index get set, how does it differ from position in list?

            // create a new cell cluster
            mCellClusters( tInterpCellIndex ) = std::make_shared< Cell_Cluster >();

            // get the ith enriched interpolation cell
            mCellClusters( tInterpCellIndex )->mInterpolationCell = tEnrichedInterpCells( i );

            // add subphase index to cluster index to subphase index data
            mSubphaseIndexToClusterIndex( i ) = mCellClusters( tInterpCellIndex )->get_xtk_interpolation_cell()->get_subphase_index();

            // base cell
            moris::mtk::Cell const *tBaseInterpCell = mCellClusters( tInterpCellIndex )->mInterpolationCell->get_base_cell();

            // ask background mesh if the base cell has children (the opposite answer to this question is the trivial flag)
            mCellClusters( tInterpCellIndex )->mTrivial = !mCutIgMesh->parent_cell_has_children( tBaseInterpCell->get_index() );

            // if it has children get a pointer to the child mesh
            if ( !mCellClusters( tInterpCellIndex )->mTrivial )
            {
                // subphase index
                moris_index tProcSubphaseIndex = mCellClusters( tInterpCellIndex )->mInterpolationCell->get_subphase_index();

                std::shared_ptr< IG_Cell_Group > tSubphaseCells = mCutIgMesh->get_subphase_ig_cells( tProcSubphaseIndex );

                // number of subphase associated with the same background cell
                Vector< moris_index > const &tParentCellSubphases = mCutIgMesh->get_parent_cell_subphases( tBaseInterpCell->get_index() );

                // vertex group
                std::shared_ptr< IG_Vertex_Group > tVertexGroupForCluster = mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tBaseInterpCell->get_index() ) );

                Vector< std::shared_ptr< IG_Cell_Group > > tVoidSubphases;
                tVoidSubphases.reserve( tParentCellSubphases.size() - 1 );

                // get the subphases in the void region
                for ( uint iVoid = 0; iVoid < tParentCellSubphases.size(); iVoid++ )
                {
                    if ( tParentCellSubphases( iVoid ) != tProcSubphaseIndex )
                    {
                        tVoidSubphases.push_back( mCutIgMesh->get_subphase_ig_cells( tParentCellSubphases( iVoid ) ) );
                    }
                }

                mCellClusters( tInterpCellIndex )->set_primary_integration_cell_group( tSubphaseCells );
                mCellClusters( tInterpCellIndex )->set_void_integration_cell_groups( tVoidSubphases );
                mCellClusters( tInterpCellIndex )->set_ig_vertex_group( tVertexGroupForCluster );
            }

            // trivial case, the base of the enriched interpolation cell becomes the primary cell
            else
            {
                mCellClusters( tInterpCellIndex )->mPrimaryIntegrationCells.push_back( mCellClusters( tInterpCellIndex )->mInterpolationCell->get_base_cell() );
                Matrix< IndexMat > tVertexIndices                     = mCellClusters( tInterpCellIndex )->mPrimaryIntegrationCells( 0 )->get_vertex_inds();
                mCellClusters( tInterpCellIndex )->mVerticesInCluster = this->get_mtk_vertices_loc_inds( tVertexIndices );
                tBaseInterpCell->get_cell_info()->get_loc_coords_of_cell( mCellClusters( tInterpCellIndex )->mLocalCoords );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_cell_clusters_new()
    {
        // trace/log this function
        Tracer tTracer( "XTK", "Enriched Integration Mesh", "setup_cell_clusters", mModel->mVerboseLevel, 1 );

        // get pointer to the enrichment
        Enrichment const *tEnrichment = mModel->mEnrichment;

        // get number of base IP cells
        uint tNumBaseIpCells = mModel->mBackgroundMesh->get_num_elems();

        // get pointer to enr. IP mesh
        Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( 0 );

        // Number of interpolation cells
        uint tNumEnrInterpCells = tEnrInterpMesh->get_num_entities( mtk::EntityRank::ELEMENT );

        // perform sanity check
        MORIS_ASSERT( tEnrichment->get_num_enr_ip_cells() == tNumEnrInterpCells,
                "Enriched_Integration_Mesh::setup_cell_clusters_new() - Enrichment and enriched IP mesh report different number of enr. IP cells" );

        // allocate memory for cell cluster member data
        mCellClusters.resize( tNumEnrInterpCells );

        // Allocate memory for subphase index to cluster index map
        Vector< moris_index > tEmptyIndexList( 0 );
        mClusterIndexToSubphaseIndices.resize( tNumEnrInterpCells );    // FIXME: this map doesn't work yet

        // reference the enriched IP cells
        Vector< Interpolation_Cell_Unzipped * > const &tEnrichedInterpCells = tEnrInterpMesh->get_enriched_interpolation_cells();

        // loop over all base IP cells
        for ( uint iIpCell = 0; iIpCell < tNumBaseIpCells; iIpCell++ )
        {
            // convert to index
            moris_index tIpCellIndex = (moris_index)iIpCell;

            // get the number of SPs associated with the current IP cell
            Vector< moris_index > const &tSPsOnCell    = mCutIgMesh->get_parent_cell_subphases( tIpCellIndex );
            uint                         tNumSPsOnCell = tSPsOnCell.size();

            // check if the IP cell gets cut
            bool tAllClustersOnCellTrivial = !mCutIgMesh->parent_cell_has_children( tIpCellIndex );

            // loop over all Enr. IP cells living on the current base IP cell
            for ( uint iEnrIpCellOnBaseIpCell = 0; iEnrIpCellOnBaseIpCell < tEnrichment->get_num_unzippings_of_base_ip_cell( iIpCell ); iEnrIpCellOnBaseIpCell++ )
            {
                // get the index of the current Enr IP cell
                moris_index tEnrIpCellIndex = tEnrichment->get_enr_ip_cell_indices_on_base_ip_cell( iIpCell )( iEnrIpCellOnBaseIpCell );

                // create a new cell cluster object
                mCellClusters( tEnrIpCellIndex ) = std::make_shared< Cell_Cluster >();

                // get the enriched interpolation cell object and give it to the cluster for access
                mCellClusters( tEnrIpCellIndex )->mInterpolationCell = tEnrichedInterpCells( tEnrIpCellIndex );

                // check if the clusters associated with the current IP cell are generally trivial
                mCellClusters( tEnrIpCellIndex )->mTrivial = tAllClustersOnCellTrivial;

                // check if the cluster is void, i.e. if the cluster corresponds to one of the SPs on the IP cell
                if ( iEnrIpCellOnBaseIpCell < tNumSPsOnCell )    // material cluster
                {
                    mCellClusters( tEnrIpCellIndex )->mVoid = false;

                    // material cluster is only trivial if it is the only subphase on base IP cell, i.e. there are no triangulated child cells
                    mCellClusters( tEnrIpCellIndex )->mTrivial = tAllClustersOnCellTrivial;
                }
                else    // void clusters
                {
                    mCellClusters( tEnrIpCellIndex )->mVoid = true;

                    // all void clusters are trivial
                    mCellClusters( tEnrIpCellIndex )->mTrivial = true;
                }

                // Cluster is trivial (integration domain is single quadrilateral cell which is either void or full)
                if ( mCellClusters( tEnrIpCellIndex )->mTrivial )
                {
                    // get the base cell
                    moris::mtk::Cell const *tBaseCell = mCellClusters( tEnrIpCellIndex )->mInterpolationCell->get_base_cell();

                    // get the parametric coordinates of the vertices and store them in the cluster
                    tBaseCell->get_cell_info()->get_loc_coords_of_cell( mCellClusters( tEnrIpCellIndex )->mLocalCoords );

                    // get the indices of the vertices associated with the IP cell
                    Matrix< IndexMat > tVertexIndices = tBaseCell->get_vertex_inds();

                    // find the mtk::cells with the indices and store them in the cluster
                    mCellClusters( tEnrIpCellIndex )->mVerticesInCluster = this->get_mtk_vertices_loc_inds( tVertexIndices );

                    // store base cell as primary or void integration cell on cluster
                    if ( mCellClusters( tEnrIpCellIndex )->mVoid )    // cluster is void
                    {
                        mCellClusters( tEnrIpCellIndex )->mVoidIntegrationCells.push_back( tBaseCell );
                    }
                    else    // cluster is full
                    {
                        mCellClusters( tEnrIpCellIndex )->mPrimaryIntegrationCells.push_back( tBaseCell );

                        // sanity check for this case
                        MORIS_ASSERT( mCellClusters( tEnrIpCellIndex )->is_full(),
                                "Enriched_Integration_Mesh::setup_cell_clusters_new() - cluster is trivial and non-void, but not recognized as full" );
                    }
                }
                else    // Cluster is non-trivial and has TRIs/TETs
                {
                    // get the group of all IG mesh vertices on the current IP cell
                    std::shared_ptr< IG_Vertex_Group > tVertexGroupForCluster =
                            mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tIpCellIndex ) );

                    // store the IG vertices on the cluster
                    mCellClusters( tEnrIpCellIndex )->set_ig_vertex_group( tVertexGroupForCluster );

                    // initialize list of void IG cell groups
                    Vector< std::shared_ptr< IG_Cell_Group > > tVoidSubphases;
                    tVoidSubphases.reserve( tNumSPsOnCell - 1 );

                    // get the only SP index on cluster
                    moris_index tPrimarySpIndex = tSPsOnCell( iEnrIpCellOnBaseIpCell );

                    // get the  IG cell group in the cluster
                    std::shared_ptr< IG_Cell_Group > tIgCellGroupsInCluster = mCutIgMesh->get_subphase_ig_cells( tPrimarySpIndex );

                    // set primary IG cells in cluster
                    mCellClusters( tEnrIpCellIndex )->set_primary_integration_cell_group( tIgCellGroupsInCluster );

                    // get the subphases in the void region
                    for ( uint iVoid = 0; iVoid < tNumSPsOnCell; iVoid++ )
                    {
                        if ( tSPsOnCell( iVoid ) != tPrimarySpIndex )
                        {
                            tVoidSubphases.push_back( mCutIgMesh->get_subphase_ig_cells( tSPsOnCell( iVoid ) ) );
                        }
                    }

                    // store the void IG cells with the cluster
                    mCellClusters( tEnrIpCellIndex )->set_void_integration_cell_groups( tVoidSubphases );

                }    // end: construction of valid clusters
            }        // end: loop over enriched IP cells associated with the IP cell
        }            // end: loop over base IP cells
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_blockset_with_cell_clusters()
    {
        // get background mesh
        moris::mtk::Mesh &tBackgroundMesh = mModel->get_background_mesh();

        // enriched interpolation mesh
        Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( 0 );

        // my proc rank
        moris_index tProcRank = par_rank();

        // get block sets (in background mesh data)
        Vector< std::string > tBlockSetsNames = tBackgroundMesh.get_set_names( mtk::EntityRank::ELEMENT );

        // for each block set constructed
        for ( uint iBS = 0; iBS < tBlockSetsNames.size(); iBS++ )
        {
            // split set into child and no child as we need to have the same type of integration cell in each set
            Vector< std::string > tChildNoChildSetNames = this->split_set_name_by_child_no_child( tBlockSetsNames( iBS ) );

            // split child and no child sets by phases
            Vector< std::string > tPhaseChildBlockSetNames   = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 0 ) );
            Vector< std::string > tPhaseNoChildBlockSetNames = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 1 ) );

            // topology enums // TODO: move outside the loop
            mtk::CellTopology tChildTopo  = mModel->get_cut_integration_mesh()->get_child_element_topology();
            mtk::CellTopology tParentTopo = mModel->get_parent_cell_topology();

            // add block set names to member data
            Vector< moris_index > tChildBlockSetOrds   = this->register_block_set_names_with_cell_topo( tPhaseChildBlockSetNames, tChildTopo );
            Vector< moris_index > tNoChildBlockSetOrds = this->register_block_set_names_with_cell_topo( tPhaseNoChildBlockSetNames, tParentTopo );

            // set block set colors
            for ( moris_index iSet = 0; iSet < (moris_index)tChildBlockSetOrds.size(); iSet++ )
            {
                this->set_block_set_colors( tChildBlockSetOrds( iSet ), { { iSet } } );
                this->set_block_set_colors( tNoChildBlockSetOrds( iSet ), { { iSet } } );
            }

            // get the IP cells in this block
            Vector< moris::mtk::Cell const * > tCellsInBlock = tBackgroundMesh.get_block_set_cells( tBlockSetsNames( iBS ) );

            // get the enriched IP cells in this block
            Vector< xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsInBlock =
                    tEnrInterpMesh->get_enriched_cells_from_base_cells( tCellsInBlock );

            // iterate through and add cluster associated with enriched cell to block set
            for ( uint iEnrIpCell = 0; iEnrIpCell < tEnrichedCellsInBlock.size(); iEnrIpCell++ )
            {
                // get access to the current UIPC
                xtk::Interpolation_Cell_Unzipped const *tEnrIpCell = tEnrichedCellsInBlock( iEnrIpCell );

                // get the bulk phase
                moris_index tBulkPhaseIndex = tEnrIpCell->get_bulkphase_index();

                // check the primary sub-phase to check whether this is a void UIPC
                moris_index tPrimarySubPhase = tEnrIpCell->get_subphase_index();

                // get cluster associated with enriched cell
                xtk::Cell_Cluster const &tCluster = this->get_xtk_cell_cluster( tEnrIpCell->get_index() );

                if ( tEnrIpCell->get_owner() == tProcRank )
                {
                    // get the index for the current blockset (depends on bulk-phase index and cut or non-cut)
                    moris_index tSetOrd = MORIS_INDEX_MAX;

                    // only register non-void clusters, i.e. clusters that have a primary sub-phase
                    if ( tPrimarySubPhase > -1 )
                    {
                        // sort into block with non-cut clusters
                        if ( tCluster.is_trivial() )
                        {
                            tSetOrd = tNoChildBlockSetOrds( tBulkPhaseIndex );
                        }

                        // sort into block with cut clusters
                        else
                        {
                            tSetOrd = tChildBlockSetOrds( tBulkPhaseIndex );
                        }

                        // add cluster to Set
                        mPrimaryBlockSetClusters( tSetOrd ).push_back( &tCluster );
                    }
                }
            }
        }

        for ( uint iBlockSet = mListOfBlocks.size(); iBlockSet < mPrimaryBlockSetClusters.size(); iBlockSet++ )
        {
            this->commit_block_set( iBlockSet );
        }

        // communicate committed block sets
        this->communicate_sets_of_type( mtk::SetType::BULK );

    }    // end function: Enriched_Integration_Mesh::setup_blockset_with_cell_clusters()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_side_set_clusters()
    {
        // get data for easy access
        Enriched_Interpolation_Mesh *tEnrInterpMesh  = mModel->mEnrichedInterpMesh( 0 );
        moris::mtk::Mesh            &tBackgroundMesh = *mModel->mBackgroundMesh;
        Integration_Mesh_Generator   tIGMeshGen;

        // rank enum for facets
        mtk::EntityRank tFacetRank = mModel->mBackgroundMesh->get_facet_rank();

        // get side sets (in background mesh data)
        Vector< std::string > tSideSetNames = tBackgroundMesh.get_set_names( tFacetRank );

        tSideSetNames = mModel->check_for_and_remove_internal_seacas_side_sets( tSideSetNames );

        // my proc rank
        moris_index tParRank = par_rank();

        // access background facet to child facet connectivity
        Vector< std::shared_ptr< Vector< moris::moris_index > > > const &tBGFacetToChildFacet =
                mCutIgMesh->get_background_facet_to_child_facet_connectivity();

        // for each side set construct
        for ( uint iSS = 0; iSS < tSideSetNames.size(); iSS++ )
        {
            // split set into child and no child as we need to have the same type of integration cell in each set
            Vector< std::string > tChildNoChildSetNames = this->split_set_name_by_child_no_child( tSideSetNames( iSS ) );

            // split child and no child sets by phases
            Vector< std::string > tPhaseChildSideSetNames   = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 0 ) );
            Vector< std::string > tPhaseNoChildSideSetNames = this->split_set_name_by_bulk_phase( tChildNoChildSetNames( 1 ) );

            // add side set names to member data
            Vector< moris_index > tChildSideSetOrds   = this->register_side_set_names( tPhaseChildSideSetNames );
            Vector< moris_index > tNoChildSideSetOrds = this->register_side_set_names( tPhaseNoChildSideSetNames );

            // set side set colors
            for ( moris_index i = 0; i < (moris_index)tChildSideSetOrds.size(); i++ )
            {
                this->set_side_set_colors( tChildSideSetOrds( i ), { { i } } );
                this->set_side_set_colors( tNoChildSideSetOrds( i ), { { i } } );
            }

            // get the base IP cells in this side set and their side ordinals
            Vector< mtk::Cell const * > tCellsInSideSet;
            Matrix< IndexMat >          tCellOrdsInSideSet;

            tBackgroundMesh.get_sideset_cells_and_ords(
                    tSideSetNames( iSS ),
                    tCellsInSideSet,
                    tCellOrdsInSideSet );

            // iterate through base IP cells in side set
            for ( uint iC = 0; iC < tCellsInSideSet.size(); iC++ )
            {
                mtk::Cell const *tBaseCell  = tCellsInSideSet( iC );
                moris_index      tSideOrd   = tCellOrdsInSideSet( iC );
                moris_index      tSideIndex = tBackgroundMesh.get_entity_connected_to_entity_loc_inds(
                        tBaseCell->get_index(),
                        mtk::EntityRank::ELEMENT,
                        tBackgroundMesh.get_facet_rank() )( tSideOrd );

                // only place cluster's related to the background cells owned by current proc in sets
                if ( tBaseCell->get_owner() == tParRank )
                {
                    // get the enriched interpolation cells associated with base cell
                    Vector< xtk::Interpolation_Cell_Unzipped const * > tEnrichedCellsOfBaseCell =
                            tEnrInterpMesh->get_enriched_cells_from_base_cell( tBaseCell );

                    // get the subphase indices associated with the current parent cell
                    Vector< moris_index > const &tSubphasesWrtParentCell = mCutIgMesh->get_parent_cell_subphases( tBaseCell->get_index() );

                    // check some things out
                    for ( uint iSP = 0; iSP < tSubphasesWrtParentCell.size(); iSP++ )
                    {
                        MORIS_ASSERT( tSubphasesWrtParentCell( iSP ) == tEnrichedCellsOfBaseCell( iSP )->get_subphase_index(),
                                "Discrepancy in subphase index." );
                    }

                    // if this is a null ptr, the bg facet does not have any child facets associated with it,
                    // this indicates that the side set is on a facet not intersected by the geometry
                    if ( tBGFacetToChildFacet( tSideIndex ) != nullptr )
                    {
                        // for the base cell, get the associated vertices
                        std::shared_ptr< IG_Vertex_Group > tVertexGroupForCluster =
                                mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tBaseCell->get_index() ) );

                        // collect all the integration cells on the current facet and the side ordinal
                        Vector< moris::mtk::Cell * > tIgCellsOnBgFacet;
                        Vector< moris_index >        tIgCellsSideOrdsOnBgFacet;
                        tIGMeshGen.collect_ig_cells_and_side_ords_on_bg_facet( mCutIgMesh, tSideIndex, tIgCellsOnBgFacet, tIgCellsSideOrdsOnBgFacet );

                        // iterate through the subphases
                        for ( uint iSP = 0; iSP < tSubphasesWrtParentCell.size(); iSP++ )
                        {
                            Vector< moris::mtk::Cell const * > tIgCellsInSubphase;
                            Vector< moris_index >              tIgCellsSideOrdsInSubphase;
                            tIgCellsInSubphase.reserve( tIgCellsOnBgFacet.size() );
                            tIgCellsInSubphase.reserve( tIgCellsSideOrdsInSubphase.size() );

                            // bulk phase of this subphase
                            moris_index tBulkPhase = mCutIgMesh->get_subphase_bulk_phase( tSubphasesWrtParentCell( iSP ) );

                            // iterate through ig cells on bg facet
                            for ( uint iIGCell = 0; iIGCell < tIgCellsOnBgFacet.size(); iIGCell++ )
                            {
                                if ( mCutIgMesh->get_ig_cell_subphase_index( tIgCellsOnBgFacet( iIGCell )->get_index() ) == tSubphasesWrtParentCell( iSP ) )
                                {
                                    tIgCellsInSubphase.push_back( tIgCellsOnBgFacet( iIGCell ) );
                                    tIgCellsSideOrdsInSubphase.push_back( tIgCellsSideOrdsOnBgFacet( iIGCell ) );
                                }
                            }

                            // convert the side ordinal to a matrix
                            Matrix< IndexMat > tSideOrdsMat( 1, tIgCellsSideOrdsInSubphase.size() );
                            for ( uint iCopy = 0; iCopy < tIgCellsSideOrdsInSubphase.size(); iCopy++ )
                            {
                                tSideOrdsMat( iCopy ) = tIgCellsSideOrdsInSubphase( iCopy );
                            }

                            if ( tIgCellsInSubphase.size() > 0 )
                            {
                                // make a side cluster
                                std::shared_ptr< Side_Cluster > tSideCluster = std::make_shared< Side_Cluster >();
                                tSideCluster->mInterpolationCell             = tEnrichedCellsOfBaseCell( iSP );
                                tSideCluster->mTrivial                       = false;
                                tSideCluster->mIntegrationCellSideOrdinals   = tSideOrdsMat;
                                tSideCluster->mAssociatedCellCluster         = &this->get_xtk_cell_cluster( *tEnrichedCellsOfBaseCell( iSP ) );
                                tSideCluster->mIntegrationCells              = tIgCellsInSubphase;

                                tSideCluster->set_ig_vertex_group( tVertexGroupForCluster );

                                // add to the integration mesh
                                mSideSets( tChildSideSetOrds( tBulkPhase ) ).push_back( tSideCluster );
                            }
                        }
                    }
                    else
                    {
                        // phase of cell
                        moris_index tBulkPhase = tEnrichedCellsOfBaseCell( 0 )->get_bulkphase_index();

                        // side set ordinal in mesh
                        moris_index tSideSetOrd = tNoChildSideSetOrds( tBulkPhase );

                        // create a new side cluster in the side set associated with this bulk phase
                        std::shared_ptr< Side_Cluster > tSideCluster = std::make_shared< Side_Cluster >();
                        mSideSets( tSideSetOrd ).push_back( tSideCluster );

                        // set trivial flag
                        tSideCluster->mTrivial = true;

                        // get the set enriched interpolation cell
                        tSideCluster->mInterpolationCell = tEnrichedCellsOfBaseCell( 0 );

                        // integration cell is the same as the interpolation cell in this case
                        tSideCluster->mIntegrationCells = { tSideCluster->mInterpolationCell->get_base_cell() };

                        // side ordinal
                        tSideCluster->mIntegrationCellSideOrdinals = Matrix< IndexMat >( { { tSideOrd } } );

                        // add vertices
                        tSideCluster->mVerticesInCluster.append( tSideCluster->mIntegrationCells( 0 )->get_vertices_on_side_ordinal( tSideOrd ) );

                        tSideCluster->mAssociatedCellCluster = &this->get_xtk_cell_cluster( *tEnrichedCellsOfBaseCell( 0 ) );
                    }
                }
            }
        }

        // build list of side set indices
        Vector< moris_index > tSideSetIndexList;
        tSideSetIndexList.reserve( mSideSets.size() - mListOfSideSets.size() );

        for ( uint Ik = mListOfSideSets.size(); Ik < mSideSets.size(); Ik++ )
        {
            tSideSetIndexList.push_back( Ik );
        }

        // commit and communicate the side sets
        this->commit_side_set( tSideSetIndexList );
        this->communicate_sets_of_type( mtk::SetType::SIDESET );

    }    // end function: Enriched_Integration_Mesh::setup_side_set_clusters()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_double_side_set_clusters()
    {
        this->setup_double_sided_interface_sides();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_color_to_set()
    {
        this->construct_color_to_set_relationship( mBlockSetColors, mColorsBlockSets );
        this->construct_color_to_set_relationship( mSideSetColors, mColorsSideSets );
        this->construct_color_to_set_relationship( mLeaderDoubleSideSetColor, mColorLeaderDoubleSideSet );
        this->construct_color_to_set_relationship( mFollowerDoubleSideSetColor, mColorFollowerDoubleSideSet );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_cluster_groups()
    {
        // Trace this function
        Tracer tTracer( "XTK", "Enriched Integration Mesh", "Setup cluster groups" );

        this->setup_cell_cluster_groups();
        this->setup_dbl_side_cluster_groups();
        this->setup_side_cluster_groups();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_cell_cluster_groups()
    {
        // determine number of B-spline meshes
        uint tNumBspMeshes = mBsplineMeshIndices.numel();
        uint tMaxDMI       = (uint)mBsplineMeshIndices.max();

        // initialize list of Clusters
        mCellClusterGroups.resize( tMaxDMI + 1 );

        // get the SPG to Cluster Index map from the enrichment
        Vector< Vector< Vector< moris_index > > > const &tSpgToClusterIndex = mModel->mEnrichment->get_SPG_to_UIPC_map();

        // establish cluster groups for every B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the discretization mesh index (DMI)
            moris_index tDMI = mBsplineMeshIndices( iBspMesh );

            // get the number of SPGs on the current B-spline mesh
            uint tNumSPGs = tSpgToClusterIndex( iBspMesh ).size();

            // initialize the list of Cell cluster groups
            mCellClusterGroups( tDMI ).resize( tNumSPGs );

            // go over Cluster groups, which correspond to the SPGs
            for ( uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
            {
                // get number of clusters in group
                uint tNumClustersInGroup = tSpgToClusterIndex( iBspMesh )( iSPG ).size();

                // collect a list of Clusters in the current group / SPG
                Vector< std::shared_ptr< mtk::Cluster > > tClustersInGroup( tNumClustersInGroup );
                for ( uint iCluster = 0; iCluster < tNumClustersInGroup; iCluster++ )
                {
                    // get the cluster index
                    moris_index tClusterIndex = tSpgToClusterIndex( iBspMesh )( iSPG )( iCluster );

                    // copy pointer to cluster to the list of clusters in the current group
                    tClustersInGroup( iCluster ) = mCellClusters( tClusterIndex );
                }

                // create and commit a new Cluster group to the list
                mCellClusterGroups( tDMI )( iSPG ) = std::make_shared< xtk::Cell_Cluster_Group >( tDMI, tClustersInGroup );

                // assign the cluster group created to all cluster which it was created from
                for ( uint iCluster = 0; iCluster < tNumClustersInGroup; iCluster++ )
                {
                    tClustersInGroup( iCluster )->set_cluster_group( tDMI, mCellClusterGroups( tDMI )( iSPG ) );
                }
            }
        }    // end for: each B-spline mesh

        // free unused memory
        mCellClusterGroups.shrink_to_fit();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_side_cluster_groups()
    {
        // determine number of B-spline meshes
        uint tNumBspMeshes = mBsplineMeshIndices.numel();
        uint tMaxDMI       = (uint)mBsplineMeshIndices.max();

        // initialize list of Clusters
        mDblSideClusterGroups.resize( tMaxDMI + 1 );

        // get access to the B-spline mesh information
        Vector< xtk::Bspline_Mesh_Info * > &tBsplineMeshInfos = mCutIgMesh->get_bspline_mesh_info();

        // get the facet connectivity
        std::shared_ptr< xtk::Facet_Based_Connectivity > tFacetConnectivity = mCutIgMesh->get_face_connectivity();

        // find the total number of side clusters for initialize storage later
        uint tNumSideClusters = 0;
        for ( uint iSideSet = 0; iSideSet < mSideSets.size(); iSideSet++ )
        {
            // skip visualization ghost sets
            std::string tSideSetName = mSideSetLabels( iSideSet );
            if ( tSideSetName.find( "Ghost" ) != std::string::npos )
            {
                continue;
            }

            // add up
            tNumSideClusters += mSideSets( iSideSet ).size();
        }

        // establish cluster group measures for every B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the discretization mesh index
            moris_index tDMI = mBsplineMeshIndices( iBspMesh );

            // TODO: this estimate needs another look at it
            // estimate the number of side cluster groups as 1.0 times the number of cell cluster groups
            uint tApproxNumSideClusterGroups = mCellClusterGroups( iBspMesh ).size();

            // reserve memory for side cluster groups
            mDblSideClusterGroups( tDMI ).reserve( tApproxNumSideClusterGroups );

            // get the information for the current B-spline mesh
            xtk::Bspline_Mesh_Info *tBsplineMeshInfo = tBsplineMeshInfos( iBspMesh );

            // get the number of SPGs on the current B-spline mesh
            uint tNumSPGs = tBsplineMeshInfo->get_num_SPGs();

            // initialize list of side clusters attached to every SPG
            Vector< Vector< std::shared_ptr< mtk::Cluster > > > tSideClustersAttachedToSpg( tNumSPGs );
            Vector< Vector< moris_index > >                     tSpgsSideClustersAreConnectingTo( tNumSPGs );
            Vector< Vector< moris_index > >                     tGlobalSideOrdinalsSideClustersAreOn( tNumSPGs );

            // reserve memory
            tSideClustersAttachedToSpg.reserve( tNumSideClusters );
            tSpgsSideClustersAreConnectingTo.reserve( tNumSideClusters );
            tGlobalSideOrdinalsSideClustersAreOn.reserve( tNumSideClusters );

            // collect side clusters from each set
            for ( uint iSideSet = 0; iSideSet < mSideSets.size(); iSideSet++ )
            {
                // skip visualization ghost sets
                std::string tSideSetName = mSideSetLabels( iSideSet );
                if ( tSideSetName.find( "Ghost" ) != std::string::npos )
                {
                    continue;
                }

                // get the number of side clusters in the current side set
                uint tNumSideClustersInSet = mSideSets( iSideSet ).size();

                // go through the side clusters on the current set and sort into the bins
                for ( uint iSideClusterOnSet = 0; iSideClusterOnSet < tNumSideClustersInSet; iSideClusterOnSet++ )
                {
                    // get access to the current side cluster
                    std::shared_ptr< xtk::Side_Cluster > tSideCluster = mSideSets( iSideSet )( iSideClusterOnSet );

                    // get the UIPC's index the side cluster sits in
                    moris_index tUipcIndex = tSideCluster->get_interpolation_cell_index();

                    // get the SPG index the UIPC is in
                    moris_index tSpgIndex = mModel->mEnrichment->get_SPG_on_UIPC( iBspMesh, tUipcIndex );

                    // the UIPC should not contain a basis extension
                    MORIS_ASSERT(
                            tSpgIndex != -1 && tSpgIndex != MORIS_INDEX_MAX,
                            "Enriched_Integration_Mesh::setup_side_cluster_groups() - "
                            "Side cluster attached to unzipped IP cell which only exists for basis extension purposes. This should not happen." );

                    // NOTE: the following few steps try to figure out the neighbor the side cluster connects to

                    // get a representative facet for the side cluster
                    const moris::mtk::Cell *tIgCell      = tSideCluster->get_cells_in_side_cluster()( 0 );
                    moris_index             tCellIndex   = tIgCell->get_index();
                    moris_index             tSideOrdinal = tSideCluster->get_cell_side_ordinals()( 0 );

                    // below, get the number of cells the current facet connects to
                    uint        tNumCellsAttachedToFacet = 0;
                    moris_index tFacetIndex              = -1;

                    // if element is coarse it will not be part of facet connectivity
                    if ( tFacetConnectivity->mCellIndexToCellOrdinal.find( tCellIndex ) == tFacetConnectivity->mCellIndexToCellOrdinal.end() )
                    {
                        tNumCellsAttachedToFacet = 1;
                    }

                    // otherwise get information through the facet connectivity
                    else
                    {
                        // locate the cell in the facet connectivity
                        moris_index tCellIndexInFacetConnectivity = tFacetConnectivity->get_cell_ordinal( tCellIndex );

                        // get the facet index
                        tFacetIndex = tFacetConnectivity->mCellToFacet( tCellIndexInFacetConnectivity )( tSideOrdinal );

                        // check if there's another cell attached to the facet
                        tNumCellsAttachedToFacet = tFacetConnectivity->mFacetToCell( tFacetIndex ).size();
                    }

                    // check that the number of cells reportedly attached to the current facet are as expected
                    MORIS_ERROR( tNumCellsAttachedToFacet == 1 || tNumCellsAttachedToFacet == 2,
                            "Enriched_Integration_Mesh::setup_side_cluster_groups() - "
                            "There are %i cells attached to a facet. This shouldn't happen. Each facet is connected to either one or two facets.",
                            tNumCellsAttachedToFacet );

                    // treat side clusters on outer mesh boundaries and at interfaces differently
                    if ( tNumCellsAttachedToFacet == 1 )    // case: outer mesh boundary
                    {
                        // compute the outward normal
                        Matrix< DDRMat > tNormal = tIgCell->compute_outward_side_normal( tSideOrdinal );

                        // match the outward normal to a global side ordinal direction
                        moris_index tGlobalSideOrdinal = xtk::match_normal_to_side_ordinal( tNormal );

                        // MORIS_ERROR( tGlobalSideOrdinal != MORIS_INDEX_MAX,
                        //         "Enriched_Integration_Mesh::setup_side_cluster_groups() - "
                        //         "Facet of side cluster is only connected to single IG-Cell #%i but not an ordinal of the global mesh block. Normal = %s",
                        //         tCellIndex,
                        //         ios::stringify_log( tNormal ).c_str() );

                        // store away this information
                        tSideClustersAttachedToSpg( tSpgIndex ).push_back( tSideCluster );
                        tSpgsSideClustersAreConnectingTo( tSpgIndex ).push_back( -1 );
                        tGlobalSideOrdinalsSideClustersAreOn( tSpgIndex ).push_back( tGlobalSideOrdinal );
                    }
                    else    // case: interface facet
                    {
                        // get the other IG cell's index
                        moris_index tNeighborCellIndex = tFacetConnectivity->mFacetToCell( tFacetIndex )( 1 )->get_index();
                        if ( tNeighborCellIndex == tCellIndex )
                        {
                            tNeighborCellIndex = tFacetConnectivity->mFacetToCell( tFacetIndex )( 0 )->get_index();
                        }

                        // find the SPG of the neighbor cell
                        moris_index tNeighborSpIndex  = mCutIgMesh->get_ig_cell_subphase_index( tNeighborCellIndex );
                        moris_index tNeighborSpgIndex = tBsplineMeshInfo->mSpToSpgMap( tNeighborSpIndex );

                        // store away this information
                        tSideClustersAttachedToSpg( tSpgIndex ).push_back( tSideCluster );
                        tSpgsSideClustersAreConnectingTo( tSpgIndex ).push_back( tNeighborSpgIndex );
                        tGlobalSideOrdinalsSideClustersAreOn( tSpgIndex ).push_back( -1 );
                    }

                }    // end for: each side cluster in set

            }        // end for: each side set

            // go through SPGs and collect side cluster groups related to each one
            for ( uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
            {
                // the the number of side clusters related to the current SPG
                uint tNumClustersOnSpg = tSideClustersAttachedToSpg( iSPG ).size();

                // find the number of outer and interface side cluster groups
                uint                            tNumSideClusterGroups = 0;
                Vector< moris_index >           tBinsForClusters( tNumClustersOnSpg, MORIS_INDEX_MAX );
                map< moris_index, moris_index > tInterfaceClusterGroups;
                map< moris_index, moris_index > tBoundaryClusterGroups;


                // establish for which side ordinals and for which neighbor SPGs side cluster groups will need to be constructed
                for ( uint iClusterOnSpg = 0; iClusterOnSpg < tNumClustersOnSpg; iClusterOnSpg++ )
                {
                    // case: is boundary cluster
                    if ( tGlobalSideOrdinalsSideClustersAreOn( iSPG )( iClusterOnSpg ) != -1 )
                    {
                        // get the side ordinal
                        moris_index tGlobalSideOrdinal = tGlobalSideOrdinalsSideClustersAreOn( iSPG )( iClusterOnSpg );

                        // check if this side ordinal will already get a group, if not mark it so
                        if ( !tBoundaryClusterGroups.key_exists( tGlobalSideOrdinal ) )
                        {
                            tBoundaryClusterGroups[ tGlobalSideOrdinal ] = tNumSideClusterGroups;
                            tBinsForClusters( iClusterOnSpg )            = tNumSideClusterGroups;
                            tNumSideClusterGroups++;
                        }
                        else
                        {
                            moris_index tBinIndex             = tBoundaryClusterGroups.find( tGlobalSideOrdinal );
                            tBinsForClusters( iClusterOnSpg ) = tBinIndex;
                        }
                    }

                    // case: is interface cluster
                    else if ( tSpgsSideClustersAreConnectingTo( iSPG )( iClusterOnSpg ) != -1 )
                    {
                        // get the neighbor SPG
                        moris_index tNeighborSPG = tSpgsSideClustersAreConnectingTo( iSPG )( iClusterOnSpg );

                        // check if this side ordinal will already get a group, if not mark it so
                        if ( !tInterfaceClusterGroups.key_exists( tNeighborSPG ) )
                        {
                            tInterfaceClusterGroups[ tNeighborSPG ] = tNumSideClusterGroups;
                            tBinsForClusters( iClusterOnSpg )       = tNumSideClusterGroups;
                            tNumSideClusterGroups++;
                        }
                        else
                        {
                            moris_index tBinIndex             = tInterfaceClusterGroups.find( tNeighborSPG );
                            tBinsForClusters( iClusterOnSpg ) = tBinIndex;
                        }
                    }

                    // something went wrong
                    else
                    {
                        MORIS_ERROR( false,
                                "Enriched_Integration_Mesh::setup_side_cluster_groups() - "
                                "Side cluster is neither marked as boundary nor as interface side cluster. "
                                "Something must have gone wrong." );
                    }

                }    // end for: find group for each cluster on the SPG

                // initialize bins to sort the clusters into
                Vector< Vector< std::shared_ptr< mtk::Cluster > > > tClusterGroups( tNumSideClusterGroups );

                // sort each of the side clusters into the bins
                for ( uint iClusterOnSpg = 0; iClusterOnSpg < tNumClustersOnSpg; iClusterOnSpg++ )
                {
                    // get the bin
                    moris_index tBinIndex = tBinsForClusters( iClusterOnSpg );

                    // sort side clusters into cluster groups
                    tClusterGroups( tBinIndex ).push_back( tSideClustersAttachedToSpg( iSPG )( iClusterOnSpg ) );

                }    // end for: sort each cluster on the SPG into the groups

                // create cluster groups from each of the bins and reversely assign that cluster group to all the side clusters in it
                for ( uint iBin = 0; iBin < tNumSideClusterGroups; iBin++ )
                {
                    // get the corresponding bulk cluster group
                    moris_index                           tUipcIndex                  = tClusterGroups( iBin )( 0 )->get_interpolation_cell_index();
                    std::shared_ptr< mtk::Cluster_Group > tAssociatedCellClusterGroup = mCellClusters( tUipcIndex )->get_cluster_group( tDMI );
                    MORIS_ASSERT( tAssociatedCellClusterGroup.get() != nullptr,
                            "Enriched_Integration_Mesh::setup_side_cluster_groups() - "
                            "Cluster corresponding to UIPC does not have a Cluster group associated with it yet." );

                    // create side cluster group
                    mDblSideClusterGroups( tDMI ).push_back(
                            std::make_shared< xtk::Side_Cluster_Group >( tDMI, tClusterGroups( iBin ), tAssociatedCellClusterGroup ) );

                    // index of the newly created Cluster Group in the list
                    uint tNewSideClusterGroupIndex = mDblSideClusterGroups( tDMI ).size() - 1;

                    // assign the cluster group created to all cluster which it was created from
                    for ( uint iCluster = 0; iCluster < tClusterGroups( iBin ).size(); iCluster++ )
                    {
                        tClusterGroups( iBin )( iCluster )->set_cluster_group( tDMI, mDblSideClusterGroups( tDMI )( tNewSideClusterGroupIndex ) );
                    }
                }

            }    // end for: each SPG

            // free unused memory
            mDblSideClusterGroups( tDMI ).shrink_to_fit();

            // log how good the memory reservation works
            MORIS_LOG_INFO(
                    "B-spline Mesh %i: Number of side cluster groups: Estimated: %i | Actual: %zu",
                    tDMI,
                    tApproxNumSideClusterGroups,
                    mDblSideClusterGroups( tDMI ).size() );

        }    // end for: each B-spline mesh

        // free unused memory
        mDblSideClusterGroups.shrink_to_fit();

    }    // end function: Enriched_Integration_Mesh::setup_side_cluster_groups()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_dbl_side_cluster_groups()
    {
        // determine number of B-spline meshes
        uint tNumBspMeshes = mBsplineMeshIndices.numel();
        uint tMaxDMI       = (uint)mBsplineMeshIndices.max();

        // initialize list of Clusters
        mDblSideClusterGroups.resize( tMaxDMI + 1 );

        // get access to the B-spline mesh information
        Vector< xtk::Bspline_Mesh_Info * > &tBsplineMeshInfos = mCutIgMesh->get_bspline_mesh_info();

        // get the facet connectivity
        std::shared_ptr< xtk::Facet_Based_Connectivity > tFacetConnectivity = mCutIgMesh->get_face_connectivity();

        // establish cluster group measures for every B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the discretization mesh index
            moris_index tDMI = mBsplineMeshIndices( iBspMesh );

            // TODO: this estimate needs another look at it
            // estimate the number of dbl side cluster groups as 1.0 times the number of cell cluster groups
            uint tApproxNumDblSideClusterGroups = mCellClusterGroups( iBspMesh ).size();

            // reserve memory for dbl side cluster groups
            mDblSideClusterGroups( tDMI ).reserve( 2 * tApproxNumDblSideClusterGroups );

            // get the information for the current B-spline mesh
            xtk::Bspline_Mesh_Info *tBsplineMeshInfo = tBsplineMeshInfos( iBspMesh );

            // get the number of SPGs on the current B-spline mesh
            uint tNumSPGs = tBsplineMeshInfo->get_num_SPGs();

            // collect side clusters from each set
            for ( uint iDblSideSet = 0; iDblSideSet < mDoubleSideSets.size(); iDblSideSet++ )
            {
                // skip visualization ghost sets
                std::string tDblSideSetName = mDoubleSideSetLabels( iDblSideSet );
                if ( tDblSideSetName.find( "ghost" ) != std::string::npos || tDblSideSetName.find( "Ghost" ) != std::string::npos )
                {
                    continue;
                }

                // get the number of side clusters in the current side set
                uint tNumSideClustersInSet = mDoubleSideSets( iDblSideSet ).size();

                // initialize list of side clusters attached to every SPG in the current set
                // NOTE: for dbl side clusters, the cluster groups will necessarily consist of elements coming from the same dbl side set.
                // NOTE: Hence, we can collect clusters set-wise to reduce memory consumption and swapping
                Vector< Vector< std::shared_ptr< mtk::Cluster > > > tLeaderSideClustersAttachedToSpg( tNumSPGs );
                Vector< Vector< std::shared_ptr< mtk::Cluster > > > tFollowerSideClustersAttachedToSpg( tNumSPGs );
                Vector< Vector< moris_index > >                     tSpgsLeaderSideClustersAreConnectingTo( tNumSPGs );
                Vector< Vector< moris_index > >                     tSpgsFollowerSideClustersAreConnectingTo( tNumSPGs );

                // reserve memory
                tLeaderSideClustersAttachedToSpg.reserve( tNumSideClustersInSet );
                tFollowerSideClustersAttachedToSpg.reserve( tNumSideClustersInSet );
                tSpgsLeaderSideClustersAreConnectingTo.reserve( tNumSideClustersInSet );
                tSpgsFollowerSideClustersAreConnectingTo.reserve( tNumSideClustersInSet );

                // go through the side clusters on the current set and sort into the bins
                for ( uint iDblSideClusterInSet = 0; iDblSideClusterInSet < tNumSideClustersInSet; iDblSideClusterInSet++ )
                {
                    // get access to the current dbl side cluster
                    std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster = mDoubleSideSets( iDblSideSet )( iDblSideClusterInSet );

                    // indices of the associated single sided clusters
                    moris_index tCurrentDblSideLeaderSideClusterIndex   = mDoubleSideSetsLeaderIndex( iDblSideSet )( iDblSideClusterInSet );
                    moris_index tCurrentDblSideFollowerSideClusterIndex = mDoubleSideSetsFollowerIndex( iDblSideSet )( iDblSideClusterInSet );

                    // get pointers to the leader and follower single sided side clusters
                    std::shared_ptr< mtk::Cluster > tLeaderSideCluster   = mDoubleSideSingleSideClusters( tCurrentDblSideLeaderSideClusterIndex );
                    std::shared_ptr< mtk::Cluster > tFollowerSideCluster = mDoubleSideSingleSideClusters( tCurrentDblSideFollowerSideClusterIndex );

                    // get the leader and follower side clusters' UIPC indices
                    moris_index tLeaderUipcIndex   = tLeaderSideCluster->get_interpolation_cell_index();
                    moris_index tFollowerUipcIndex = tFollowerSideCluster->get_interpolation_cell_index();

                    // get the SPG indices
                    moris_index tLeaderSpgIndex   = mModel->mEnrichment->get_SPG_on_UIPC( iBspMesh, tLeaderUipcIndex );
                    moris_index tFollowerSpgIndex = mModel->mEnrichment->get_SPG_on_UIPC( iBspMesh, tFollowerUipcIndex );

                    // the UIPCs should not contain a basis extension
                    MORIS_ASSERT(
                            tLeaderSpgIndex != -1 && tLeaderSpgIndex != MORIS_INDEX_MAX && tFollowerSpgIndex != -1 && tFollowerSpgIndex != MORIS_INDEX_MAX,
                            "Enriched_Integration_Mesh::setup_dbl_side_cluster_groups() - "
                            "Side cluster attached to unzipped IP cell which only exists for basis extension purposes. This should not happen." );

                    // store the clusters based on which SPG they come from and which SPG they connect to
                    tLeaderSideClustersAttachedToSpg( tLeaderSpgIndex ).push_back( tLeaderSideCluster );
                    tFollowerSideClustersAttachedToSpg( tFollowerSpgIndex ).push_back( tFollowerSideCluster );
                    tSpgsLeaderSideClustersAreConnectingTo( tLeaderSpgIndex ).push_back( tFollowerSpgIndex );
                    tSpgsFollowerSideClustersAreConnectingTo( tFollowerSpgIndex ).push_back( tLeaderSpgIndex );

                }    // end for: each dbl side cluster on set

                // go through each SPG and group side clusters attached to them and connecting to the same neighbor SPG
                for ( uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
                {
                    // get the number of leader and follower side clusters connected to this SPG
                    uint tNumLeaderSideClustersOnSPG   = tLeaderSideClustersAttachedToSpg( iSPG ).size();
                    uint tNumFollowerSideClustersOnSPG = tFollowerSideClustersAttachedToSpg( iSPG ).size();

                    // initialize maps that list neighbor SPGs for a given SPG
                    uint                            tNumLeaderSideClusterGroups   = 0;
                    uint                            tNumFollowerSideClusterGroups = 0;
                    Vector< moris_index >           tBinsForLeaderClusters( tNumLeaderSideClustersOnSPG, MORIS_INDEX_MAX );
                    Vector< moris_index >           tBinsForFollowerClusters( tNumFollowerSideClustersOnSPG, MORIS_INDEX_MAX );
                    map< moris_index, moris_index > tLeaderClusterGroupMap;
                    map< moris_index, moris_index > tFollowerClusterGroupMap;

                    // collect leader side cluster groups
                    for ( uint iLeaderSideClusterOnSPG = 0; iLeaderSideClusterOnSPG < tNumLeaderSideClustersOnSPG; iLeaderSideClusterOnSPG++ )
                    {
                        // get the neighbor SPG for the current side cluster
                        moris_index tNeighborSPG = tSpgsLeaderSideClustersAreConnectingTo( iSPG )( iLeaderSideClusterOnSPG );

                        // check if the neighbor SPG has already been found, if not list it
                        if ( !tLeaderClusterGroupMap.key_exists( tNeighborSPG ) )
                        {
                            tLeaderClusterGroupMap[ tNeighborSPG ]            = tNumLeaderSideClusterGroups;
                            tBinsForLeaderClusters( iLeaderSideClusterOnSPG ) = tNumLeaderSideClusterGroups;
                            tNumLeaderSideClusterGroups++;
                        }

                        // if it has been found before just get the correct bin and associated with the cluster
                        else
                        {
                            moris_index tBinIndex                             = tLeaderClusterGroupMap.find( tNeighborSPG );
                            tBinsForLeaderClusters( iLeaderSideClusterOnSPG ) = tBinIndex;
                        }
                    }

                    // collect follower side cluster groups
                    for ( uint iFollowerSideClusterOnSPG = 0; iFollowerSideClusterOnSPG < tNumFollowerSideClustersOnSPG; iFollowerSideClusterOnSPG++ )
                    {
                        // get the neighbor SPG for the current side cluster
                        moris_index tNeighborSPG = tSpgsFollowerSideClustersAreConnectingTo( iSPG )( iFollowerSideClusterOnSPG );

                        // check if the neighbor SPG has already been found, if not list it
                        if ( !tFollowerClusterGroupMap.key_exists( tNeighborSPG ) )
                        {
                            tFollowerClusterGroupMap[ tNeighborSPG ]              = tNumFollowerSideClusterGroups;
                            tBinsForFollowerClusters( iFollowerSideClusterOnSPG ) = tNumFollowerSideClusterGroups;
                            tNumFollowerSideClusterGroups++;
                        }

                        // if it has been found before just get the correct bin and associated with the cluster
                        else
                        {
                            moris_index tBinIndex                                 = tFollowerClusterGroupMap.find( tNeighborSPG );
                            tBinsForFollowerClusters( iFollowerSideClusterOnSPG ) = tBinIndex;
                        }
                    }

                    // initialize bins to sort the clusters into
                    Vector< Vector< std::shared_ptr< mtk::Cluster > > > tLeaderClusterGroups( tNumLeaderSideClusterGroups );
                    Vector< Vector< std::shared_ptr< mtk::Cluster > > > tFollowerClusterGroups( tNumFollowerSideClusterGroups );
                    tLeaderClusterGroups.reserve( tNumLeaderSideClustersOnSPG );
                    tFollowerClusterGroups.reserve( tNumFollowerSideClustersOnSPG );

                    // sort each of the leader side clusters into their respective bins
                    for ( uint iLeaderSideClusterOnSPG = 0; iLeaderSideClusterOnSPG < tNumLeaderSideClustersOnSPG; iLeaderSideClusterOnSPG++ )
                    {
                        // get the bin
                        moris_index tBinIndex = tBinsForLeaderClusters( iLeaderSideClusterOnSPG );

                        // sort side clusters into cluster groups
                        tLeaderClusterGroups( tBinIndex ).push_back( tLeaderSideClustersAttachedToSpg( iSPG )( iLeaderSideClusterOnSPG ) );

                    }    // end for: sort each leader cluster on the SPG into the groups

                    // sort each of the follower side clusters into their respective bins
                    for ( uint iFollowerSideClusterOnSPG = 0; iFollowerSideClusterOnSPG < tNumFollowerSideClustersOnSPG; iFollowerSideClusterOnSPG++ )
                    {
                        // get the bin
                        moris_index tBinIndex = tBinsForFollowerClusters( iFollowerSideClusterOnSPG );

                        // sort side clusters into cluster groups
                        tFollowerClusterGroups( tBinIndex ).push_back( tFollowerSideClustersAttachedToSpg( iSPG )( iFollowerSideClusterOnSPG ) );

                    }    // end for: sort each follower cluster on the SPG into the groups

                    // create cluster groups from each of the leader bins and reversely assign that cluster group to all the side clusters in it
                    for ( uint iLeaderBin = 0; iLeaderBin < tNumLeaderSideClusterGroups; iLeaderBin++ )
                    {
                        // get the corresponding bulk cluster group
                        moris_index                           tUipcIndex                  = tLeaderClusterGroups( iLeaderBin )( 0 )->get_interpolation_cell_index();
                        std::shared_ptr< mtk::Cluster_Group > tAssociatedCellClusterGroup = mCellClusters( tUipcIndex )->get_cluster_group( tDMI );
                        MORIS_ASSERT( tAssociatedCellClusterGroup.get() != nullptr,
                                "Enriched_Integration_Mesh::setup_dbl_side_cluster_groups() - "
                                "Cluster corresponding to UIPC does not have a Cluster group associated with it yet." );

                        // create side cluster group
                        mDblSideClusterGroups( tDMI ).push_back(
                                std::make_shared< xtk::Side_Cluster_Group >( tDMI, tLeaderClusterGroups( iLeaderBin ), tAssociatedCellClusterGroup ) );

                        // index of the newly created Cluster Group in the list
                        uint tNewSideClusterGroupIndex = mDblSideClusterGroups( tDMI ).size() - 1;

                        // assign the cluster group created to all cluster which it was created from
                        for ( uint iCluster = 0; iCluster < tLeaderClusterGroups( iLeaderBin ).size(); iCluster++ )
                        {
                            tLeaderClusterGroups( iLeaderBin )( iCluster )->set_cluster_group( tDMI, mDblSideClusterGroups( tDMI )( tNewSideClusterGroupIndex ) );
                        }
                    }    // end for: each leader cluster group constructed on the current SPG

                    // create cluster groups from each of the follower bins and reversely assign that cluster group to all the side clusters in it
                    for ( uint iFollowerBin = 0; iFollowerBin < tNumFollowerSideClusterGroups; iFollowerBin++ )
                    {
                        // get the corresponding bulk cluster group
                        moris_index                           tUipcIndex                  = tFollowerClusterGroups( iFollowerBin )( 0 )->get_interpolation_cell_index();
                        std::shared_ptr< mtk::Cluster_Group > tAssociatedCellClusterGroup = mCellClusters( tUipcIndex )->get_cluster_group( tDMI );
                        MORIS_ASSERT( tAssociatedCellClusterGroup.get() != nullptr,
                                "Enriched_Integration_Mesh::setup_dbl_side_cluster_groups() - "
                                "Cluster corresponding to UIPC does not have a Cluster group associated with it yet." );

                        // create side cluster group
                        mDblSideClusterGroups( tDMI ).push_back(
                                std::make_shared< xtk::Side_Cluster_Group >( tDMI, tFollowerClusterGroups( iFollowerBin ), tAssociatedCellClusterGroup ) );

                        // index of the newly created Cluster Group in the list
                        uint tNewSideClusterGroupIndex = mDblSideClusterGroups( tDMI ).size() - 1;

                        // assign the cluster group created to all cluster which it was created from
                        for ( uint iCluster = 0; iCluster < tFollowerClusterGroups( iFollowerBin ).size(); iCluster++ )
                        {
                            tFollowerClusterGroups( iFollowerBin )( iCluster )->set_cluster_group( tDMI, mDblSideClusterGroups( tDMI )( tNewSideClusterGroupIndex ) );
                        }
                    }    // end for: each follower cluster group constructed on the current SPG

                }        // end for: each SPG on current B-spline mesh

            }            // end for: each dbl sided side set

            // free unused memory
            mDblSideClusterGroups( tDMI ).shrink_to_fit();

            // log how good the memory reservation works
            MORIS_LOG_INFO(
                    "B-spline Mesh %i: Number of dbl-side cluster groups: Estimated: %i | Actual: %zu",
                    tDMI,
                    tApproxNumDblSideClusterGroups,
                    mDblSideClusterGroups( tDMI ).size() );

        }    // end for: each B-spline mesh

        // free unused memory
        mDblSideClusterGroups.shrink_to_fit();

    }    // end function: Enriched_Integration_Mesh::setup_dbl_side_cluster_groups()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::visualize_cluster_measures()
    {
        //----------------------------------------------------------------
        // compute and visualize cell cluster volumes

        // set a field name for the cluster volumes
        std::string tClusterVolumeFieldName = "ClusterVolume";

        // create a list of Field indices for the Cluster group volume
        moris_index tClusterVolumeFieldIndex = this->create_field( tClusterVolumeFieldName, mtk::EntityRank::ELEMENT );

        // initialize arrays to store fields
        Matrix< DDRMat > tClusterVolumes( 1, this->get_num_elems(), -1.0 );

        // get the number of cell cluster groups
        uint tNumCellClusters = mCellClusters.size();

        // go over Clusters
        for ( uint iCluster = 0; iCluster < tNumCellClusters; iCluster++ )
        {
            // get the cell cluster
            std::shared_ptr< mtk::Cluster > tCluster = mCellClusters( iCluster );

            // compute the cluster volume
            real tClusterVolume = tCluster->compute_cluster_cell_measure();

            // get the cells in cluster
            Vector< mtk::Cell const * > const &tIgCellsInCluster = tCluster->get_primary_cells_in_cluster();

            for ( uint iIgCell = 0; iIgCell < tIgCellsInCluster.size(); iIgCell++ )
            {
                tClusterVolumes( tIgCellsInCluster( iIgCell )->get_index() ) = tClusterVolume;
            }
        }    // end for: each cluster

        // commit field data to exo
        this->add_field_data( tClusterVolumeFieldIndex, mtk::EntityRank::ELEMENT, tClusterVolumes );

        //----------------------------------------------------------------
        // compute and visualize side cluster interface area/length

        // Get side set names
        // Vector< std::string > tSideSetNames = this->get_set_names( this->get_facet_rank() );

        // set a field name for the cluster volumes
        std::string tSideClusterSizeFieldName = "SideClusterSize";

        // get the total number of side clusters and facets therein
        uint tNumSideSets = mSideSets.size();
        for ( uint iSideSet = 0; iSideSet < tNumSideSets; iSideSet++ )
        {
            // create a list of Field indices for the Cluster group volume
            moris_index tSideClusterSizeFieldIndex = this->create_field( tSideClusterSizeFieldName, this->get_facet_rank(), iSideSet );

            // count the ig cells on the side set
            uint tNumIgCellsOnSideSet = 0;
            for ( auto iSideCluster : mSideSets( iSideSet ) )
            {
                tNumIgCellsOnSideSet += iSideCluster->get_cells_in_side_cluster().size();
            }

            // initialize storage for field data
            Matrix< DDRMat > tSideClusterSizes( 1, tNumIgCellsOnSideSet, -1.0 );

            // go over the clusters in the
            uint tSideElemIndexInSet = 0;
            for ( auto iSideCluster : mSideSets( iSideSet ) )
            {
                // get the cluster measure for this cluster
                real tSideClusterSize = iSideCluster->compute_cluster_cell_side_measure();

                // get the number of side elements
                uint tNumSideElemsInCluster = iSideCluster->get_cells_in_side_cluster().size();

                // assign the same cluster value to all elements in cluster
                for ( uint iSideElem = 0; iSideElem < tNumSideElemsInCluster; iSideElem++ )
                {
                    tSideClusterSizes( tSideElemIndexInSet ) = tSideClusterSize;
                    tSideElemIndexInSet++;
                }
            }

            // commit field data for the current side set
            this->add_field_data( tSideClusterSizeFieldIndex, this->get_facet_rank(), tSideClusterSizes, iSideSet );

        }    // end for: each side set

    }        // end function: Enriched_Integration_Mesh::visualize_cluster_measures()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::visualize_cluster_group_measures( const bool aWriteBsplineClusterInfo )
    {
        // determine number of B-spline meshes
        uint tNumBspMeshes = mBsplineMeshIndices.numel();

        //----------------------------------------------------------------
        // compute and visualize cell cluster groups volumes

        // create a list of Field indices for the Cluster group volume
        Vector< moris_index > tClusterGroupVolumeFieldIndices( tNumBspMeshes );

        // initialize list of field names
        std::string           tVolumeFieldName = "BSplineClusterVolume_B";
        Vector< std::string > tClusterGroupVolumeFields( tNumBspMeshes );

        // initialize arrays to store fields
        Matrix< DDRMat >           tDummy( 1, this->get_num_elems(), -1.0 );
        Vector< Matrix< DDRMat > > tClusterGroupVolumes( tNumBspMeshes, tDummy );

        // compute and output cluster group measures for every B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // create list of field names to write to exodus
            tClusterGroupVolumeFields( iBspMesh )       = tVolumeFieldName + std::to_string( iBspMesh );
            tClusterGroupVolumeFieldIndices( iBspMesh ) = this->create_field( tClusterGroupVolumeFields( iBspMesh ), mtk::EntityRank::ELEMENT );

            // get the number of cell cluster groups
            uint tNumCellClusters = mCellClusters.size();

            // go over Cluster groups, which correspond to the SPGs
            for ( uint iCluster = 0; iCluster < tNumCellClusters; iCluster++ )
            {
                // get the cell cluster
                std::shared_ptr< mtk::Cluster > tCluster = mCellClusters( iCluster );

                // get the cells in cluster
                Vector< mtk::Cell const * > const &tIgCellsInCluster = tCluster->get_primary_cells_in_cluster();

                // skip cluster group volume computation if there are no primary IG cells in the cluster
                if ( tIgCellsInCluster.size() == 0 )
                {
                    continue;
                }

                // skip clusters with no cluster group
                real tClusterGroupVolume = std::numeric_limits< real >::quiet_NaN();
                if ( tCluster->has_cluster_group( iBspMesh ) )
                {
                    tClusterGroupVolume = tCluster->compute_cluster_group_cell_measure( iBspMesh );
                }

                for ( uint iIgCell = 0; iIgCell < tIgCellsInCluster.size(); iIgCell++ )
                {
                    moris_index tIgCellIndex                         = tIgCellsInCluster( iIgCell )->get_index();
                    tClusterGroupVolumes( iBspMesh )( tIgCellIndex ) = tClusterGroupVolume;
                }
            }    // end for: each cluster

            // commit field data to exo
            this->add_field_data( tClusterGroupVolumeFieldIndices( iBspMesh ), mtk::EntityRank::ELEMENT, tClusterGroupVolumes( iBspMesh ) );

        }    // end for: each B-spline mesh

        //----------------------------------------------------------------
        // compute and visualize side cluster groups side length/surface

        // Get side set names
        // Vector< std::string > tSideSetNames = this->get_set_names( this->get_facet_rank() );

        // get the names of the fields for every B-spline mesh
        Vector< std::string > tSideClusterGroupSizeFieldNames( tNumBspMeshes );
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            tSideClusterGroupSizeFieldNames( iBspMesh ) = "SideClusterSize_B" + std::to_string( iBspMesh );
        }

        // get the total number of side clusters and facets therein
        uint tNumSideSets = mSideSets.size();
        for ( uint iSideSet = 0; iSideSet < tNumSideSets; iSideSet++ )
        {
            // skip visualization ghost sets
            std::string tSideSetName = mSideSetLabels( iSideSet );
            if ( tSideSetName.find( "Ghost" ) != std::string::npos )
            {
                continue;
            }

            // count the ig cells on the side set
            uint tNumIgCellsOnSideSet = 0;
            for ( auto iSideCluster : mSideSets( iSideSet ) )
            {
                tNumIgCellsOnSideSet += iSideCluster->get_cells_in_side_cluster().size();
            }

            // compute and output cluster group measures for every B-spline mesh
            for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
            {
                // get the field name for the current B-spline mesh
                std::string tSideClusterGroupSizeFieldName = tSideClusterGroupSizeFieldNames( iBspMesh );

                // create a list of Field indices for the Cluster group volume
                moris_index tSideClusterGroupSizeFieldIndex = this->create_field( tSideClusterGroupSizeFieldName, this->get_facet_rank(), iSideSet );

                // initialize storage for field data
                Matrix< DDRMat > tSideClusterGroupSizes( 1, tNumIgCellsOnSideSet, -1.0 );

                // go over the clusters in the
                uint tSideElemIndexInSet = 0;
                for ( auto iSideCluster : mSideSets( iSideSet ) )
                {
                    // get the number of side elements
                    uint tNumSideElemsInCluster = iSideCluster->get_cells_in_side_cluster().size();

                    // if there are no facets in side cluster, skip the computation of facet areas/edge lengths
                    if ( tNumSideElemsInCluster == 0 )
                    {
                        continue;
                    }

                    // get the cluster measure for this cluster
                    real tSideClusterGroupSize = iSideCluster->compute_cluster_group_cell_side_measure( iBspMesh );

                    // assign the same cluster value to all elements in cluster
                    for ( uint iSideElem = 0; iSideElem < tNumSideElemsInCluster; iSideElem++ )
                    {
                        tSideClusterGroupSizes( tSideElemIndexInSet ) = tSideClusterGroupSize;
                        tSideElemIndexInSet++;
                    }
                }

                // commit field data for the current side set
                this->add_field_data( tSideClusterGroupSizeFieldIndex, this->get_facet_rank(), tSideClusterGroupSizes, iSideSet );

            }    // end for: each B-spline mesh

        }        // end for: each side set

        //----------------------------------------------------------------
        // Generate and write SPG fields

        // skip this step if it is not requested
        if ( !aWriteBsplineClusterInfo )
        {
            return;
        }

        // arrays to store the field indices where the SPG IDs and SPG indices fields will be stored for visualization
        Vector< moris::moris_index > tSpgIndexFieldIndices( tNumBspMeshes );
        Vector< moris::moris_index > tSpgIdFieldIndices( tNumBspMeshes );

        // create list with information where the write the output fields
        std::string           tSpgIndexFieldName = "SPG_Indices_B";
        std::string           tSpgIdFieldName    = "SPG_IDs_B";
        Vector< std::string > tSpgIndexFieldNames( tNumBspMeshes );
        Vector< std::string > tSpgIdFieldNames( tNumBspMeshes );
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            tSpgIndexFieldNames( iBspMesh )   = tSpgIndexFieldName + std::to_string( iBspMesh );
            tSpgIdFieldNames( iBspMesh )      = tSpgIdFieldName + std::to_string( iBspMesh );
            tSpgIndexFieldIndices( iBspMesh ) = this->create_field( tSpgIndexFieldNames( iBspMesh ), mtk::EntityRank::ELEMENT );
            tSpgIdFieldIndices( iBspMesh )    = this->create_field( tSpgIdFieldNames( iBspMesh ), mtk::EntityRank::ELEMENT );
        }

        // initialize list that holds the SPG indices for every IG cell for every B-spline mesh
        Vector< Matrix< DDRMat > > tSpgIndices( tNumBspMeshes, tDummy );
        Vector< Matrix< DDRMat > > tSpgIds( tNumBspMeshes, tDummy );

        // retrieve the IG cell's SPG indices for every B-spline mesh
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the B-spline mesh info
            Bspline_Mesh_Info *tBsplineMeshInfo = mModel->mCutIntegrationMesh->get_bspline_mesh_info()( iBspMesh );

            // get the integration cells on every SPG and assign
            for ( uint iSPG = 0; iSPG < tBsplineMeshInfo->get_num_SPGs(); iSPG++ )
            {
                // get the the integration cell belonging to the current SPG
                Vector< moris_index > const &tIgCellsInSPG = tBsplineMeshInfo->mSubphaseGroups( iSPG )->get_ig_cell_indices_in_group();

                // get the SPG's index
                moris_id tSpgId = tBsplineMeshInfo->get_id_for_spg_index( iSPG );

                // assign current SPG index to all IG cells in SPG
                for ( uint iIgCell = 0; iIgCell < tIgCellsInSPG.size(); iIgCell++ )
                {
                    tSpgIndices( iBspMesh )( tIgCellsInSPG( iIgCell ) ) = (real)iSPG;
                    tSpgIds( iBspMesh )( tIgCellsInSPG( iIgCell ) )     = (real)tSpgId;
                }
            }

            // commit SPG field data to exo
            this->add_field_data( tSpgIndexFieldIndices( iBspMesh ), mtk::EntityRank::ELEMENT, tSpgIndices( iBspMesh ) );
            this->add_field_data( tSpgIdFieldIndices( iBspMesh ), mtk::EntityRank::ELEMENT, tSpgIds( iBspMesh ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_double_sided_interface_sides()
    {
        this->declare_interface_double_side_sets();

        this->create_interface_double_side_sets_and_clusters();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::declare_interface_double_side_sets()
    {
        uint tNumBulkPhases = mModel->get_geom_engine()->get_num_bulk_phase();

        Vector< std::string > tDoubleInterfaceSideNames;

        mBulkPhaseToDblSideIndex.resize( tNumBulkPhases, tNumBulkPhases );
        mBulkPhaseToDblSideIndex.fill( MORIS_INDEX_MAX );

        moris_index tCount = 0;

        Vector< Matrix< IndexMat > > tInterfaceLeaderSideColors;
        Vector< Matrix< IndexMat > > tInterfaceFollowerSideColors;

        for ( moris::moris_index iP0 = 0; iP0 < (moris_index)tNumBulkPhases; iP0++ )
        {
            for ( moris::moris_index iP1 = iP0 + 1; iP1 < (moris_index)tNumBulkPhases; iP1++ )
            {

                std::string tInterfaceSideSetName = this->get_dbl_interface_side_set_name( iP0, iP1 );

                tDoubleInterfaceSideNames.push_back( tInterfaceSideSetName );
                tInterfaceLeaderSideColors.push_back( { { iP0 } } );
                tInterfaceFollowerSideColors.push_back( { { iP1 } } );

                mBulkPhaseToDblSideIndex( iP0, iP1 ) = tCount;
                mBulkPhaseToDblSideIndex( iP1, iP0 ) = tCount;
                tCount++;
            }
        }

        Vector< moris_index > tDblSideSetOrds = this->register_double_side_set_names( tDoubleInterfaceSideNames );

        // set interface side set colors
        for ( moris_index iSS = 0; iSS < (moris_index)tDblSideSetOrds.size(); iSS++ )
        {
            this->set_double_side_set_colors( tDblSideSetOrds( iSS ), tInterfaceLeaderSideColors( iSS ), tInterfaceFollowerSideColors( iSS ) );
        }
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::get_dbl_side_set_index(
            moris_index aPhase0,
            moris_index aPhase1 )
    {
        MORIS_ASSERT(
                aPhase0 < aPhase1,
                "Enriched_Integration_Mesh::get_dbl_side_set_index() - "
                "Double side sets are defined from low phase index to high." );

        MORIS_ASSERT(
                mDoubleSideSetLabels( mBulkPhaseToDblSideIndex( aPhase0, aPhase1 ) ) == this->get_dbl_interface_side_set_name( aPhase0, aPhase1 ),
                "Enriched_Integration_Mesh::get_dbl_side_set_index() - "
                "Interface double side set not showing up in correct index" );

        return mBulkPhaseToDblSideIndex( aPhase0, aPhase1 );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_interface_double_side_sets_and_clusters()
    {
        Tracer tTracer( "XTK", "Enriched Integration Mesh", "create_interface_double_side_sets_and_clusters", mModel->mVerboseLevel, 1 );

        // tool for generating double sided interface
        Integration_Mesh_Generator tIGMeshGen;

        // access interpolation mesh
        Enriched_Interpolation_Mesh *tEnrInterpMesh = mModel->mEnrichedInterpMesh( 0 );

        // get the enriched interpolation cell
        Vector< Interpolation_Cell_Unzipped * > &tEnrIpCells = tEnrInterpMesh->get_enriched_interpolation_cells();

        Vector< Vector< std::shared_ptr< IG_Cell_Double_Side_Group > > > const &tDoubleSidedInterface = mCutIgMesh->get_bulk_phase_to_bulk_phase_dbl_side_interface();

        // for a subphase to subphase side cluster - value in map is the location in tSideClusters
        Vector< IndexMap > tSubphaseToSubphaseSideClusterIndex( mCutIgMesh->get_num_subphases() );

        Vector< std::shared_ptr< xtk::Side_Cluster > > tSideClusters;
        Vector< Vector< moris_index > >                tSideClusterSideOrdinals;

        // get access to the map linking SPs to their corresponding cluster/UIPC indices
        Vector< moris_index > const &tSubphaseIndexToEnrIpCellIndex = mModel->mEnrichment->get_subphase_to_UIPC_map();

        for ( uint iBP0 = 0; iBP0 < tDoubleSidedInterface.size(); iBP0++ )
        {
            for ( uint iBP1 = 0; iBP1 < tDoubleSidedInterface.size(); iBP1++ )
            {
                if ( tDoubleSidedInterface( iBP0 )( iBP1 ) != nullptr )
                {
                    // tDoubleSidedInterface(iBP0)(iBP1)->print();

                    // IndexMap tSubphaseIndexToClusterInterfaceOrd;

                    // access pointer
                    std::shared_ptr< IG_Cell_Double_Side_Group > tDblSideGroup = tDoubleSidedInterface( iBP0 )( iBP1 );

                    // iterate through the facet pairs
                    for ( uint iDblFacet = 0; iDblFacet < tDoubleSidedInterface( iBP0 )( iBP1 )->mLeaderIgCells.size(); iDblFacet++ )
                    {
                        moris_index tLeaderSubphaseIndex   = mCutIgMesh->get_ig_cell_subphase_index( tDblSideGroup->mLeaderIgCells( iDblFacet )->get_index() );
                        moris_index tFollowerSubphaseIndex = mCutIgMesh->get_ig_cell_subphase_index( tDblSideGroup->mFollowerIgCells( iDblFacet )->get_index() );

                        moris_index tLeaderUIPCIndex   = tSubphaseIndexToEnrIpCellIndex( tLeaderSubphaseIndex );
                        moris_index tFollowerUIPCIndex = tSubphaseIndexToEnrIpCellIndex( tFollowerSubphaseIndex );

                        moris_index tLeaderCellSideOrdinal   = tDblSideGroup->mLeaderIgCellSideOrdinals( iDblFacet );
                        moris_index tFollowerCellSideOrdinal = tDblSideGroup->mFollowerIgCellSideOrdinals( iDblFacet );

                        moris_index tLeaderBulkPhase   = mCutIgMesh->get_subphase_bulk_phase( tLeaderSubphaseIndex );
                        moris_index tFollowerBulkPhase = mCutIgMesh->get_subphase_bulk_phase( tFollowerSubphaseIndex );

                        moris_index tLeaderOwner = tEnrIpCells( tLeaderUIPCIndex )->get_owner();

                        if ( tLeaderBulkPhase < tFollowerBulkPhase && tLeaderOwner == moris::par_rank() )
                        {

                            auto tLeaderClusterIter = tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).find( tFollowerSubphaseIndex );

                            // do we have a double side cluster between these subphases yet?

                            // if not we construct the two side clusters
                            if ( tLeaderClusterIter == tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).end() )
                            {
                                // index of double side set
                                moris_index tDoubleSideSetIndex = this->get_dbl_side_set_index( tLeaderBulkPhase, tFollowerBulkPhase );

                                // if not construct one from leader to follower
                                moris_index tNewClusterIndex = tSideClusters.size();
                                tSideClusters.push_back( std::make_shared< xtk::Side_Cluster >() );
                                tSideClusterSideOrdinals.push_back( Vector< moris_index >() );

                                // add leader side cluster
                                mDoubleSideSetsLeaderIndex( tDoubleSideSetIndex ).push_back( mDoubleSideSingleSideClusters.size() );
                                mDoubleSideSingleSideClusters.push_back( tSideClusters( tNewClusterIndex ) );

                                tSideClusters( tNewClusterIndex )->mInterpolationCell     = tEnrIpCells( tLeaderUIPCIndex );
                                tSideClusters( tNewClusterIndex )->mTrivial               = false;
                                tSideClusters( tNewClusterIndex )->mAssociatedCellCluster = &this->get_xtk_cell_cluster( *tSideClusters( tNewClusterIndex )->mInterpolationCell );

                                // leader vertex group
                                std::shared_ptr< IG_Vertex_Group > tLeaderVertexGroup =
                                        mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tEnrIpCells( tLeaderUIPCIndex )->get_base_cell()->get_index() ) );

                                tSideClusters( tNewClusterIndex )->set_ig_vertex_group( tLeaderVertexGroup );

                                // add to the map
                                tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex )[ tFollowerSubphaseIndex ] = tNewClusterIndex;

                                // construct one from follower to leader
                                tNewClusterIndex = tSideClusters.size();
                                tSideClusters.push_back( std::make_shared< xtk::Side_Cluster >() );
                                tSideClusterSideOrdinals.push_back( Vector< moris_index >() );

                                // add leader side cluster to double side set
                                mDoubleSideSetsFollowerIndex( tDoubleSideSetIndex ).push_back( mDoubleSideSingleSideClusters.size() );
                                mDoubleSideSingleSideClusters.push_back( tSideClusters( tNewClusterIndex ) );

                                tSideClusters( tNewClusterIndex )->mInterpolationCell     = tEnrIpCells( tFollowerUIPCIndex );
                                tSideClusters( tNewClusterIndex )->mTrivial               = false;
                                tSideClusters( tNewClusterIndex )->mAssociatedCellCluster = &this->get_xtk_cell_cluster( *tSideClusters( tNewClusterIndex )->mInterpolationCell );

                                // follower vertex group
                                std::shared_ptr< IG_Vertex_Group > tFollowerVertexGroup =
                                        mCutIgMesh->get_vertex_group( mCutIgMesh->get_parent_cell_group_index( tEnrIpCells( tFollowerUIPCIndex )->get_base_cell()->get_index() ) );

                                tSideClusters( tNewClusterIndex )->set_ig_vertex_group( tFollowerVertexGroup );

                                MORIS_ASSERT( tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).find( tLeaderSubphaseIndex ) == tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).end(),
                                        "Enriched_Integration_Mesh::create_interface_double_side_sets_and_clusters() - Nonconcurrent leader follower interface construction occurred." );

                                // add to the map
                                tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex )[ tLeaderSubphaseIndex ] = tNewClusterIndex;

                                // create double side set
                                std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster = std::make_shared< mtk::Double_Side_Cluster >(
                                        tSideClusters( tNewClusterIndex - 1 ).get(),
                                        tSideClusters( tNewClusterIndex ).get(),
                                        tSideClusters( tNewClusterIndex - 1 )->get_vertices_in_cluster() );

                                mDoubleSideClusters.push_back( tDblSideCluster );
                                mDoubleSideSets( tDoubleSideSetIndex ).push_back( tDblSideCluster );
                            }

                            // get the relevant side cluster indices
                            moris_index tLeaderToFollowerSideClusterIndex = tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).find( tFollowerSubphaseIndex )->second;
                            moris_index tFollowerToLeaderSideClusterIndex = tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).find( tLeaderSubphaseIndex )->second;

                            MORIS_ASSERT(
                                    tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).find( tFollowerSubphaseIndex ) != tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).end(),
                                    "Enriched_Integration_Mesh::create_interface_double_side_sets_and_clusters() - Side cluster not constructed" );
                            MORIS_ASSERT(
                                    tSubphaseToSubphaseSideClusterIndex( tFollowerSubphaseIndex ).find( tLeaderSubphaseIndex ) != tSubphaseToSubphaseSideClusterIndex( tLeaderSubphaseIndex ).end(),
                                    "Enriched_Integration_Mesh::create_interface_double_side_sets_and_clusters() - Side cluster not constructed" );

                            // place these in a coincident struct which we will convert to matrix and store in the cluster later
                            tSideClusterSideOrdinals( tLeaderToFollowerSideClusterIndex ).push_back( tLeaderCellSideOrdinal );
                            tSideClusterSideOrdinals( tFollowerToLeaderSideClusterIndex ).push_back( tFollowerCellSideOrdinal );

                            // add integration cell pointers
                            tSideClusters( tLeaderToFollowerSideClusterIndex )->mIntegrationCells.push_back( tDblSideGroup->mLeaderIgCells( iDblFacet ) );
                            tSideClusters( tFollowerToLeaderSideClusterIndex )->mIntegrationCells.push_back( tDblSideGroup->mFollowerIgCells( iDblFacet ) );
                        }
                    }
                }
            }
        }

        // convert the cells of side ordinals to matrix and add to clusters
        for ( uint iSC = 0; iSC < tSideClusterSideOrdinals.size(); iSC++ )
        {
            Matrix< IndexMat > tSideOrds( 1, tSideClusterSideOrdinals( iSC ).size() );

            for ( uint iSide = 0; iSide < tSideClusterSideOrdinals( iSC ).size(); iSide++ )
            {
                tSideOrds( iSide ) = tSideClusterSideOrdinals( iSC )( iSide );
            }

            tSideClusters( iSC )->mIntegrationCellSideOrdinals = tSideOrds;
        }

        // build list of side set indices
        Vector< moris_index > tDoubleSideSetIndexList;
        tDoubleSideSetIndexList.reserve( mDoubleSideSets.size() - mListOfDoubleSideSets.size() );

        // construct double side set interfaces
        for ( uint Ik = mListOfDoubleSideSets.size(); Ik < mDoubleSideSets.size(); Ik++ )
        {
            tDoubleSideSetIndexList.push_back( Ik );
        }

        // set double side set indices
        this->commit_double_side_set( tDoubleSideSetIndexList );
        this->communicate_sets_of_type( mtk::SetType::DOUBLE_SIDED_SIDESET );
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::add_side_to_cluster(
            std::shared_ptr< xtk::Side_Cluster > aSideCluster,
            moris_index                          aCellIndex,
            moris_index                          aSideOrdinal )
    {
        // number of current sides in cluster
        uint tNumCurrentSides = aSideCluster->get_num_sides_in_cluster();

        aSideCluster->mIntegrationCells.push_back( &this->get_mtk_cell( aCellIndex ) );

        // add sides
        aSideCluster->mIntegrationCellSideOrdinals.resize( 1, tNumCurrentSides + 1 );
        aSideCluster->mIntegrationCellSideOrdinals( tNumCurrentSides ) = aSideOrdinal;
    }

    //------------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Integration_Mesh::split_set_name_by_bulk_phase( std::string aBaseName )
    {
        uint                  tNumPhases = mModel->mGeometryEngine->get_num_bulk_phase();
        Vector< std::string > tSetNames( tNumPhases );
        for ( uint i = 0; i < tNumPhases; i++ )
        {
            tSetNames( i ) = aBaseName + "_p" + std::to_string( i );
        }

        return tSetNames;
    }

    //------------------------------------------------------------------------------

    Vector< std::string >
    Enriched_Integration_Mesh::split_set_name_by_child_no_child( std::string aBaseName )
    {
        Vector< std::string > tSetNames( 2 );
        tSetNames( 0 ) = aBaseName + "_c";
        tSetNames( 1 ) = aBaseName + "_n";
        return tSetNames;
    }

    //------------------------------------------------------------------------------
    Vector< moris_index >
    Enriched_Integration_Mesh::register_vertex_set_names( Vector< std::string > const &aVertexSetNames )
    {
        uint tNumSetsToRegister = aVertexSetNames.size();

        // block set ordinals
        Vector< moris_index > tVertexSetOrds( tNumSetsToRegister );

        // iterate and add sets
        for ( uint i = 0; i < tNumSetsToRegister; i++ )
        {
            tVertexSetOrds( i ) = mVertexSetNames.size();

            mVertexSetNames.push_back( aVertexSetNames( i ) );
            MORIS_ASSERT( mVertexSetLabelToOrd.find( aVertexSetNames( i ) ) == mVertexSetLabelToOrd.end(),
                    "Duplicate vertex set in mesh" );

            mVertexSetLabelToOrd[ aVertexSetNames( i ) ] = tVertexSetOrds( i );
        }

        mVerticesInVertexSet.resize( mVerticesInVertexSet.size() + tNumSetsToRegister );
        mVertexSetColors.resize( mVertexSetColors.size() + tNumSetsToRegister );

        return tVertexSetOrds;
    }

    //------------------------------------------------------------------------------

    Vector< moris_index >
    Enriched_Integration_Mesh::register_block_set_names_with_cell_topo(
            Vector< std::string > const &aBlockSetNames,
            mtk::CellTopology            aBlockTopology )
    {
        uint tNumSetsToRegister = aBlockSetNames.size();

        // block set ordinals
        Vector< moris_index > tBlockSetOrds( tNumSetsToRegister );

        // iterate and add sets
        for ( uint i = 0; i < tNumSetsToRegister; i++ )
        {
            tBlockSetOrds( i ) = mBlockSetNames.size();

            mBlockSetNames.push_back( aBlockSetNames( i ) );
            mBlockSetTopology.push_back( aBlockTopology );
            MORIS_ASSERT( mBlockSetLabelToOrd.find( aBlockSetNames( i ) ) == mBlockSetLabelToOrd.end(), "Duplicate block set in mesh" );
            mBlockSetLabelToOrd[ aBlockSetNames( i ) ] = tBlockSetOrds( i );
        }

        mPrimaryBlockSetClusters.resize( mPrimaryBlockSetClusters.size() + tNumSetsToRegister );
        mBlockSetColors.resize( mBlockSetColors.size() + tNumSetsToRegister );
        return tBlockSetOrds;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::set_block_set_colors(
            moris_index const        &aBlockSetIndex,
            Matrix< IndexMat > const &aBlockSetColors )
    {
        MORIS_ASSERT( moris::isempty( mBlockSetColors( aBlockSetIndex ) ), "Attempting to overwrite colors of a block set" );

        mBlockSetColors( aBlockSetIndex ) = aBlockSetColors;
    }

    //------------------------------------------------------------------------------

    Vector< moris_index >
    Enriched_Integration_Mesh::register_side_set_names( Vector< std::string > const &aSideSetNames )
    {
        uint tNumSetsToRegister = aSideSetNames.size();

        // block set ordinals
        Vector< moris_index > tSideSetOrds( tNumSetsToRegister );

        // iterate and add sets
        for ( uint i = 0; i < tNumSetsToRegister; i++ )
        {
            tSideSetOrds( i ) = mSideSetLabels.size();

            mSideSetLabels.push_back( aSideSetNames( i ) );
            MORIS_ASSERT( mSideSideSetLabelToOrd.find( aSideSetNames( i ) ) == mSideSideSetLabelToOrd.end(),
                    "Duplicate block set in mesh" );

            mSideSideSetLabelToOrd[ aSideSetNames( i ) ] = tSideSetOrds( i );
        }

        mSideSets.resize( mSideSets.size() + tNumSetsToRegister );
        mSideSetColors.resize( mSideSetColors.size() + tNumSetsToRegister );

        return tSideSetOrds;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::set_side_set_colors(
            moris_index const        &aSideSetIndex,
            Matrix< IndexMat > const &aSideSetColors )
    {

        MORIS_ASSERT( moris::isempty( mSideSetColors( aSideSetIndex ) ),
                "Attempting to overwrite colors of a side set" );

        mSideSetColors( aSideSetIndex ) = aSideSetColors;
    }

    //------------------------------------------------------------------------------

    Vector< moris_index >
    Enriched_Integration_Mesh::register_double_side_set_names( Vector< std::string > const &aDblSideSetNames )
    {
        uint tNumSetsToRegister = aDblSideSetNames.size();

        // block set ordinals
        Vector< moris_index > tDblSideSetOrds( tNumSetsToRegister );

        // get number of already existing double sided side sets (i.e. interface side sets)
        uint tNumInterfaceSets = mDoubleSideSets.size();

        // iterate and add double side sets
        for ( uint iDSS = 0; iDSS < tNumSetsToRegister; iDSS++ )
        {
            // build map associating each ghost side set with a double side set index
            tDblSideSetOrds( iDSS ) = tNumInterfaceSets + iDSS;

            mDoubleSideSetLabels.push_back( aDblSideSetNames( iDSS ) );
            MORIS_ASSERT( mDoubleSideSetLabelToOrd.find( aDblSideSetNames( iDSS ) ) == mDoubleSideSetLabelToOrd.end(), "Duplicate double side set in mesh" );
            mDoubleSideSetLabelToOrd[ aDblSideSetNames( iDSS ) ] = tDblSideSetOrds( iDSS );
        }

        // update dbl. side set info in enriched integration mesh
        mDoubleSideSets.resize( mDoubleSideSets.size() + tNumSetsToRegister );
        mDoubleSideSetsLeaderIndex.resize( mDoubleSideSetsLeaderIndex.size() + tNumSetsToRegister );
        mDoubleSideSetsFollowerIndex.resize( mDoubleSideSetsFollowerIndex.size() + tNumSetsToRegister );
        mLeaderDoubleSideSetColor.resize( mLeaderDoubleSideSetColor.size() + tNumSetsToRegister );
        mFollowerDoubleSideSetColor.resize( mFollowerDoubleSideSetColor.size() + tNumSetsToRegister );

        // return list containing the dbl. side set indices for new ghost dbl. side sets
        return tDblSideSetOrds;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::set_double_side_set_colors(
            moris_index const        &aDblSideSetIndex,
            Matrix< IndexMat > const &aLeaderSideColors,
            Matrix< IndexMat > const &aFollowerSideColors )
    {
        MORIS_ASSERT( moris::isempty( mLeaderDoubleSideSetColor( aDblSideSetIndex ) ),
                "Attempting to overwrite colors of a leader side of double side set" );

        MORIS_ASSERT( moris::isempty( mFollowerDoubleSideSetColor( aDblSideSetIndex ) ),
                "Attempting to overwrite colors of a follower side of double side set" );

        mLeaderDoubleSideSetColor( aDblSideSetIndex )   = aLeaderSideColors;
        mFollowerDoubleSideSetColor( aDblSideSetIndex ) = aFollowerSideColors;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::setup_interface_side_sets()
    {
        this->declare_interface_side_sets();

        this->create_interface_side_sets_and_clusters();
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::declare_interface_side_sets()
    {
        uint tNumBulkPhases = mModel->get_geom_engine()->get_num_bulk_phase();

        Vector< std::string >        tInterfaceSideNames;
        Vector< Matrix< IndexMat > > tInterfaceSideColors;
        for ( moris::moris_index iP0 = 0; iP0 < (moris_index)tNumBulkPhases; iP0++ )
        {
            for ( moris::moris_index iP1 = 0; iP1 < (moris_index)tNumBulkPhases; iP1++ )
            {
                if ( iP1 != iP0 )
                {
                    std::string tInterfaceSideSetName = get_interface_side_set_name( 0, iP0, iP1 );

                    tInterfaceSideNames.push_back( tInterfaceSideSetName );
                    tInterfaceSideColors.push_back( { { iP0 } } );
                }
            }
        }

        Vector< moris_index > tSideSetOrds = this->register_side_set_names( tInterfaceSideNames );

        // set interface side set colors
        for ( moris_index iSS = 0; iSS < (moris_index)tSideSetOrds.size(); iSS++ )
        {
            this->set_side_set_colors( tSideSetOrds( iSS ), tInterfaceSideColors( iSS ) );
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_interface_side_sets_and_clusters()
    {
        uint tNumBulkPhases = mModel->get_geom_engine()->get_num_bulk_phase();

        // iterate through bulk phases
        for ( moris::moris_index iBP0 = 0; iBP0 < (moris_index)tNumBulkPhases; iBP0++ )
        {
            // iterate through bulk phase +1 (high to low
            for ( moris::moris_index iBP1 = iBP0 + 1; iBP1 < (moris_index)tNumBulkPhases; iBP1++ )
            {
                // create the interface side sets
                this->create_interface_side_sets_from_interface_double_side_set( iBP0, iBP1 );
            }
        }

        // build list of side set indices
        Vector< moris_index > tSideSetIndexList;
        tSideSetIndexList.reserve( mSideSets.size() - mListOfSideSets.size() );

        for ( uint Ik = mListOfSideSets.size(); Ik < mSideSets.size(); Ik++ )
        {
            tSideSetIndexList.push_back( Ik );
        }

        // commit and communicate the side sets
        this->commit_side_set( tSideSetIndexList );
        this->communicate_sets_of_type( mtk::SetType::SIDESET );

    }    // end function: Enriched_Integration_Mesh::create_interface_side_sets_and_clusters()

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::construct_color_to_set_relationship(
            Vector< Matrix< IndexMat > > const &aSetColors,
            Vector< Vector< moris_index > >    &aColorToSetIndex )
    {
        moris_index tMaxColor = 0;
        for ( uint i = 0; i < aSetColors.size(); i++ )
        {
            moris_index tLocMax = aSetColors( i ).max();
            if ( tLocMax > tMaxColor )
            {
                tMaxColor = tLocMax;
            }
        }

        // size
        aColorToSetIndex.clear();
        aColorToSetIndex.resize( tMaxColor + 1 );

        for ( uint i = 0; i < aSetColors.size(); i++ )
        {
            for ( uint iC = 0; iC < aSetColors( i ).numel(); iC++ )
            {
                aColorToSetIndex( aSetColors( i )( iC ) ).push_back( (moris_index)i );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::create_interface_side_sets_from_interface_double_side_set(
            moris_index const &aBulkphase0,
            moris_index const &aBulkphase1 )
    {
        // get the side set names one going from 0->1 and another 1->0
        std::string tISideNameLeaderToFollower = this->get_interface_side_set_name( 0, aBulkphase0, aBulkphase1 );
        std::string tISideNameFollowerToLeader = this->get_interface_side_set_name( 0, aBulkphase1, aBulkphase0 );

        // get the corresponding indices
        moris_index tISideIndexLeaderToFollower = this->get_side_set_index( tISideNameLeaderToFollower );
        moris_index tISideIndexFollowerToLeader = this->get_side_set_index( tISideNameFollowerToLeader );

        // get the double side set index
        moris_index tDblSideSetIndex = this->get_dbl_side_set_index( aBulkphase0, aBulkphase1 );

        Vector< std::shared_ptr< mtk::Double_Side_Cluster > > &tDblSideClusters = mDoubleSideSets( tDblSideSetIndex );

        // place the clusters in the two side sets
        for ( uint i = 0; i < tDblSideClusters.size(); i++ )
        {
            // get the index
            moris_index tLeaderIndex   = mDoubleSideSetsLeaderIndex( tDblSideSetIndex )( i );
            moris_index tFollowerIndex = mDoubleSideSetsFollowerIndex( tDblSideSetIndex )( i );

            mSideSets( tISideIndexLeaderToFollower ).push_back( mDoubleSideSingleSideClusters( tLeaderIndex ) );
            mSideSets( tISideIndexFollowerToLeader ).push_back( mDoubleSideSingleSideClusters( tFollowerIndex ) );
        }
    }

    //------------------------------------------------------------------------------

    Vector< moris_index >
    Enriched_Integration_Mesh::declare_interface_vertex_sets()
    {
        // number of geometries in the mesh
        uint tNumGeometries = mModel->get_geom_engine()->get_number_of_geometries();

        // allocate a cell of strings
        Vector< std::string > tInterfaceVertexSetNames( tNumGeometries );

        // base set name (interface vertex geometry #)
        std::string tSetNameBase = "iv_g_";

        for ( uint i = 0; i < tNumGeometries; i++ )
        {
            // add the vertex set to the cell
            tInterfaceVertexSetNames( i ) = std::string( tSetNameBase + std::to_string( i ) );
        }

        // register vertex sets
        Vector< moris_index > tVertexSetOrds = this->register_vertex_set_names( tInterfaceVertexSetNames );

        // make the geometric index the color
        for ( uint i = 0; i < tNumGeometries; i++ )
        {
            // the color of the interface node sets is the geometric index
            this->set_vertex_set_color( tVertexSetOrds( i ), Matrix< IndexMat >( { { (moris_index)i } } ) );
        }

        return tVertexSetOrds;
    }

    //------------------------------------------------------------------------------

    void
    Enriched_Integration_Mesh::set_vertex_set_color(
            moris_index const        &aVertexSetIndex,
            Matrix< IndexMat > const &aVertexSetColors )
    {
        MORIS_ASSERT( moris::isempty( mVertexSetColors( aVertexSetIndex ) ),
                "Attempting to overwrite colors of a side set" );

        mVertexSetColors( aVertexSetIndex ) = aVertexSetColors;
    }

    //------------------------------------------------------------------------------

    uint
    Enriched_Integration_Mesh::get_num_elements_in_side_set( const std::string &aSideSetName )
    {
        // find the side set with the requested name
        auto tIter = mSideSideSetLabelToOrd.find( aSideSetName );
        MORIS_ERROR( tIter != mSideSideSetLabelToOrd.end(),
                "XTK::Enriched_Integration_Mesh::get_num_elements_in_side_set() - "
                "Side set with name '%s' does not exist in mesh.",
                aSideSetName.c_str() );
        uint tSideSetOrdinal = tIter->second;

        // count up number of elements in side clusters
        uint tNumElemsInSideSet = 0;
        for ( auto iCluster : mSideSets( tSideSetOrdinal ) )
        {
            tNumElemsInSideSet += iCluster->get_num_primary_cells();
        }

        // return this information
        return tNumElemsInSideSet;
    }

    //------------------------------------------------------------------------------

    bool
    Enriched_Integration_Mesh::field_exists(
            const std::string        aLabel,
            const mtk::EntityRank    aEntityRank,
            const moris::moris_index aSetOrdinal )
    {
        moris::moris_index tIndex = this->get_entity_rank_field_index( aEntityRank );

        // if global list
        if ( aSetOrdinal == MORIS_INDEX_MAX )
        {
            std::unordered_map< std::string, moris_index > const &tGlobalLabelMap = mGlobalSetFieldLabelToIndex( tIndex );
            return tGlobalLabelMap.find( aLabel ) != tGlobalLabelMap.end();
        }
        else
        {
            if ( (moris_index)mSetWiseFieldLabelToIndex( tIndex ).size() <= aSetOrdinal )
            {
                return false;
            }
            else
            {
                std::unordered_map< std::string, moris_index > const &tSetLabelMap = mSetWiseFieldLabelToIndex( tIndex )( aSetOrdinal );
                return tSetLabelMap.find( aLabel ) != tSetLabelMap.end();
            }
        }
    }

    //------------------------------------------------------------------------------

    moris_index
    Enriched_Integration_Mesh::get_entity_rank_field_index( const mtk::EntityRank aEntityRank )
    {
        // MORIS_ERROR( aEntityRank == mtk::EntityRank::NODE || aEntityRank == mtk::EntityRank::ELEMENT,
        //         "Enriched_Integration_Mesh::get_entity_rank_field_index() - Only node and cell fields are supported" );

        moris_index tIndex = MORIS_INDEX_MAX;

        if ( aEntityRank == mtk::EntityRank::NODE )
        {
            tIndex = 0;
        }

        else if ( aEntityRank == mtk::EntityRank::ELEMENT )
        {
            tIndex = 1;
        }

        else if ( aEntityRank == mtk::EntityRank::FACE || aEntityRank == mtk::EntityRank::EDGE )
        {
            tIndex = 2;
        }

        return tIndex;
    }
}    // namespace moris::xtk
