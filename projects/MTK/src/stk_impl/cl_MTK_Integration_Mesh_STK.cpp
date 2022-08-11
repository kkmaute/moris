/*
 * cl_MTK_Integration_Mesh_STK.cpp
 *
 *  Created on: Apr 29, 2019
 *      Author: doble
 */

#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Cell_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"

#include "cl_MTK_Block.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Cell_Info.hpp"
namespace moris
{
    namespace mtk
    {
        // ----------------------------------------------------------------------------

        Integration_Mesh_STK::Integration_Mesh_STK( std::shared_ptr< Mesh_Data_STK > aSTKMeshData )
                : Mesh_Core_STK( aSTKMeshData )
        {
        }

        // ----------------------------------------------------------------------------

        Integration_Mesh_STK::Integration_Mesh_STK(
                std::string  aFileName,
                MtkMeshData *aSuppMeshData,
                const bool   aCreateFacesAndEdges )
                : Mesh_Core_STK( aFileName, aSuppMeshData, aCreateFacesAndEdges )
        {
        }

        // ----------------------------------------------------------------------------

        Integration_Mesh_STK::Integration_Mesh_STK( MtkMeshData &aMeshData )
                : Mesh_Core_STK( aMeshData )
        {
            this->setup_cell_clusters();

            this->setup_blockset_with_cell_clusters_trivial();

            //    this->setup_side_set_clusters_trivial();

            this->collect_all_sets();
        }

        // ----------------------------------------------------------------------------

        Integration_Mesh_STK::Integration_Mesh_STK(
                MtkMeshData        &aMeshData,
                Interpolation_Mesh *aInterpMesh )
                : Mesh_Core_STK( aMeshData )
        {
            // setup cells and cell clusters
            this->setup_cell_clusters( *aInterpMesh, aMeshData.CellClusterInput );
            this->setup_blockset_with_cell_clusters();

            // setup side set clusters
            this->setup_side_set_clusters( *aInterpMesh, aMeshData.SideClusterInput );

            if ( aMeshData.DoubleSideClusterInput != nullptr )
            {
                this->setup_double_side_set_clusters( *aInterpMesh, aMeshData.DoubleSideClusterInput );
            }

            this->collect_all_sets();
        }

        // ----------------------------------------------------------------------------

        Integration_Mesh_STK::Integration_Mesh_STK(
                Interpolation_Mesh &aInterpMesh,
                Cell_Cluster_Input *aCellClusterInput )
        {
            MORIS_ERROR( aInterpMesh.get_mesh_type() == MeshType::STK,
                    "operator= between an interpolation and integration mesh only valid between stk meshes" );

            Interpolation_Mesh_STK *tInterpolationSTK =
                    dynamic_cast< Interpolation_Mesh_STK * >( &aInterpMesh );

            // get the shared data from the stk interpolation mesh
            mSTKMeshData = tInterpolationSTK->get_stk_data_shared_pointer();

            this->setup_cell_clusters( aInterpMesh, aCellClusterInput );

            this->setup_blockset_with_cell_clusters();

            // setup side set clusters
            this->setup_side_set_clusters( aInterpMesh, nullptr );

            this->collect_all_sets();
        }

        // ----------------------------------------------------------------------------

        Integration_Mesh_STK::~Integration_Mesh_STK()
        {

            for ( auto p : mListofBlocks )
            {
                delete p;
            }
            mListofBlocks.clear();

            for ( auto p : mListofSideSets )
            {
                delete p;
            }
            mListofSideSets.clear();

            for ( auto p : mListofDoubleSideSets )
            {
                delete p;
            }
            mListofDoubleSideSets.clear();
        }

        // ----------------------------------------------------------------------------

        Cell_Cluster const &
        Integration_Mesh_STK::get_cell_cluster( Cell const &aInterpCell ) const
        {
            return mCellClusters( aInterpCell.get_index() );
        }

        // ----------------------------------------------------------------------------

        Cell_Cluster const &
        Integration_Mesh_STK::get_cell_cluster( moris_index aInterpCellIndex ) const
        {
            MORIS_ASSERT( aInterpCellIndex < (moris_index)mCellClusters.size(),
                    "Interpolation Cell index out of bounds" );

            return mCellClusters( aInterpCellIndex );
        }

        moris::Cell< std::string >
        Integration_Mesh_STK::get_block_set_names() const
        {
            return mPrimaryBlockSetNames;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< Cluster const * >
        Integration_Mesh_STK::get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const
        {
            MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mPrimaryBlockSetNames.size(),
                    "Requested block set ordinal out of bounds." );

            moris::Cell< moris::moris_index > const &tClusterIndsInSet =
                    mPrimaryBlockSetClusters( aBlockSetOrdinal );

            moris::Cell< Cluster const * > tClusterInSet( tClusterIndsInSet.size() );

            for ( moris::uint i = 0; i < tClusterIndsInSet.size(); i++ )
            {
                tClusterInSet( i ) = &this->get_cell_cluster( tClusterIndsInSet( i ) );
            }

            return tClusterInSet;
        }
        std::string
        Integration_Mesh_STK::get_block_set_label( moris_index aBlockSetOrdinal ) const
        {
            MORIS_ASSERT( aBlockSetOrdinal < (moris_index)mListofBlocks.size(),
                    "Block set ordinal out of bounds" );

            return mPrimaryBlockSetNames( aBlockSetOrdinal );
        }
        moris_index
        Integration_Mesh_STK::get_block_set_index( std::string aBlockSetLabel ) const
        {
            auto tIter = mBlockSetLabelToOrd.find( aBlockSetLabel );

            MORIS_ERROR( tIter != mBlockSetLabelToOrd.end(), "block side set label not found" );

            return tIter->second;
        }
        // ----------------------------------------------------------------------------

        moris::Cell< Cluster const * >
        Integration_Mesh_STK::get_side_set_cluster( moris_index aSideSetOrdinal ) const
        {
            MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSets.size(),
                    "Side set ordinal out of bounds" );

            moris::uint tNumSideClustersInSet = mSideSets( aSideSetOrdinal ).size();

            moris::Cell< Cluster const * > tSideClustersInSet( tNumSideClustersInSet );

            for ( moris::uint i = 0; i < tNumSideClustersInSet; i++ )
            {
                tSideClustersInSet( i ) = &mSideSets( aSideSetOrdinal )( i );
            }

            return tSideClustersInSet;
        }

        // ----------------------------------------------------------------------------

        uint
        Integration_Mesh_STK::get_num_side_sets() const
        {
            return mSideSets.size();
        }

        // ----------------------------------------------------------------------------

        std::string
        Integration_Mesh_STK::get_side_set_label( moris_index aSideSetOrdinal ) const
        {
            MORIS_ASSERT( aSideSetOrdinal < (moris_index)mSideSetLabels.size(),
                    "Side set ordinal out of bounds" );

            return mSideSetLabels( aSideSetOrdinal );
        }

        // ----------------------------------------------------------------------------

        moris_index
        Integration_Mesh_STK::get_side_set_index( std::string aSideSetLabel ) const
        {
            auto tIter = mSideSideSetLabelToOrd.find( aSideSetLabel );

            MORIS_ERROR( tIter != mSideSideSetLabelToOrd.end(),
                    "side side set label not found" );

            return tIter->second;
        }

        // ----------------------------------------------------------------------------

        uint
        Integration_Mesh_STK::get_num_double_sided_sets() const
        {
            return mDoubleSideSetLabels.size();
        }

        // ----------------------------------------------------------------------------

        std::string
        Integration_Mesh_STK::get_double_sided_set_label( moris_index aSideSetOrdinal ) const
        {
            MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(),
                    "Double side set ordinal out of bounds" );

            return mDoubleSideSetLabels( aSideSetOrdinal );
        }

        // ----------------------------------------------------------------------------

        moris_index
        Integration_Mesh_STK::get_double_sided_set_index( std::string aDoubleSideSetLabel ) const
        {
            auto tIter = mDoubleSideSetLabelToOrd.find( aDoubleSideSetLabel );

            MORIS_ERROR( tIter != mDoubleSideSetLabelToOrd.end(),
                    "double side set label not found" );

            return tIter->second;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< Cluster const * >
        Integration_Mesh_STK::get_double_side_set_cluster( moris_index aSideSetOrdinal ) const
        {
            MORIS_ASSERT( aSideSetOrdinal < (moris_index)mDoubleSideSetLabels.size(),
                    "Double side set ordinal out of bounds" );

            return mDoubleSideSets( aSideSetOrdinal );
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh_STK::setup_cell_clusters()
        {
            moris::uint tNumClusters = mSTKMeshData->mMtkCells.size();

            mCellClusters.resize( tNumClusters );

            // iterate through cells
            for ( moris::uint Ik = 0; Ik < tNumClusters; Ik++ )
            {
                moris::mtk::Cell *tCell = &mSTKMeshData->mMtkCells( Ik );

                mCellClusters( Ik ).set_interpolation_cell( tCell );

                mCellClusters( Ik ).add_primary_integration_cell( tCell );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh_STK::setup_blockset_with_cell_clusters_trivial()
        {
            moris::Cell< std::string > mPrimaryBlockSetNames = this->get_set_names( EntityRank::ELEMENT );

            // iterate and create map
            mListofBlocks.resize( mPrimaryBlockSetNames.size(), nullptr );

            for ( moris::uint Ik = 0; Ik < mPrimaryBlockSetNames.size(); Ik++ )
            {
                Matrix< IndexMat > tCellIndices =
                        this->get_set_entity_loc_inds( EntityRank::ELEMENT, mPrimaryBlockSetNames( Ik ) );

                moris::Cell< Cluster const * > tCellClusters( tCellIndices.numel() );

                for ( moris::uint Ii = 0; Ii < tCellIndices.numel(); Ii++ )
                {
                    tCellClusters( Ii ) = &mCellClusters( tCellIndices( Ii ) );
                }

                mListofBlocks( Ik ) = new moris::mtk::Block(
                        mPrimaryBlockSetNames( Ik ),
                        tCellClusters,
                        { { 0 } },
                        this->get_spatial_dim() );

                MORIS_ASSERT( mBlockSetLabelToOrd.find( mPrimaryBlockSetNames( Ik ) ) == mBlockSetLabelToOrd.end(),
                        "Duplicate block set in mesh" );

                mBlockSetLabelToOrd[ mPrimaryBlockSetNames( Ik ) ] = Ik;
            }
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh_STK::setup_side_set_clusters_trivial()
        {
            enum EntityRank tSideSetRank = this->get_facet_rank();

            moris::Cell< std::string > aSideSetNames = this->get_set_names( tSideSetRank );

            mSideSets.resize( aSideSetNames.size() );

            // copy strings labels
            mSideSetLabels.append( aSideSetNames );

            // add to map
            for ( moris::uint i = 0; i < aSideSetNames.size(); i++ )
            {
                MORIS_ASSERT( mSideSideSetLabelToOrd.find( mSideSetLabels( i ) ) == mSideSideSetLabelToOrd.end(),
                        "Duplicate side set label detected." );

                mSideSideSetLabelToOrd[ mSideSetLabels( i ) ] = i;
            }

            // iterate through block sets
            for ( moris::uint i = 0; i < aSideSetNames.size(); i++ )
            {
                // get the cells and side ordinals from the mesh for this side set
                moris::Cell< mtk::Cell const * > tCellsInSet( 0 );
                moris::Matrix< moris::IndexMat > tSideOrdsInSet( 0, 0 );
                this->get_sideset_cells_and_ords( aSideSetNames( i ), tCellsInSet, tSideOrdsInSet );

                // loop over cells in the side set and make sure they have all been included
                for ( moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++ )
                {
                    mSideSets( i ).push_back(
                            Side_Cluster_STK(
                                    tCellsInSet( iIGCell ),
                                    tCellsInSet( iIGCell ),
                                    tCellsInSet( iIGCell )->get_vertices_on_side_ordinal( tSideOrdsInSet( iIGCell ) ),
                                    tSideOrdsInSet( iIGCell ) ) );
                }
            }

            mListofSideSets.resize( mSideSets.size(), nullptr );

            for ( moris::uint Ik = 0; Ik < mListofSideSets.size(); Ik++ )
            {
                mListofSideSets( Ik ) = new moris::mtk::Side_Set(
                        aSideSetNames( Ik ),
                        this->get_side_set_cluster( Ik ),
                        { { 0 } },
                        this->get_spatial_dim() );
            }
        }

        // ----------------------------------------------------------------------------

        void
        Integration_Mesh_STK::setup_cell_clusters(
                Interpolation_Mesh &aInterpMesh,
                Cell_Cluster_Input *aCellClusterInput )
        {
            // Number of interpolation cells
            moris::uint tNumInterpCells = aInterpMesh.get_num_entities( EntityRank::ELEMENT );

            mCellClusters.resize( tNumInterpCells );

            moris::mtk::Cell_Info_Factory tFactory;

            if ( aCellClusterInput != nullptr )
            {
                // iterate through cells
                for ( moris::uint i = 0; i < tNumInterpCells; i++ )
                {
                    moris_id tCellId =
                            aInterpMesh.get_glb_entity_id_from_entity_loc_index( (moris_index)i, EntityRank::ELEMENT );

                    moris_index tClusterIndex = aCellClusterInput->get_cluster_index( tCellId );

                    if ( tClusterIndex != MORIS_INDEX_MAX )
                    {
                        // mark as nontrivial
                        mCellClusters( i ).mark_as_nontrivial();

                        mCellClusters( i ).set_interpolation_cell( aCellClusterInput->get_interp_cell( tClusterIndex ) );

                        moris::Matrix< IndexMat > const *tPrimaryCellIds = aCellClusterInput->get_primary_cell_ids( tClusterIndex );
                        mCellClusters( i ).add_primary_integration_cell( this->get_cell_pointers_from_ids( *tPrimaryCellIds ) );

                        moris::Matrix< IndexMat > const *tVoidCellIds = aCellClusterInput->get_void_cell_ids( tClusterIndex );
                        mCellClusters( i ).add_void_integration_cell( this->get_cell_pointers_from_ids( *tVoidCellIds ) );

                        moris::Matrix< IndexMat > const *tVertexIds = aCellClusterInput->get_vertex_in_cluster_ids( tClusterIndex );

                        mCellClusters( i ).add_vertex_to_cluster( this->get_vertex_pointers_from_ids( *tVertexIds ) );

                        moris::Matrix< DDRMat > const *tVertexLocalCoords = aCellClusterInput->get_vertex_local_coords_wrt_interpolation_cell( tClusterIndex );
                        mCellClusters( i ).add_vertex_local_coordinates_wrt_interp_cell( *tVertexLocalCoords );
                    }
                    else
                    {
                        // an assumption is made here that if the cell does not have a cluster then it is related 1 to 1
                        // with the interpolation cell with the same id

                        // interpolation cell
                        mtk::Cell const *tInterpCell = &aInterpMesh.get_mtk_cell( (moris_index)i );
                        mCellClusters( i ).set_interpolation_cell( tInterpCell );

                        // integration cell (only primary cells here)
                        moris_index      tIntegCellIndex = this->get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );
                        mtk::Cell const *tPrimaryCell    = &this->get_mtk_cell( tIntegCellIndex );
                        mCellClusters( i ).add_primary_integration_cell( tPrimaryCell );

                        moris::Cell< moris::mtk::Vertex * >       tVertexIds = tPrimaryCell->get_vertex_pointers();
                        moris::Cell< moris::mtk::Vertex const * > tConstVertexPtrs( tVertexIds.size() );
                        for ( moris::uint iV = 0; iV < tVertexIds.size(); iV++ )
                        {
                            tConstVertexPtrs( iV ) = tVertexIds( iV );
                        }

                        // interpolation order of the interpolation cell
                        enum Interpolation_Order tInterpOrder = tInterpCell->get_interpolation_order();

                        // geometry
                        enum Geometry_Type tGeomType = tInterpCell->get_geometry_type();

                        // Cell info
                        moris::mtk::Cell_Info *tCellInfo = tFactory.create_cell_info( tGeomType, tInterpOrder );

                        Matrix< DDRMat > tXi;
                        tCellInfo->get_loc_coords_of_cell( tXi );

                        mCellClusters( i ).add_vertex_to_cluster( tConstVertexPtrs );

                        mCellClusters( i ).add_vertex_local_coordinates_wrt_interp_cell( tXi );

                        delete tCellInfo;
                    }
                }
            }
            else
            {
                for ( moris::uint i = 0; i < tNumInterpCells; i++ )
                {

                    moris_id tCellId = aInterpMesh.get_glb_entity_id_from_entity_loc_index( (moris_index)i, EntityRank::ELEMENT );

                    // interpolation cell
                    mtk::Cell const *tInterpCell = &aInterpMesh.get_mtk_cell( (moris_index)i );
                    mCellClusters( i ).set_interpolation_cell( tInterpCell );

                    // integration cell (only primary cells here)
                    moris_index      tIntegCellIndex = this->get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );
                    mtk::Cell const *tPrimaryCell    = &this->get_mtk_cell( tIntegCellIndex );
                    mCellClusters( i ).add_primary_integration_cell( tPrimaryCell );

                    moris::Cell< moris::mtk::Vertex * >       tVertexIds = tPrimaryCell->get_vertex_pointers();
                    moris::Cell< moris::mtk::Vertex const * > tConstVertexPtrs( tVertexIds.size() );
                    for ( moris::uint iV = 0; iV < tVertexIds.size(); iV++ )
                    {
                        tConstVertexPtrs( iV ) = tVertexIds( iV );
                    }

                    // interpolation order of the interpolation cell
                    enum Interpolation_Order tInterpOrder = tInterpCell->get_interpolation_order();

                    // geometry
                    enum Geometry_Type tGeomType = tInterpCell->get_geometry_type();

                    // Cell info
                    moris::mtk::Cell_Info *tCellInfo = tFactory.create_cell_info( tGeomType, tInterpOrder );

                    Matrix< DDRMat > tXi;
                    tCellInfo->get_loc_coords_of_cell( tXi );

                    mCellClusters( i ).add_vertex_to_cluster( tConstVertexPtrs );

                    mCellClusters( i ).add_vertex_local_coordinates_wrt_interp_cell( tXi );

                    delete tCellInfo;
                }
            }
        }

        void
        Integration_Mesh_STK::setup_blockset_with_cell_clusters()
        {
            // construct integration to cell cluster relationship
            moris::Cell< moris::moris_index > tPrimaryIntegrationCellToClusterIndex(
                    this->get_num_entities( EntityRank::ELEMENT ), MORIS_INDEX_MAX );

            // iterate through clusters
            for ( moris::uint i = 0; i < mCellClusters.size(); i++ )
            {
                Cell_Cluster_STK const &tCellCluster = mCellClusters( i );

                moris::Cell< moris::mtk::Cell const * > const &tPrimaryCells =
                        tCellCluster.get_primary_cells_in_cluster();

                // iterate through primary cells
                for ( moris::uint j = 0; j < tCellCluster.get_num_primary_cells(); j++ )
                {
                    moris::moris_index tCellIndex = tPrimaryCells( j )->get_index();

                    MORIS_ASSERT( tPrimaryIntegrationCellToClusterIndex( tCellIndex ) == MORIS_INDEX_MAX,
                            "Integration cell can only appear as a primary cell in one cell cluster" );

                    tPrimaryIntegrationCellToClusterIndex( tCellIndex ) = (moris_index)i;
                }
            }

            // get all block sets from the mesh
            moris::Cell< std::string > tBlockSetNames = this->get_set_names( EntityRank::ELEMENT );

            mPrimaryBlockSetClusters.resize( tBlockSetNames.size() );
            mPrimaryBlockSetNames = tBlockSetNames;

            moris::Cell< uint > tSetsToRemove;

            for ( moris::uint i = 0; i < tBlockSetNames.size(); i++ )
            {

                moris::Matrix< moris::IndexMat > tCellsInSet =
                        this->get_set_entity_loc_inds( EntityRank::ELEMENT, tBlockSetNames( i ) );

                bool tSetHasCluster = false;

                // add the cells in set's primary cluster to member data (make unique later)
                for ( moris::uint j = 0; j < tCellsInSet.numel(); j++ )
                {
                    moris::moris_index tCellIndex    = tCellsInSet( j );
                    moris::moris_index tClusterIndex = tPrimaryIntegrationCellToClusterIndex( tCellIndex );

                    if ( tClusterIndex != MORIS_INDEX_MAX )
                    {
                        tSetHasCluster = true;
                        mPrimaryBlockSetClusters( i ).push_back( tClusterIndex );
                    }
                }

                // if there are no elements which are a primary cell in one cluster for this block set, then we remove it.
                if ( !tSetHasCluster )
                {
                    tSetsToRemove.push_back( i );
                }

                // remove duplicates
                else
                {
                    moris::unique( mPrimaryBlockSetClusters( i ) );
                }
            }

            // remove block sets which had not primary clusters
            for ( moris::uint i = tSetsToRemove.size(); i > 0; i-- )
            {
                //         mPrimaryBlockSetClusters.erase(tSetsToRemove(i-1));
                //         mPrimaryBlockSetNames.erase(tSetsToRemove(i-1));
            }

            // iterate and create map
            mListofBlocks.resize( mPrimaryBlockSetClusters.size(), nullptr );

            moris::Cell< std::string > tBSNames = this->get_block_set_names();

            for ( moris::uint Ik = 0; Ik < mListofBlocks.size(); Ik++ )
            {
                mListofBlocks( Ik ) = new moris::mtk::Block(
                        tBSNames( Ik ),
                        this->get_cell_clusters_in_set( Ik ),
                        { { 0 } },
                        this->get_spatial_dim() );

                MORIS_ASSERT( mBlockSetLabelToOrd.find( mPrimaryBlockSetNames( Ik ) ) == mBlockSetLabelToOrd.end(),
                        "Duplicate block set in mesh" );

                mBlockSetLabelToOrd[ mPrimaryBlockSetNames( Ik ) ] = Ik;
            }
        }

        void
        Integration_Mesh_STK::setup_side_set_clusters(
                Interpolation_Mesh &aInterpMesh,
                Side_Cluster_Input *aSideClusterInput )
        {

            enum EntityRank tSideSetRank = this->get_facet_rank();

            moris::Cell< std::string > aSideSetNames = this->get_set_names( tSideSetRank );

            mSideSets.resize( aSideSetNames.size() );

            // copy strings labels
            mSideSetLabels.append( aSideSetNames );

            // add to map
            for ( moris::uint i = 0; i < aSideSetNames.size(); i++ )
            {
                MORIS_ASSERT( mSideSideSetLabelToOrd.find( mSideSetLabels( i ) ) == mSideSideSetLabelToOrd.end(),
                        "Duplicate side set label detected." );

                mSideSideSetLabelToOrd[ mSideSetLabels( i ) ] = i;
            }

            // iterate through block sets
            for ( moris::uint i = 0; i < aSideSetNames.size(); i++ )
            {
                // get the cells and side ordinals from the mesh for this side set
                moris::Cell< mtk::Cell const * > tCellsInSet( 0 );
                moris::Matrix< moris::IndexMat > tSideOrdsInSet( 0, 0 );
                this->get_sideset_cells_and_ords( aSideSetNames( i ), tCellsInSet, tSideOrdsInSet );

                moris::mtk::Cell_Info_Factory tFactory;

                // figure out which integration cells are in the side cluster input. these are assumed
                // the only non-trivial ones, all others will be marked as trivial
                if ( aSideClusterInput != nullptr )
                {
                    moris::moris_index tSideClusterOrd = aSideClusterInput->get_side_label_ordinal( aSideSetNames( i ) );

                    if ( tSideClusterOrd == MORIS_INDEX_MAX )
                    {
                        // loop over cells in the side set and make sure they have all been included
                        for ( moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++ )
                        {
                            // integration cell id
                            moris_id tCellId = tCellsInSet( iIGCell )->get_id();

                            // interpolation cell index
                            moris_index tCellIndex = aInterpMesh.get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );

                            // get the mtk cell
                            moris::mtk::Cell *tInterpCell = &aInterpMesh.get_mtk_cell( tCellIndex );

                            // interpolation order of the interpolation cell
                            enum Interpolation_Order tInterpOrder = tInterpCell->get_interpolation_order();

                            // geometry
                            enum Geometry_Type tGeomType = tInterpCell->get_geometry_type();

                            // Cell info
                            moris::mtk::Cell_Info *tCellInfo = tFactory.create_cell_info( tGeomType, tInterpOrder );

                            // get local coordinates on the side
                            Matrix< DDRMat > tXi;
                            tCellInfo->get_loc_coord_on_side_ordinal( tSideOrdsInSet( iIGCell ), tXi );

                            mSideSets( i ).push_back( Side_Cluster_STK( true,
                                    tInterpCell,
                                    { tCellsInSet( iIGCell ) },
                                    { { tSideOrdsInSet( iIGCell ) } },
                                    tCellsInSet( iIGCell )->get_vertices_on_side_ordinal( tSideOrdsInSet( iIGCell ) ),
                                    tXi ) );
                            delete tCellInfo;
                        }
                    }
                    else
                    {
                        // access side set cluster data
                        Side_Set_Cluster_Data const &tSideSetClusterData = aSideClusterInput->get_cluster_data( tSideClusterOrd );

                        // mark integration cells which are not trivial
                        std::unordered_map< moris_index, bool > tIntegrationCellsInSideSet;

                        // collect integration cells in clusters for this side set
                        uint tNumClusters = tSideSetClusterData.get_num_cell_clusters();

                        for ( moris::uint iC = 0; iC < tNumClusters; iC++ )
                        {
                            // get data from the side set cluster data

                            // cell ids and side ords
                            moris::Matrix< moris::IdMat > const    *tCellIdsAndOrds = tSideSetClusterData.get_integration_cell_ids_and_side_ords( iC );
                            moris::Cell< moris::mtk::Cell const * > tCellPointers   = this->get_cell_pointers_from_ids( tCellIdsAndOrds->get_column( 0 ) );
                            moris::Matrix< moris::IdMat >           tSideOrds       = tCellIdsAndOrds->get_column( 1 );

                            // vertex pointers
                            moris::Matrix< moris::IdMat > const      *tVertexInCluster = tSideSetClusterData.get_vertex_in_cluster_ids( iC );
                            moris::Cell< moris::mtk::Vertex const * > tVertices        = this->get_vertex_pointers_from_ids( *tVertexInCluster );

                            // vertices in cluster
                            // vertex parametric coordinate relative to
                            mSideSets( i ).push_back( Side_Cluster_STK(
                                    false,
                                    tSideSetClusterData.get_interp_cell( iC ),
                                    tCellPointers,
                                    tSideOrds,
                                    tVertices,
                                    *tSideSetClusterData.get_vertex_local_coords_wrt_interpolation_cell( iC ) ) );

                            // mark all integration cells in this cluster
                            for ( moris::uint iIGCell = 0; iIGCell < tCellIdsAndOrds->n_rows(); iIGCell++ )
                            {
                                tIntegrationCellsInSideSet[ ( *tCellIdsAndOrds )( iIGCell, 0 ) ] = true;
                            }
                        }

                        // loop over cells in the side set and make sure they have all been included
                        //            for(moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++)
                        //            {
                        //                moris_id tCellId = tCellsInSet(iIGCell)->get_id();
                        //
                        //                if(tIntegrationCellsInSideSet.find(tCellId) == tIntegrationCellsInSideSet.end())
                        //                {
                        //                    // interpolation cell index
                        //                    moris_index tCellIndex = aInterpMesh.get_loc_entity_ind_from_entity_glb_id(tCellId,EntityRank::ELEMENT);
                        //
                        //                    // construct a trivial side cluster
                        //                    moris::mtk::Cell* tInterpCell = &aInterpMesh.get_mtk_cell(tCellIndex);
                        //                    mSideSets(i).push_back(Side_Cluster_STK(tInterpCell,tCellsInSet(iIGCell), tCellsInSet(iIGCell)->get_vertices_on_side_ordinal(tSideOrdsInSet(iIGCell)), tSideOrdsInSet(iIGCell)));
                        //
                        //                }
                        //            }
                    }
                }
                // all trivial case
                else
                {
                    // loop over cells in the side set and make sure they have all been included
                    for ( moris::uint iIGCell = 0; iIGCell < tCellsInSet.size(); iIGCell++ )
                    {
                        moris_id tCellId = tCellsInSet( iIGCell )->get_id();

                        // interpolation cell index
                        moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id( tCellId, EntityRank::ELEMENT );

                        // construct a trivial side cluster
                        moris::mtk::Cell_STK *tInterpCell = &mSTKMeshData->mMtkCells( tCellIndex );

                        // get the cell info of the interp cell
                        mtk::Cell_Info const *tCellInfo = tInterpCell->get_cell_info();

                        // get local coordinates on the side
                        Matrix< DDRMat > tXi;
                        tCellInfo->get_loc_coord_on_side_ordinal( tSideOrdsInSet( iIGCell ), tXi );

                        mSideSets( i ).push_back(
                                Side_Cluster_STK(
                                        true,
                                        tInterpCell,
                                        { tCellsInSet( iIGCell ) },
                                        { { tSideOrdsInSet( iIGCell ) } },
                                        tCellsInSet( iIGCell )->get_vertices_on_side_ordinal( tSideOrdsInSet( iIGCell ) ),
                                        tXi ) );
                    }
                }
            }

            mListofSideSets.resize( mSideSets.size(), nullptr );

            for ( moris::uint Ik = 0; Ik < mListofSideSets.size(); Ik++ )
            {
                mListofSideSets( Ik ) = new moris::mtk::Side_Set(
                        aSideSetNames( Ik ),
                        this->get_side_set_cluster( Ik ),
                        { { 0 } },
                        this->get_spatial_dim() );
            }
        }
        // ----------------------------------------------------------------------------

        void
        Integration_Mesh_STK::setup_double_side_set_clusters(
                Interpolation_Mesh        &aInterpMesh,
                Double_Side_Cluster_Input *aDoubleSideClusterInput )
        {
            moris::Cell< std::string > const &tDoubleSideSetLabels = aDoubleSideClusterInput->get_double_side_set_labels();

            // copy strings labels
            mDoubleSideSetLabels.append( tDoubleSideSetLabels );

            // copy the map
            mDoubleSideSetLabelToOrd = aDoubleSideClusterInput->mSideSetLabelToOrd;

            // resize member data
            mDoubleSideSets.resize( tDoubleSideSetLabels.size() );

            // iterate through double side sets
            for ( moris::uint i = 0; i < tDoubleSideSetLabels.size(); i++ )
            {
                // access left and right side clusters
                Side_Set_Cluster_Data const &tLeftClusterData  = aDoubleSideClusterInput->mLeftSideClusters( i );
                Side_Set_Cluster_Data const &tRightClusterData = aDoubleSideClusterInput->mRightSideClusters( i );

                MORIS_ASSERT( tLeftClusterData.get_num_cell_clusters() == tRightClusterData.get_num_cell_clusters(),
                        "Mismatch between left sides and right sides when constructing double sided cluster" );

                // number of side clusters in this side set
                moris_index tNumClusters = tLeftClusterData.get_num_cell_clusters();

                for ( moris::uint iC = 0; iC < (uint)tNumClusters; iC++ )
                {
                    // construct a single side cluster for the left

                    bool        tLeftTrivial = tLeftClusterData.is_trivial( iC );
                    moris_index tLeftIndex   = mDoubleSideSetSideClusters.size();

                    // if trivial
                    if ( tLeftTrivial )
                    {
                        // Get interpolation cell of the left
                        moris::mtk::Cell_STK const *tInterpCell =
                                reinterpret_cast< moris::mtk::Cell_STK const * >( tLeftClusterData.get_interp_cell( iC ) );

                        // integration cell and side ordinals
                        moris::Matrix< IndexMat > const *tIntegCellId =
                                tLeftClusterData.get_integration_cell_ids_and_side_ords( iC );

                        MORIS_ASSERT( tIntegCellId->numel() == 2,
                                "more than one integration cell in interpolation cluster" );

                        moris_index tCellIndex =
                                this->get_loc_entity_ind_from_entity_glb_id( ( *tIntegCellId )( 0 ), EntityRank::ELEMENT );

                        moris::mtk::Cell const *tIntegCell = &this->get_mtk_cell( tCellIndex );
                        moris_index             tSideOrd   = ( *tIntegCellId )( 1 );

                        moris::Cell< moris::mtk::Vertex const * > tVerticesOnSide = tIntegCell->get_vertices_on_side_ordinal( tSideOrd );

                        // get the cell info of the interp cell
                        mtk::Cell_Info const *tCellInfo = tInterpCell->get_cell_info();

                        // get local coordinates on the side
                        Matrix< DDRMat > tXi;
                        tCellInfo->get_loc_coord_on_side_ordinal( tSideOrd, tXi );

                        mDoubleSideSetSideClusters.push_back(
                                Side_Cluster_STK(
                                        true,
                                        tInterpCell,
                                        { tIntegCell },
                                        { { tSideOrd } },
                                        tVerticesOnSide,
                                        tXi ) );
                    }
                    // if not trivial
                    else
                    {
                        MORIS_ERROR( 0, "Non-trivial not implemented" );
                    }

                    bool        tRightTrivial = tRightClusterData.is_trivial( iC );
                    moris_index tRightIndex   = mDoubleSideSetSideClusters.size();
                    // if trivial
                    if ( tRightTrivial )
                    {
                        // Get interpolation cell of the left
                        moris::mtk::Cell_STK const *tInterpCell =
                                reinterpret_cast< moris::mtk::Cell_STK const * >( tRightClusterData.get_interp_cell( iC ) );

                        // integration cell and side ordinals
                        moris::Matrix< IndexMat > const *tIntegCellId = tRightClusterData.get_integration_cell_ids_and_side_ords( iC );

                        MORIS_ASSERT( tIntegCellId->numel() == 2, "more than one integration cell in interpolation cluster" );

                        // integration cell
                        moris_index tCellIndex =
                                this->get_loc_entity_ind_from_entity_glb_id( ( *tIntegCellId )( 0 ), EntityRank::ELEMENT );
                        moris_index             tSideOrd   = ( *tIntegCellId )( 1 );
                        moris::mtk::Cell const *tIntegCell = &this->get_mtk_cell( tCellIndex );

                        moris::Cell< moris::mtk::Vertex const * > tVerticesOnSide = tIntegCell->get_vertices_on_side_ordinal( tSideOrd );

                        // get the cell info of the interp cell
                        mtk::Cell_Info const *tCellInfo = tInterpCell->get_cell_info();

                        // get local coordinates on the side
                        Matrix< DDRMat > tXi;
                        tCellInfo->get_loc_coord_on_side_ordinal( tSideOrd, tXi );

                        mDoubleSideSetSideClusters.push_back(
                                Side_Cluster_STK(
                                        true,
                                        tInterpCell,
                                        { tIntegCell },
                                        { { tSideOrd } },
                                        tVerticesOnSide,
                                        tXi ) );
                    }
                    // if not trivial
                    else
                    {
                        MORIS_ERROR( 0, "Non-trivial not implemented" );
                    }

                    // construct left to right vertex pairing
                    moris::Matrix< moris::IdMat >            *tVertexPairing = aDoubleSideClusterInput->mVertexPairing( i )( iC );
                    moris::uint                               tNumVertPairs  = tVertexPairing->n_rows();
                    moris::Cell< moris::mtk::Vertex const * > tVertexLeftToRightPair( tNumVertPairs );

                    // left vertex pointers
                    moris::Cell< moris::mtk::Vertex const * > tLeftVertPtrs =
                            this->get_vertex_pointers_from_ids( tVertexPairing->get_column( 0 ) );
                    moris::Cell< moris::mtk::Vertex const * > tRightVertPtrs =
                            this->get_vertex_pointers_from_ids( tVertexPairing->get_column( 1 ) );

                    for ( moris::uint iV = 0; iV < tNumVertPairs; iV++ )
                    {
                        // get left vertex index
                        moris_index tLeftClustIndex =
                                mDoubleSideSetSideClusters( tLeftIndex ).get_vertex_cluster_index( tLeftVertPtrs( iV ) );

                        tVertexLeftToRightPair( tLeftClustIndex ) = tRightVertPtrs( iV );
                    }

                    // construct the double side cluster
                    mDoubleSideSets( i ).push_back(
                            new Double_Side_Cluster(
                                    &mDoubleSideSetSideClusters( tLeftIndex ),
                                    &mDoubleSideSetSideClusters( tRightIndex ),
                                    tVertexLeftToRightPair ) );
                }
            }
            mListofDoubleSideSets.resize( mDoubleSideSets.size(), nullptr );

            for ( moris::uint Ik = 0; Ik < mListofDoubleSideSets.size(); Ik++ )
            {
                mListofDoubleSideSets( Ik ) = new moris::mtk::Double_Side_Set(
                        mDoubleSideSetLabels( Ik ),
                        this->get_double_side_set_cluster( Ik ),
                        { { 0 } },
                        this->get_spatial_dim() );
            }
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::mtk::Cell const * >
        Integration_Mesh_STK::get_cell_pointers_from_ids( moris::Matrix< moris::IdMat > const &aCellIds ) const
        {

            moris::Cell< moris::mtk::Cell const * > tCellPtrs( aCellIds.numel() );

            for ( moris::uint i = 0; i < aCellIds.numel(); i++ )
            {
                moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id( aCellIds( i ), EntityRank::ELEMENT );
                tCellPtrs( i )         = &this->get_mtk_cell( tCellIndex );
            }

            return tCellPtrs;
        }

        // ----------------------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex const * >
        Integration_Mesh_STK::get_vertex_pointers_from_ids( moris::Matrix< moris::IdMat > const &aVertexIds ) const
        {
            moris::Cell< moris::mtk::Vertex const * > tVertexPtrs( aVertexIds.numel() );

            for ( moris::uint i = 0; i < aVertexIds.numel(); i++ )
            {
                moris_index tCellIndex = this->get_loc_entity_ind_from_entity_glb_id( aVertexIds( i ), EntityRank::NODE );
                tVertexPtrs( i )       = &this->get_mtk_vertex( tCellIndex );
            }

            return tVertexPtrs;
        }    // ----------------------------------------------------------------------------

    }    // namespace mtk
}    // namespace moris
