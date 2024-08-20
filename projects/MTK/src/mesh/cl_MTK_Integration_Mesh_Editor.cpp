/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Mesh_Editor.cpp
 *
 */

#include "cl_MTK_Integration_Mesh_Editor.hpp"
#include "cl_Tracer.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Mesh_DataBase_IP.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Side_Cluster_DataBase.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Set_Communicator.hpp"
#include "cl_Tracer.hpp"

namespace moris::mtk
{
    // ----------------------------------------------------------------------------

    Integration_Mesh_Editor::Integration_Mesh_Editor(
            moris::mtk::Integration_Mesh*               aMTKMesh,
            moris::mtk::Interpolation_Mesh_DataBase_IP* aIPMeshDataBase,
            bool                                        aCheckMesh )
            : mInputMesh( aMTKMesh )
            , mIPMeshDataBase( aIPMeshDataBase )
            , mCheckMesh( aCheckMesh )
    {
        mIGMeshInfo              = new Integration_Mesh_Info;
        mIGMeshInfo->mVertices   = mInputMesh->get_all_vertices();
        mIGMeshInfo->mSpatialDim = mInputMesh->get_spatial_dim();
    }

    // ----------------------------------------------------------------------------

    Integration_Mesh_Editor::~Integration_Mesh_Editor()
    {
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_vertex_data()
    {
        this->initialize_vertex_data();

        for ( uint iVertex = 0; iVertex < mIGMeshInfo->mVertices.size(); iVertex++ )
        {
            // Assign the coordinates of the lagrange IP mesh
            Matrix< DDRMat > tCoords = mIGMeshInfo->mVertices( iVertex )->get_coords();

            // write the matrix data in the memory slot
            std::copy( tCoords.data(), tCoords.data() + mIGMeshInfo->mSpatialDim, mOutputMesh->mVertexCoordinates.colptr( iVertex ) );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::initialize_vertex_data()
    {
        mOutputMesh->mVertexCoordinates.set_size( mIGMeshInfo->mSpatialDim, mIGMeshInfo->mVertices.size() );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_cell_data()
    {
        // initialize the size of the connectivity info
        this->initialize_cell_data();

        // loop over to populate the index of the vertices in the cell to vertex connectivity
        for ( size_t iCell = 0; iCell < mInputMesh->get_num_elems(); iCell++ )
        {
            // get the mtk cell from the mesh
            mtk::Cell const & tCell = mInputMesh->get_mtk_cell( iCell );

            // get the vertex pointers to extract the
            Vector< Vertex* > tVertexPointers = tCell.get_vertex_pointers();

            // use the transform to transfer vertex to index and write them in the correct spot of the connectivity info
            std::transform( tVertexPointers.begin(),
                    tVertexPointers.end(),
                    mIGMeshInfo->mCellToVertexIndicies.begin() + mOutputMesh->mCellToVertexOffSet( iCell ),
                    []( Vertex*& aVertex ) -> moris_index { return aVertex->get_index(); } );
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::initialize_cell_data()
    {
        // get number of cells
        uint tNumCells = mInputMesh->get_num_elems();

        // determine size of the offset value
        mOutputMesh->mCellToVertexOffSet.resize( tNumCells + 1 );

        // the first offset is always zero
        size_t tOffSet                        = 0;
        mOutputMesh->mCellToVertexOffSet( 0 ) = tOffSet;

        // loop over the cells to fill out the offset
        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            tOffSet += mInputMesh->get_mtk_cell( iCell ).get_number_of_vertices();
            mOutputMesh->mCellToVertexOffSet( iCell + 1 ) = tOffSet;
        }

        // reserve enough space ( resize has to be used because we will be using std::transform on this)
        mIGMeshInfo->mCellToVertexIndicies.resize( tOffSet );
    }

    // ----------------------------------------------------------------------------
    void
    Integration_Mesh_Editor::initialize_cell_cluster_data()
    {
        //  get number of clusters,  there is a 1 to 1 relation between IP cells and clusters
        uint tNumCellClusters = mIPMeshDataBase->get_num_elems();

        // determine the size of the connectivity
        mOutputMesh->mCellClusterToPrimaryIGCellOffSet.resize( tNumCellClusters + 1 );
        mOutputMesh->mCellClusterToVoidIGCellOffset.resize( tNumCellClusters + 1 );
        mOutputMesh->mCellClusterToVertexOffset.resize( tNumCellClusters + 1 );

        // declare and initialize the offset values
        size_t tPrimaryCellOffset     = 0;
        size_t tVoidCellOffset        = 0;
        size_t tVertexOffset          = 0;
        size_t tLocalCoordOffset      = 0;
        size_t tNumNonTrivialClusters = 0;
        // size_t tNonTrivialVertexOffSet = 0;

        // first offset value if always zero
        mOutputMesh->mCellClusterToPrimaryIGCellOffSet( 0 ) = 0;
        mOutputMesh->mCellClusterToVoidIGCellOffset( 0 )    = 0;
        mOutputMesh->mCellClusterToVertexOffset( 0 )        = 0;

        // loop over the clusters to get the number of the connected entities
        for ( moris_index iCluster = 0; iCluster < (moris_index)tNumCellClusters; iCluster++ )
        {
            Cell_Cluster const & tCluster = mInputMesh->get_cell_cluster( iCluster );
            // increment the entties
            tPrimaryCellOffset += tCluster.get_num_primary_cells();
            tVoidCellOffset += tCluster.get_num_void_cells();
            tVertexOffset += tCluster.get_num_vertices_in_cluster();

            // this is to determine if the cluster is trivial not to store the local coordinates as they can be obtained from the cell info
            bool tIsTrivial = tCluster.is_trivial();

            // assign the offset coefficients
            mOutputMesh->mCellClusterToPrimaryIGCellOffSet( iCluster + 1 ) = tPrimaryCellOffset;
            mOutputMesh->mCellClusterToVoidIGCellOffset( iCluster + 1 )    = tVoidCellOffset;
            mOutputMesh->mCellClusterToVertexOffset( iCluster + 1 )        = tVertexOffset;

            // if it is trivial then add so we know the size
            if ( !tIsTrivial )
            {

                // MORIS_ASSERT( tCluster.get_num_vertices_in_cluster() == tCluster.get_vertices_local_coordinates_wrt_interp_cell().n_rows(), "DimMistach" );

                tLocalCoordOffset += mInputMesh->get_cell_cluster( iCluster ).get_num_vertices_in_cluster();

                tNumNonTrivialClusters++;
            }
        }

        // Resize the connectivity info for the cluster based on the right size
        // need to call resize since we will be using transform on these cells
        mIGMeshInfo->mCellClusterToPrimaryIGCellIndices.resize( tPrimaryCellOffset );
        mIGMeshInfo->mCellClusterToVoidIGCellIndices.resize( tVoidCellOffset );
        mIGMeshInfo->mCellClusterToVertexIndices.resize( tVertexOffset );

        // set the size of the matrix pointer
        mOutputMesh->mCellClusterVertexCoords = new moris::Matrix< moris::DDRMat >( tLocalCoordOffset, mIGMeshInfo->mSpatialDim );

        // reserve enough space on the map
        mOutputMesh->mCellClusterIndexToRowNumber.reserve( tNumNonTrivialClusters );
    }

    // ----------------------------------------------------------------------------
    void
    Integration_Mesh_Editor::generate_cell_cluster_data()
    {
        // initialize the correct sizes
        this->initialize_cell_cluster_data();

        //  get number of clusters,  there is a 1 to 1 relation between IP cells and clusters
        uint tNumCellClusters = mIPMeshDataBase->get_num_elems();

        // This is to keep track of the incremt of the coordinate matrix size
        // we make the important assumption that clusters are ordered consecutively
        size_t tLocalCoordsOffset = 0;

        // loop over the cluster to get the index of the cells and append them to the connectivity
        for ( moris_index iCluster = 0; iCluster < (moris_index)tNumCellClusters; iCluster++ )
        {
            // get the cell cluster
            Cell_Cluster const & tCluster = mInputMesh->get_cell_cluster( iCluster );

            // get the primary , void and vertex pointer of the cells to extact their index
            Vector< moris::mtk::Cell const * > const & tPrimaryCells = tCluster.get_primary_cells_in_cluster();
            Vector< moris::mtk::Cell const * > const & tVoidCells    = tCluster.get_void_cells_in_cluster();
            Vector< moris::mtk::Vertex const * >       tVertex       = tCluster.get_vertices_in_cluster();

            // fill out the index of the connectivity
            std::transform(
                    tPrimaryCells.cbegin(),
                    tPrimaryCells.cend(),
                    mIGMeshInfo->mCellClusterToPrimaryIGCellIndices.begin() + mOutputMesh->mCellClusterToPrimaryIGCellOffSet( iCluster ),
                    []( moris::mtk::Cell const * aCell ) {
                        return aCell->get_index();
                    } );

            std::transform(
                    tVoidCells.cbegin(),
                    tVoidCells.cend(),
                    mIGMeshInfo->mCellClusterToVoidIGCellIndices.begin() + mOutputMesh->mCellClusterToVoidIGCellOffset( iCluster ),
                    []( moris::mtk::Cell const * aCell ) {
                        return aCell->get_index();
                    } );

            std::transform(
                    tVertex.begin(),
                    tVertex.end(),
                    mIGMeshInfo->mCellClusterToVertexIndices.begin() + mOutputMesh->mCellClusterToVertexOffset( iCluster ),
                    []( moris::mtk::Vertex const * aVertex ) {
                        return aVertex->get_index();
                    } );

            // if the cluster is non-trivial (i.e. cut base IP cell)
            if ( !tCluster.is_trivial() )
            {
                // if non trivial store the coordinates
                moris::Matrix< moris::DDRMat > tVertexCoords = mInputMesh->get_cell_cluster( iCluster ).get_vertices_local_coordinates_wrt_interp_cell();

                // transfer the coordinates to the right block in the matrix
                mOutputMesh->mCellClusterVertexCoords->operator()( { tLocalCoordsOffset, tLocalCoordsOffset + tVertexCoords.n_rows() - 1 }, { 0, mIGMeshInfo->mSpatialDim - 1 } ) = tVertexCoords.matrix_data();

                // add data to the map
                mOutputMesh->mCellClusterIndexToRowNumber[ iCluster ] = (moris_index)tLocalCoordsOffset;

                // increase the count of non-trivial ones
                tLocalCoordsOffset += tVertexCoords.n_rows();
            }
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::initialize_side_cluster_data()
    {
        // get the side sets from th mesh
        Vector< moris::mtk::Side_Set* > const & tSideSets = mInputMesh->get_side_sets();

        // initialize side set info offSet
        mIGMeshInfo->mSideSetToSideClusterOffset.resize( tSideSets.size() + 1 );

        // initialize number of side clusters
        mOutputMesh->mNumSideClusters = 0;

        // initialize the set counter
        size_t iCounterSet = 0;

        mIGMeshInfo->mSideSetToSideClusterOffset( iCounterSet ) = mOutputMesh->mNumSideClusters;

        // loop over the side sets to sum up the individual clusters
        for ( const auto& iSet : tSideSets.data() )
        {
            mOutputMesh->mNumSideClusters += iSet->get_num_clusters_on_set();

            mIGMeshInfo->mSideSetToSideClusterOffset( iCounterSet + 1 ) = mOutputMesh->mNumSideClusters;

            iCounterSet++;
        }

        // Based on the number of side cluster initialize the offset cells
        mOutputMesh->mSideClusterToPrimaryIGCellOffset.resize( mOutputMesh->mNumSideClusters + 1 );
        mOutputMesh->mSideClusterToVoidIGCellOffset.resize( mOutputMesh->mNumSideClusters + 1 );
        mOutputMesh->mSideClusterToVertexOffSet.resize( mOutputMesh->mNumSideClusters + 1 );

        // initialize the first entry
        size_t tPrimaryCellOffset = 0;
        size_t tVoidCellOffset    = 0;
        size_t tVertexOffSet      = 0;

        // set the first entry (always zero)
        mOutputMesh->mSideClusterToPrimaryIGCellOffset( 0 ) = tPrimaryCellOffset;
        mOutputMesh->mSideClusterToVertexOffSet( 0 )        = tVertexOffSet;

        uint iCounter = 0;

        // loop over sets
        for ( auto const & iSet : tSideSets.data() )
        {
            Vector< Cluster const * > tSideClusters = iSet->get_clusters_on_set();

            // loop over clusters to populate them with the index
            for ( auto iCluster : tSideClusters )
            {
                tPrimaryCellOffset += iCluster->get_num_primary_cells();
                // tVoidCellOffset += iCluster->get_num_void_cells();
                tVertexOffSet += iCluster->get_num_vertices_in_cluster();

                mOutputMesh->mSideClusterToPrimaryIGCellOffset( iCounter + 1 ) = tPrimaryCellOffset;
                mOutputMesh->mSideClusterToVoidIGCellOffset( iCounter + 1 )    = tVoidCellOffset;
                mOutputMesh->mSideClusterToVertexOffSet( iCounter + 1 )        = tVertexOffSet;

                iCounter++;
            }
        }

        // use resize instead of reserve because we will use transfrom function
        mIGMeshInfo->mSideClusterToPrimaryIGCellIndices.resize( tPrimaryCellOffset );
        mIGMeshInfo->mSideClusterToVoidIGCellIndices.resize( tVoidCellOffset );
        mIGMeshInfo->mSideClusterToVertexIndices.resize( tVertexOffSet );

        // we will use .insert on this so reserve is enough
        mOutputMesh->mSideClusterToPrimaryIGCellSideOrd.reserve( tPrimaryCellOffset );
        mOutputMesh->mSideClusterToIPCell.reserve( mOutputMesh->mNumSideClusters );

        // reserve enough space for the map
        mIGMeshInfo->mPreviousSideClusterToNewSideCluster.reserve( mOutputMesh->mNumSideClusters );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_side_cluster_data()
    {
        this->initialize_side_cluster_data();

        // get the sets
        Vector< moris::mtk::Side_Set* > const & tSideSets = mInputMesh->get_side_sets();

        // initialize the cluster counter
        uint iCounter = 0;

        for ( auto&& iSet : tSideSets.data() )
        {
            Vector< Cluster const * > tSideClusters = iSet->get_clusters_on_set();

            // loop over the clusters
            for ( const auto& iCluster : tSideClusters )
            {
                // add the cluster address to the map
                mIGMeshInfo->mPreviousSideClusterToNewSideCluster[ iCluster ] = iCounter;

                // get the member data
                moris::Matrix< moris::IndexMat >           tSideOrdinal  = iCluster->get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCells = iCluster->get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertex       = iCluster->get_vertices_in_cluster();
                uint                                       tIPCellIndex  = iCluster->get_interpolation_cell_index();

                // convert them to indices and append them to the big cell of data
                mOutputMesh->mSideClusterToPrimaryIGCellSideOrd.insert( mOutputMesh->mSideClusterToPrimaryIGCellOffset( iCounter ), tSideOrdinal.begin(), tSideOrdinal.end() );

                std::transform(
                        tPrimaryCells.cbegin(),
                        tPrimaryCells.cend(),
                        mIGMeshInfo->mSideClusterToPrimaryIGCellIndices.begin() + mOutputMesh->mSideClusterToPrimaryIGCellOffset( iCounter ),
                        []( moris::mtk::Cell const * aCell ) {
                            return aCell->get_index();
                        } );

                std::transform(
                        tVertex.begin(),
                        tVertex.end(),
                        mIGMeshInfo->mSideClusterToVertexIndices.begin() + mOutputMesh->mSideClusterToVertexOffSet( iCounter ),
                        []( moris::mtk::Vertex const * aVertex ) {
                            return aVertex->get_index();
                        } );

                mOutputMesh->mSideClusterToIPCell.push_back( tIPCellIndex );

                iCounter++;
            }
        }
    }

    // ---------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_double_sided_cluster_data()
    {
        this->initialize_double_sided_cluster_data();

        // get the old data of the side sets
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // counter for the double sided cluster
        size_t iCounter = 0;

        // loop over the double sided set
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // if ghost, skip
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // get the clusters
            Vector< Cluster const * > tClusters = iSet->get_clusters_on_set();

            // loop over the clusters
            for ( const auto& iCluster : tClusters )
            {
                // leader side cluster and follower side cluster
                moris::mtk::Cluster const & tFollowerSideCluster = iCluster->get_follower_side_cluster();
                moris::mtk::Cluster const & tLeaderSideCluster   = iCluster->get_leader_side_cluster();

                // get the index of side clusters based on the address of the old clusters
                moris_index tLeaderSideIndex   = mIGMeshInfo->mPreviousSideClusterToNewSideCluster[ &tLeaderSideCluster ];
                moris_index tFollowerSideIndex = mIGMeshInfo->mPreviousSideClusterToNewSideCluster[ &tFollowerSideCluster ];

                // make a pair showing index of leader and follower side for the double sided cluster
                mIGMeshInfo->mDoubleSidedClusterToNewSideClusterIndex( iCounter ) = std::make_pair( tLeaderSideIndex, tFollowerSideIndex );

                // get the leader vertex pairs
                Vector< moris::mtk::Vertex const * > const & tVertexPairs = iCluster->get_leader_vertex_pairs();

                // transfer of vertices to the index
                std::transform( tVertexPairs.cbegin(),
                        tVertexPairs.cend(),
                        mIGMeshInfo->mDoubleSidedClusterToPairedVerticesIndex.begin() + mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( iCounter ),
                        []( Vertex const * aVertex ) -> moris_index { return aVertex->get_index(); } );

                // increment the count
                iCounter++;
            }
        }
    }

    //----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::initialize_double_sided_cluster_data()
    {
        // first determine total number of double sided clusters
        size_t tNumDoubleSidedCluster = 0;

        // get the double side set
        Vector< moris::mtk::Double_Side_Set* > const & tSideSets = mInputMesh->get_double_side_sets();

        // set counter
        uint iCounterSet = 0;

        // initialize the offset for double sideset to double side cluster
        mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset.resize( tSideSets.size() + 1 );

        mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( 0 ) = iCounterSet;

        for ( const auto& iSet : tSideSets.data() )
        {
            // if it is a ghost skip
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // update number of double sided cluster
            tNumDoubleSidedCluster += iSet->get_num_clusters_on_set();

            mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( iCounterSet + 1 ) = tNumDoubleSidedCluster;

            iCounterSet++;
        }

        // Propely size the connectivity
        mIGMeshInfo->mDoubleSidedClusterToVertexOffSet.resize( tNumDoubleSidedCluster + 1 );
        mIGMeshInfo->mDoubleSidedClusterToNewSideClusterIndex.resize( tNumDoubleSidedCluster );

        // loop over to assign the connectivity values
        size_t tDoubleSidedClusterToVertexPairOffset = 0;

        // counter for the dbl sided cluster
        size_t iCounter = 0;

        // dbl side cluster to vertex offset
        mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( 0 ) = 0;

        // loop over the sets
        for ( const auto& iSet : tSideSets.data() )
        {
            // if ghost, skip it
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // clusters on the set
            Vector< Cluster const * > tClusters = iSet->get_clusters_on_set();

            // loop over the clusters
            for ( mtk::Cluster const * iCluster : tClusters )
            {
                // update number of vertices
                tDoubleSidedClusterToVertexPairOffset += iCluster->get_leader_vertex_pairs().size();

                // populate the offset
                mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( iCounter + 1 ) = tDoubleSidedClusterToVertexPairOffset;

                iCounter++;
            }
        }

        // reserve enough space for indices
        mIGMeshInfo->mDoubleSidedClusterToPairedVerticesIndex.resize( tDoubleSidedClusterToVertexPairOffset );
    }

    //----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::initialize_ghost_data()
    {
        // get double sided sets
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // number of ghost dbl sided sets
        size_t tNumberOfGhostSideSets = 0;

        // number of clusters
        size_t tNumGhostDblSideClusters = 0;

        // loop over double side sets
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                tNumGhostDblSideClusters += iSet->get_num_clusters_on_set();
                tNumberOfGhostSideSets++;
            }
        };

        mIGMeshInfo->mGhostToVertexOffset.resize( tNumGhostDblSideClusters + 1 );

        mIGMeshInfo->mGhostToVertexOffset( 0 ) = 0;

        mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset.resize( tNumberOfGhostSideSets + 1 );

        uint iCounterSet = 0;
        uint iCounter    = 0;

        // number of double sided clusters for connectivity of sideset to side clusters
        uint tNumDoubleSidedCluster = 0;
        uint tNumGhostVertices      = 0;

        // coordinates initialize
        size_t tLocalCoordsOffSet     = 0;
        size_t tNumNonTrivialClusters = 0;

        mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( 0 ) = iCounterSet;

        // loop over double sided sets
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // make if it is a ghost
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                // clusters on the set
                Vector< Cluster const * > tClusters = iSet->get_clusters_on_set();

                // update double sided clusters number
                tNumDoubleSidedCluster += tClusters.size();

                // offset of the set to the cluster
                mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( iCounterSet + 1 ) = tNumDoubleSidedCluster;

                // loop over the clusters
                for ( mtk::Cluster const * iCluster : tClusters )
                {
                    // number of leader vertex pair size
                    tNumGhostVertices += iCluster->get_leader_vertex_pairs().size();

                    // populate the offset of the ghost to leader vertex pair
                    mIGMeshInfo->mGhostToVertexOffset( iCounter + 1 ) = tNumGhostVertices;

                    // increment the double side cluster count
                    iCounter++;

                    moris::mtk::Cluster const & tLeaderSideCluster   = iCluster->get_leader_side_cluster();
                    moris::mtk::Cluster const & tFollowerSideCluster = iCluster->get_follower_side_cluster();

                    // this is to identify the coordinates
                    if ( !tLeaderSideCluster.is_trivial() )
                    {
                        tLocalCoordsOffSet += tLeaderSideCluster.get_num_vertices_in_cluster();
                        tNumNonTrivialClusters++;
                    }

                    if ( !tFollowerSideCluster.is_trivial() )
                    {
                        tLocalCoordsOffSet += tFollowerSideCluster.get_num_vertices_in_cluster();
                        tNumNonTrivialClusters++;
                    }
                }

                // increment the set counter
                iCounterSet++;
            }
        }

        // set the size of the coordinate matrix
        mOutputMesh->mSecondaryClusterVertexCoords = new moris::Matrix< moris::DDRMat >( tLocalCoordsOffSet, mIGMeshInfo->mSpatialDim );
        mOutputMesh->mSecondaryClusterIndexToRowNumber.reserve( tNumNonTrivialClusters );

        // reserve the ip cell
        mIGMeshInfo->mGhostLeaderToIPCellIndex.resize( tNumDoubleSidedCluster );
        mIGMeshInfo->mGhostFollowerToIPCellIndex.resize( tNumDoubleSidedCluster );

        // reserve the ig cell
        mIGMeshInfo->mGhostLeaderToIGCellIndex.resize( tNumDoubleSidedCluster );
        mIGMeshInfo->mGhostFollowerToIGCellIndex.resize( tNumDoubleSidedCluster );

        // reserve space for side ordinals
        mOutputMesh->mGhostLeaderFollowerOrd.reserve( 2 * tNumDoubleSidedCluster );

        // reserve space for trivial
        mOutputMesh->mGhostLeaderFollowerIsTrivial.reserve( 2 * tNumDoubleSidedCluster );

        // vertex pair
        mIGMeshInfo->mGhostToVertexPairIndices.resize( tNumGhostVertices );

        // This used to work in leader branch,  after definition of the ghost is changed this needed to be modifed
        // since the ghost is trivial
        // int tOffset = 0;
        // if ( mIGMeshInfo->mGhostToVertexOffset.size() > 1 )
        // {
        //     tOffset = mIGMeshInfo->mGhostToVertexOffset( 1 ) - mIGMeshInfo->mGhostToVertexOffset( 0 );
        // }

        // // generate the offset data trivially
        // std::generate( mOutputMesh->mGhostLeaderFollowerVertexOffSet.begin() + 1, mOutputMesh->mGhostLeaderFollowerVertexOffSet.end(), [n = 0, &tOffset]() mutable { return n += tOffset; } );

        // leader vertex pairs
        mOutputMesh->mGhostLeaderFollowerVertexOffSet.resize( 2 * tNumGhostDblSideClusters + 1 );

        mOutputMesh->mGhostLeaderFollowerVertexOffSet( 0 ) = 0;

        size_t tGhostLeaderFollowerVertexInd = 0;

        size_t iSingleSideCounter = 1;

        // loop over double sided sets
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // make if it is a ghost
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                // clusters on the set
                Vector< Cluster const * > tClusters = iSet->get_clusters_on_set();

                // loop over the clusters
                for ( mtk::Cluster const * iCluster : tClusters )
                {
                    moris::mtk::Cluster const & tLeaderSideCluster   = iCluster->get_leader_side_cluster();
                    moris::mtk::Cluster const & tFollowerSideCluster = iCluster->get_follower_side_cluster();

                    // this is to identify the coordinates
                    tGhostLeaderFollowerVertexInd += tLeaderSideCluster.get_num_vertices_in_cluster();

                    mOutputMesh->mGhostLeaderFollowerVertexOffSet( iSingleSideCounter++ ) = tGhostLeaderFollowerVertexInd;

                    tGhostLeaderFollowerVertexInd += tFollowerSideCluster.get_num_vertices_in_cluster();

                    mOutputMesh->mGhostLeaderFollowerVertexOffSet( iSingleSideCounter++ ) = tGhostLeaderFollowerVertexInd;
                }
            }
        }

        // leader and follower side cluster vertices
        mIGMeshInfo->mGhostLeaderFollowerVertexIndices.resize( tGhostLeaderFollowerVertexInd );
    }

    // ---------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_ghost_data()
    {
        // initialize the ghost data
        this->initialize_ghost_data();

        // get the double sided sets
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // initialize the counters
        size_t      iDblSideCluster    = 0;
        size_t      iSingleSideCounter = 0;
        size_t      tLocalCoordsOffset = 0;
        moris_index iClusterIndex      = 0;

        // loop over the double sided sets
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // if it is not a ghost skip
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                // get clusters on the ghost set
                Vector< Cluster const * > tClusters = iSet->get_clusters_on_set();

                // loop over the double sided clusters
                for ( const auto& iCluster : tClusters )
                {
                    // get leader and follower side clusters
                    moris::mtk::Cluster const & tLeaderSideCluster   = iCluster->get_leader_side_cluster();
                    moris::mtk::Cluster const & tFollowerSideCluster = iCluster->get_follower_side_cluster();

                    // transfer ip cell data
                    mIGMeshInfo->mGhostLeaderToIPCellIndex( iDblSideCluster )   = tLeaderSideCluster.get_interpolation_cell_index();
                    mIGMeshInfo->mGhostFollowerToIPCellIndex( iDblSideCluster ) = tFollowerSideCluster.get_interpolation_cell_index();

                    // transfer ig cell data
                    mIGMeshInfo->mGhostLeaderToIGCellIndex( iDblSideCluster )   = tLeaderSideCluster.get_primary_cells_in_cluster()( 0 )->get_index();
                    mIGMeshInfo->mGhostFollowerToIGCellIndex( iDblSideCluster ) = tFollowerSideCluster.get_primary_cells_in_cluster()( 0 )->get_index();

                    // transfer side ordinal
                    mOutputMesh->mGhostLeaderFollowerOrd.push_back( tLeaderSideCluster.get_cell_side_ordinals()( 0 ) );
                    mOutputMesh->mGhostLeaderFollowerOrd.push_back( tFollowerSideCluster.get_cell_side_ordinals()( 0 ) );

                    // get the vertex data of each cluster
                    Vector< moris::mtk::Vertex const * > tGhostLeaderVertices   = tLeaderSideCluster.get_vertices_in_cluster();
                    Vector< moris::mtk::Vertex const * > tGhostFollowerVertices = tFollowerSideCluster.get_vertices_in_cluster();

                    // get the leader vertex pairs
                    Vector< moris::mtk::Vertex const * > const & tVertexPairs = iCluster->get_leader_vertex_pairs();

                    // transfer vertex pointers to indices
                    std::transform( tVertexPairs.cbegin(),
                            tVertexPairs.cend(),
                            mIGMeshInfo->mGhostToVertexPairIndices.begin() + mIGMeshInfo->mGhostToVertexOffset( iDblSideCluster ),
                            []( Vertex const * aVertex ) -> moris_index { return aVertex->get_index(); } );

                    // transfer leader vertex pointers to the idices
                    std::transform( tGhostLeaderVertices.begin(),
                            tGhostLeaderVertices.end(),
                            mIGMeshInfo->mGhostLeaderFollowerVertexIndices.begin() + mOutputMesh->mGhostLeaderFollowerVertexOffSet( iSingleSideCounter++ ),
                            []( Vertex const * aVertex ) -> moris_index { return aVertex->get_index(); } );

                    // transfer leader vertex pointers to the idices
                    std::transform( tGhostFollowerVertices.begin(),
                            tGhostFollowerVertices.end(),
                            mIGMeshInfo->mGhostLeaderFollowerVertexIndices.begin() + mOutputMesh->mGhostLeaderFollowerVertexOffSet( iSingleSideCounter++ ),
                            []( Vertex const * aVertex ) -> moris_index { return aVertex->get_index(); } );

                    // increment the double sided cluster count
                    iDblSideCluster++;

                    // get trivial information of the side clusters
                    mOutputMesh->mGhostLeaderFollowerIsTrivial.push_back( tLeaderSideCluster.is_trivial() );
                    mOutputMesh->mGhostLeaderFollowerIsTrivial.push_back( tFollowerSideCluster.is_trivial() );

                    // this is to store coordinates
                    if ( !tLeaderSideCluster.is_trivial() )
                    {
                        // get the coordinates
                        moris::Matrix< moris::DDRMat > tVertexCoords = tLeaderSideCluster.get_vertices_local_coordinates_wrt_interp_cell();

                        // put the matrix in the right spot
                        mOutputMesh->mSecondaryClusterVertexCoords->operator()( { tLocalCoordsOffset, tLocalCoordsOffset + tVertexCoords.n_rows() - 1 }, { 0, mIGMeshInfo->mSpatialDim - 1 } ) = tVertexCoords.matrix_data();

                        // make a map between cluster index and row number
                        mOutputMesh->mSecondaryClusterIndexToRowNumber[ iClusterIndex ] = (moris_index)tLocalCoordsOffset;

                        tLocalCoordsOffset += tVertexCoords.n_rows();
                    }

                    iClusterIndex++;

                    if ( !tFollowerSideCluster.is_trivial() )
                    {
                        // get the coordinates
                        moris::Matrix< moris::DDRMat > tVertexCoords = tLeaderSideCluster.get_vertices_local_coordinates_wrt_interp_cell();

                        // put the matrix in the right spot
                        mOutputMesh->mSecondaryClusterVertexCoords->operator()( { tLocalCoordsOffset, tLocalCoordsOffset + tVertexCoords.n_rows() - 1 }, { 0, mIGMeshInfo->mSpatialDim - 1 } ) = tVertexCoords.matrix_data();

                        // make a map between cluster index and row number
                        mOutputMesh->mSecondaryClusterIndexToRowNumber[ iClusterIndex ] = (moris_index)tLocalCoordsOffset;

                        tLocalCoordsOffset += tVertexCoords.n_rows();
                    }

                    iClusterIndex++;
                }
            }
        }
    }

    // ---------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_block_set_data()
    {
        this->initialize_block_set_data();

        // get the block sets
        Vector< moris::mtk::Block_Set* > const & tBlockSets = mInputMesh->get_block_sets();

        // set the counter
        uint iCounter = 0;

        // loop over the the block sets
        for ( auto const & iSet : tBlockSets.data() )
        {
            // get the clusters
            Vector< Cluster const * > tClustersOnBlock = iSet->get_clusters_on_set();

            // generate the indices of the cell cluster
            std::transform( tClustersOnBlock.begin(),
                    tClustersOnBlock.end(),
                    mIGMeshInfo->mBlockSetToCellClusterIndex.begin() + mIGMeshInfo->mBlockSetToCellClusterOffSet( iCounter ),
                    []( Cluster const * aCluster ) -> moris_index {
                        return aCluster->get_interpolation_cell().get_index();
                    } );

            // increment the counter
            iCounter++;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::initialize_block_set_data()
    {
        // get tge block sets from the mesh
        Vector< moris::mtk::Block_Set* > const & tBlockSets = mInputMesh->get_block_sets();

        // resize the connectivity
        mIGMeshInfo->mBlockSetToCellClusterOffSet.resize( tBlockSets.size() + 1 );

        // the first entry is always zero
        mIGMeshInfo->mBlockSetToCellClusterOffSet( 0 ) = 0;

        // initialize the counters
        uint   iCounter                = 0;
        size_t tLocalCellClusterOffSet = 0;

        // loop over block sets
        for ( auto const & iSet : tBlockSets.data() )
        {
            // get the number of the clusters
            tLocalCellClusterOffSet += iSet->get_num_clusters_on_set();

            // connectivity
            mIGMeshInfo->mBlockSetToCellClusterOffSet( iCounter + 1 ) = tLocalCellClusterOffSet;

            // increment the iCounter
            iCounter++;
        }

        // resize because of the transfrom function
        mIGMeshInfo->mBlockSetToCellClusterIndex.resize( tLocalCellClusterOffSet );

        MORIS_ASSERT( tLocalCellClusterOffSet = mIPMeshDataBase->get_num_elems(), "Number Of Cell Clusters are not equal to IP cells( cell_clusters)" );
    }

    // ----------------------------------------------------------------------------

    Integration_Mesh_DataBase_IG*
    Integration_Mesh_Editor::perform()
    {
        Tracer tTracer( "MTK", "IG DataBase", "Build" );

        // initialize te pointer data
        mOutputMesh = new Integration_Mesh_DataBase_IG();

        mOutputMesh->mIPMesh = mIPMeshDataBase;

        mOutputMesh->mIGMeshInfo = mIGMeshInfo;

        // generate raw data
        this->generate_mesh_data();

        // create the mesh objects( vertices, cells, ...)
        this->create_mesh();

        // check that input and output mesh are the same (it is active in debug)
        this->check_input_output_mesh();

        return mOutputMesh;
    }

    // ----------------------------------------------------------------------------

    uint
    Integration_Mesh_Editor::get_num_side_clusters()
    {
        return mOutputMesh->mNumSideClusters;
    }

    // ----------------------------------------------------------------------------

    moris::Memory_Map
    Integration_Mesh_Editor::get_memory_usage()
    {
        moris::Memory_Map tMemoryMap;
        return tMemoryMap;
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::free_memory()
    {
        delete mOutputMesh->mIGMeshInfo;
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_mesh()
    {
        mInputMesh->delete_visualization_sets();

        // create the vertecies
        this->create_vertices();

        // create the map ( will be used in gen)
        this->create_vertex_glb_id_to_loc_vertex_ind_map();

        // create cells
        this->create_cells();

        // create cell clusters
        this->create_cell_clusters();

        // create side clusters
        this->create_side_clusters();

        // create block set
        this->create_block_sets();

        // create side sets
        this->create_side_sets();

        // create double sided clusters
        this->create_double_sided_clusters();

        // create double sided sets
        this->create_double_sided_sets();

        // create ghost clusters ( side cluster and double sided cluster)
        this->create_ghost_clusters();

        // create the ghost sets
        this->create_ghost_sets();

        // set up cell toplogy map that will be used in the collect set function
        this->set_blockset_topology();

        // collect all sets to clean up the lists
        mOutputMesh->collect_all_sets();
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::generate_mesh_data()
    {
        this->generate_vertex_data();

        this->generate_cell_data();

        this->generate_cell_cluster_data();

        this->generate_side_cluster_data();

        this->generate_block_set_data();

        this->generate_double_sided_cluster_data();

        this->generate_ghost_data();
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_vertices()
    {
        // Allocate size for the vertices and their id and owner list
        mOutputMesh->mVertices.reserve( mIGMeshInfo->mVertices.size() );
        mOutputMesh->mVertexIdList.reserve( mIGMeshInfo->mVertices.size() );
        mOutputMesh->mVertexOwnerList.reserve( mIGMeshInfo->mVertices.size() );

        // counter for the vertices
        uint iCounter = 0;

        // loop over old vertices to transfer data
        for ( const auto& iVertex : mIGMeshInfo->mVertices )
        {
            // Assert to ensure consecutive vertices
            MORIS_ASSERT( iVertex->get_index() == (moris_index)iCounter, "Index alignment issue in vertices" );

            // constrcut the vertex and put it in the list
            mOutputMesh->mVertices.emplace_back( iCounter, mOutputMesh );

            // increase the count
            iCounter++;

            // add id and owner of the vertex
            mOutputMesh->mVertexIdList.push_back( iVertex->get_id() );

            // FIXME: owner is not implemented in the leader branch
            mOutputMesh->mVertexOwnerList.push_back( iVertex->get_index() );
        }
    }

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_cells()
    {
        // create pointers with the parent class for the connectivity
        mOutputMesh->mCellToVertices.reserve( mIGMeshInfo->mCellToVertexIndicies.size() );

        // populate the cell to vertex connectivity list
        for ( const int& iVertex : mIGMeshInfo->mCellToVertexIndicies )
        {
            mOutputMesh->mCellToVertices.push_back( &mOutputMesh->mVertices( iVertex ) );
        }

        // get number of cells and reserve enough space for cells
        uint tNumCells = mInputMesh->get_num_elems();

        // reserve enough space in the cells, id list and owner
        mOutputMesh->mCells.reserve( tNumCells );
        mOutputMesh->mCellIdList.reserve( tNumCells );
        mOutputMesh->mCellOwnerList.reserve( tNumCells );
        mOutputMesh->mCellInfoList.reserve( tNumCells );

        // loop over the old cells to create new cells
        // FIXME: in order to get the cell-info we need to find a better way to avoid two temporaries,  two optimization should be made
        for ( size_t iCell = 0; iCell < tNumCells; iCell++ )
        {
            moris::mtk::Cell& tCell = mInputMesh->get_mtk_cell( iCell );
            mOutputMesh->mCells.emplace_back( tCell,
                    iCell,
                    mOutputMesh );

            // add the id and owner of the cell
            mOutputMesh->mCellIdList.push_back( tCell.get_id() );
            mOutputMesh->mCellOwnerList.push_back( tCell.get_owner() );
            mOutputMesh->mCellInfoList.push_back( tCell.get_cell_info_sp() );
        }
    }

    //----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_cell_clusters()
    {
        // create pointers with the parent class for the connectivity of the cell cluster
        mOutputMesh->mCellClusterToPrimaryIGCell.reserve( mIGMeshInfo->mCellClusterToPrimaryIGCellIndices.size() );
        mOutputMesh->mCellClusterToVoidIGCell.reserve( mIGMeshInfo->mCellClusterToVoidIGCellIndices.size() );
        mOutputMesh->mCellClusterToVeretx.reserve( mIGMeshInfo->mCellClusterToVertexIndices.size() );

        // populate the cell cluster to primary cell connectivity list
        for ( const int& iPrimaryCell : mIGMeshInfo->mCellClusterToPrimaryIGCellIndices )
        {
            mOutputMesh->mCellClusterToPrimaryIGCell.push_back( &mOutputMesh->mCells( iPrimaryCell ) );
        }

        // populate the cell cluster to void cell connectivity list
        for ( const int& iVoidCell : mIGMeshInfo->mCellClusterToVoidIGCellIndices )
        {
            mOutputMesh->mCellClusterToVoidIGCell.push_back( &mOutputMesh->mCells( iVoidCell ) );
        }

        // populate the cell cluster to vertices connectivity list
        for ( const int& iVertex : mIGMeshInfo->mCellClusterToVertexIndices )
        {
            mOutputMesh->mCellClusterToVeretx.push_back( &mOutputMesh->mVertices( iVertex ) );
        }

        // get number of cells and reserve enough space for cells
        uint tNumCells = mIPMeshDataBase->get_num_elems();

        // reserve enough space
        mOutputMesh->mCellClusters.reserve( tNumCells );
        mOutputMesh->mCellClusterIsTrivial.reserve( tNumCells );

        // loop over the old cells to create new cells
        for ( size_t iCell = 0; iCell < tNumCells; iCell++ )
        {
            // constrcut the cell cluster
            mtk::Cell_Cluster_DataBase tCellCluster( (moris_index)iCell, mOutputMesh );

            // add the cell cluster to the list
            mOutputMesh->mCellClusters.push_back( tCellCluster );
            mOutputMesh->mCellClusterIsTrivial.push_back( mInputMesh->get_cell_cluster( iCell ).is_trivial() );
        }
    }

    //----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_side_clusters()
    {
        // create pointers with the parent class for the connectivity of the cell cluster
        mOutputMesh->mSideClusterToPrimaryIGCell.reserve( mIGMeshInfo->mSideClusterToPrimaryIGCellIndices.size() );
        // mSideClusterToVoidIGCell.reserve();
        mOutputMesh->mSideClusterToVeretx.reserve( mIGMeshInfo->mSideClusterToVertexIndices.size() );

        // populate the cell cluster to primary cell connectivity list
        for ( const int& iPrimaryCell : mIGMeshInfo->mSideClusterToPrimaryIGCellIndices )
        {
            mOutputMesh->mSideClusterToPrimaryIGCell.push_back( &mOutputMesh->mCells( iPrimaryCell ) );
        }

        // populate the cell cluster to void cell connectivity list
        // for ( const int& iVoidCell : mIGDataBase->mCellClusterToVoidIGCellIndices )
        // {
        //     mCellClusterToVoidIGCell.push_back( &mCells( iVoidCell ) );
        // }

        // populate the cell cluster to vertices connectivity list
        for ( const int& iVertex : mIGMeshInfo->mSideClusterToVertexIndices )
        {
            mOutputMesh->mSideClusterToVeretx.push_back( &mOutputMesh->mVertices( iVertex ) );
        }

        // get number of cells and reserve enough space for cells
        uint tNumSideClusters = mOutputMesh->mNumSideClusters;

        // reserve enough space for the side clusters
        // mSideClusters.reserve( tNumSideClusters );
        mOutputMesh->mSideClusters.resize( tNumSideClusters );
        mOutputMesh->mSideClusterIsTrivial.reserve( tNumSideClusters );

        // loop over the side clusters to create them
        for ( size_t iCluster = 0; iCluster < tNumSideClusters; iCluster++ )
        {
            // constrcut the side cluster
            mtk::Side_Cluster_DataBase tSideCluster( (moris_index)iCluster, mOutputMesh );

            // save the side cluster
            // mSideClusters.push_back( tSideCluster );
            mOutputMesh->mSideClusters( iCluster ) = tSideCluster;

            // call the realted the cell cluster to determine if it is trivial
            mOutputMesh->mSideClusterIsTrivial.push_back( mOutputMesh->mCellClusters( mOutputMesh->mSideClusterToIPCell( iCluster ) ).is_trivial() );
        }
    }

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_block_sets()
    {
        // call the block sets from the old mesh
        Vector< moris::mtk::Block_Set* > const & tBlockSets = mInputMesh->get_block_sets();

        // resize the new list of block sets for the new mesh
        mOutputMesh->mListOfBlocks.resize( tBlockSets.size() );

        // counter for block sets
        size_t iCounter = 0;

        // loop over the old sets and create new sets
        for ( const auto& iSet : tBlockSets.data() )
        {
            // create the list of clusters inside the block sets
            Vector< Cluster const * > aBlockSetClusters;

            // resize the cluster list on the set
            aBlockSetClusters.resize( mIGMeshInfo->mBlockSetToCellClusterOffSet( iCounter + 1 ) - mIGMeshInfo->mBlockSetToCellClusterOffSet( iCounter ) );

            // transform the indices of the cell clusters to pointers
            std::transform( mIGMeshInfo->mBlockSetToCellClusterIndex.begin() + mIGMeshInfo->mBlockSetToCellClusterOffSet( iCounter ),
                    mIGMeshInfo->mBlockSetToCellClusterIndex.begin() + mIGMeshInfo->mBlockSetToCellClusterOffSet( iCounter + 1 ),
                    aBlockSetClusters.begin(),
                    [ this ]( moris_index mClusterIndex ) -> Cluster const * {
                        return &mOutputMesh->mCellClusters( mClusterIndex );
                    } );

            // construct the new block set
            mOutputMesh->mListOfBlocks( iCounter ) = new mtk::Block_Set(
                    iSet->get_set_name(),
                    aBlockSetClusters,
                    iSet->get_set_colors(),
                    mIGMeshInfo->mSpatialDim );

            // increment the count
            iCounter++;
        }

        // communicate block sets
        mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfBlocks );
    }

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_side_sets()
    {
        // call the old side sets
        Vector< moris::mtk::Side_Set* > const & tSideSets = mInputMesh->get_side_sets();

        // resize the new side set list
        mOutputMesh->mListOfSideSets.resize( tSideSets.size() );

        // counter for the side set number
        size_t iCounter = 0;

        // loop over the  old sets to create the new sets
        for ( const auto& iSet : tSideSets.data() )
        {
            // create the list of clusters inside the side set
            Vector< Cluster const * > aSideSetClusters;

            // starting index of the side clusters
            uint tStartIndex = mIGMeshInfo->mSideSetToSideClusterOffset( iCounter );

            // resize the cluster list on the set
            aSideSetClusters.resize( mIGMeshInfo->mSideSetToSideClusterOffset( iCounter + 1 ) - mIGMeshInfo->mSideSetToSideClusterOffset( iCounter ) );

            // create a cell with the indices of the side clusters within the set
            Vector< moris_index > tIndexOfSideClusters( mIGMeshInfo->mSideSetToSideClusterOffset( iCounter + 1 ) - mIGMeshInfo->mSideSetToSideClusterOffset( iCounter ) );

            // populate the indices of the side clusters,  they are consecutive by construction in the IG data base
            std::iota( tIndexOfSideClusters.begin(), tIndexOfSideClusters.end(), tStartIndex );

            // convert indicies to pointers
            std::transform( tIndexOfSideClusters.begin(),
                    tIndexOfSideClusters.end(),
                    aSideSetClusters.begin(),
                    [ this ]( moris_index mClusterIndex ) -> Cluster const * { return &mOutputMesh->mSideClusters( mClusterIndex ); } );

            // construct the new side set
            mOutputMesh->mListOfSideSets( iCounter ) =
                    new moris::mtk::Side_Set( iSet->get_set_name(),
                            aSideSetClusters,
                            iSet->get_set_colors(),
                            mIGMeshInfo->mSpatialDim );

            // increment the set count by 1
            iCounter++;

        }    // end for: each side set

        // communicate information across all sets
        mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfSideSets );

    }    // end function: Integration_Mesh_Editor::create_side_sets()

    //-----------------------------------------------------------------------------
    void
    Integration_Mesh_Editor::create_double_sided_clusters()
    {
        // get number of cells and reserve enough space for cells
        uint tNumDblSideClusters = mIGMeshInfo->mDoubleSidedClusterToNewSideClusterIndex.size();

        mOutputMesh->mDblSideClusters.reserve( tNumDblSideClusters );

        // loop over the double sided cluster numers
        for ( size_t iDblCluster = 0; iDblCluster < tNumDblSideClusters; iDblCluster++ )
        {
            // get the leader and follower side cluster of the new mesh which will constrcut the double sided cluster
            mtk::Cluster const * tLeaderCluster   = &mOutputMesh->mSideClusters( mIGMeshInfo->mDoubleSidedClusterToNewSideClusterIndex( iDblCluster ).first );
            mtk::Cluster const * tFollowerCluster = &mOutputMesh->mSideClusters( mIGMeshInfo->mDoubleSidedClusterToNewSideClusterIndex( iDblCluster ).second );

            // initialize the vertex pairs of the double sided clusters with the right size
            Vector< mtk::Vertex const * > tVertexPair( mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( iDblCluster + 1 ) - mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( iDblCluster ) );

            // transform the vertex pairs indices to the pointers
            std::transform( mIGMeshInfo->mDoubleSidedClusterToPairedVerticesIndex.begin() + mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( iDblCluster ),
                    mIGMeshInfo->mDoubleSidedClusterToPairedVerticesIndex.begin() + mIGMeshInfo->mDoubleSidedClusterToVertexOffSet( iDblCluster + 1 ),
                    tVertexPair.begin(),
                    [ this ]( moris_index aVertexIndex ) -> Vertex const * { return &mOutputMesh->mVertices( aVertexIndex ); } );

            // constrcut the double sided cluster
            mtk::Double_Side_Cluster tDblSideCluster( tLeaderCluster, tFollowerCluster, tVertexPair );

            // save the constructed double sided cluster
            mOutputMesh->mDblSideClusters.push_back( tDblSideCluster );
        }
    }

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_double_sided_sets()
    {
        // get the old double sided sets
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // resize the new double sided sets
        mOutputMesh->mListOfDoubleSideSets.resize( tDoubleSideSets.size() );

        // counter for the DoubleSidedSets
        size_t iCounter = 0;

        // loop over the old double sided sets
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // if the double sided cluster is ghost skip it
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // create the list of clusters inside the set
            Vector< Cluster const * > aSideSetClusters;

            // starting index of the side clusters
            uint tStartIndex = mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( iCounter );

            // resize the cluster list on the set
            aSideSetClusters.resize( mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( iCounter + 1 ) - mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( iCounter ) );

            // create a cell with the indices of the side clusters within the set
            Vector< moris_index > tIndexOfSideClusters( mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( iCounter + 1 ) - mIGMeshInfo->mDobleSideSetToDoubleSidedClusterOffset( iCounter ) );

            // populate the indices of the side clusters,  they are consecutive by construction in the IG data base
            std::iota( tIndexOfSideClusters.begin(), tIndexOfSideClusters.end(), tStartIndex );

            // transform with a unary operations indicies to pointers
            std::transform( tIndexOfSideClusters.begin(),
                    tIndexOfSideClusters.end(),
                    aSideSetClusters.begin(),
                    [ this ]( moris_index mClusterIndex ) -> Cluster const * { return &mOutputMesh->mDblSideClusters( mClusterIndex ); } );

            // construct the double sided set
            mOutputMesh->mListOfDoubleSideSets( iCounter ) = new moris::mtk::Double_Side_Set( iSet->get_set_name(),
                    aSideSetClusters,
                    iSet->get_set_colors(),
                    mIGMeshInfo->mSpatialDim );

            // increment the set count by 1
            iCounter++;

        }    // end for: each double sided side set

        // communicate information across all sets
        mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfDoubleSideSets );

    }    // end function: Integration_Mesh_Editor::create_double_sided_sets()

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_ghost_clusters()
    {
        // get number of the double sided clusters
        uint tNumGhostClusters = mIGMeshInfo->mGhostLeaderToIPCellIndex.size();

        // reserve enough space for IP and IG cell list
        mOutputMesh->mGhostLeaderFollowerIPCellList.reserve( 2 * tNumGhostClusters );
        mOutputMesh->mGhostLeaderFollowerIGCellList.reserve( 2 * tNumGhostClusters );

        // all the vertices on leader and follower side
        // leader vertices, follower vertices for a ghost double sided ghost cluster
        mOutputMesh->mGhostLeaderFollowerToVertex.reserve( mIGMeshInfo->mGhostLeaderFollowerVertexIndices.size() );

        // populate the cluster to vertices,  convert indices to pointers
        for ( uint iVertex = 0; iVertex < mIGMeshInfo->mGhostLeaderFollowerVertexIndices.size(); iVertex++ )
        {
            mOutputMesh->mGhostLeaderFollowerToVertex.push_back( &mOutputMesh->mVertices( mIGMeshInfo->mGhostLeaderFollowerVertexIndices( iVertex ) ) );
        }

        // loop over the number of the clusters to add leader and follower in order
        for ( uint iCell = 0; iCell < tNumGhostClusters; iCell++ )
        {
            // add ip cell index of the leader and follower
            mOutputMesh->mGhostLeaderFollowerIPCellList.push_back( mIGMeshInfo->mGhostLeaderToIPCellIndex( iCell ) );
            mOutputMesh->mGhostLeaderFollowerIPCellList.push_back( mIGMeshInfo->mGhostFollowerToIPCellIndex( iCell ) );

            // add the ig cell pointer of the leader and then salve
            mOutputMesh->mGhostLeaderFollowerIGCellList.push_back( &mOutputMesh->get_mtk_cell( mIGMeshInfo->mGhostLeaderToIGCellIndex( iCell ) ) );
            mOutputMesh->mGhostLeaderFollowerIGCellList.push_back( &mOutputMesh->get_mtk_cell( mIGMeshInfo->mGhostFollowerToIGCellIndex( iCell ) ) );
        }

        // propely resize the ghost leader, follower and double sided cluster
        mOutputMesh->mGhostLeader.resize( tNumGhostClusters );
        mOutputMesh->mGhostFollower.resize( tNumGhostClusters );
        mOutputMesh->mGhostDblSidedSet.resize( tNumGhostClusters );

        // get the starting index of the ghost clusters
        moris_index tGhostOverallClusterIndex = mOutputMesh->mNumSideClusters;

        // loop over number of double sided clusters to create the new ghost clusters
        for ( uint iCluster = 0; iCluster < tNumGhostClusters; iCluster++ )
        {
            // constrcut the leader side cluster
            mtk::Side_Cluster_DataBase tLeaderSideCluster(
                    tGhostOverallClusterIndex++,
                    mOutputMesh );

            // save the leader cluster
            mOutputMesh->mGhostLeader( iCluster ) = tLeaderSideCluster;

            // constrcut the follower side cluster
            mtk::Side_Cluster_DataBase tFollowerSideCluster(
                    tGhostOverallClusterIndex++,
                    mOutputMesh );

            // save the follower cluster
            mOutputMesh->mGhostFollower( iCluster ) = tFollowerSideCluster;

            // initialize the vertexpair list for the double sided cluster
            Vector< mtk::Vertex const * > tVertexPair( mIGMeshInfo->mGhostToVertexOffset( iCluster + 1 ) - mIGMeshInfo->mGhostToVertexOffset( iCluster ) );

            // transform indices to pointers
            std::transform( mIGMeshInfo->mGhostToVertexPairIndices.begin() + mIGMeshInfo->mGhostToVertexOffset( iCluster ),
                    mIGMeshInfo->mGhostToVertexPairIndices.begin() + mIGMeshInfo->mGhostToVertexOffset( iCluster + 1 ),
                    tVertexPair.begin(),
                    [ this ]( moris_index aVertexIndex ) -> Vertex const * { return &mOutputMesh->mVertices( aVertexIndex ); } );

            // constrcut the double sided clusters
            mtk::Double_Side_Cluster tDblSideCluster( &mOutputMesh->mGhostLeader( iCluster ), &mOutputMesh->mGhostFollower( iCluster ), tVertexPair );

            // save the double sided clusters
            mOutputMesh->mGhostDblSidedSet( iCluster ) = tDblSideCluster;
        }
    }

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_ghost_sets()
    {
        // get the double sided sets for the mesh
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // find the start index of the double sided ghost to be added to the list
        size_t iCounterGlobal = tDoubleSideSets.size() - mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset.size() + 1;

        // local counter to identify the connectivity and other info
        // the reason for this is that ghost info is stored severalty and
        // to augment the mesh we need two indices
        size_t iCounterLocal = 0;

        // loop over the old double sided sets to create new ones
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // if it is not a ghost continue
            if ( !std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // create the list of clusters inside the set
            Vector< Cluster const * > aSideSetClusters;

            // starting index of the side clusters
            uint tStartIndex = mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( iCounterLocal );

            // resize the cluster list on the set
            aSideSetClusters.resize( mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( iCounterLocal + 1 ) - mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( iCounterLocal ) );

            // create a cell with the indices of the side clusters within the set
            Vector< moris_index > tIndexOfSideClusters( mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( iCounterLocal + 1 ) - mIGMeshInfo->mDoubleSidedGhostToSideClusterGhostOffset( iCounterLocal ) );

            // populate the indices of the side clusters,  they are consecutive by construction in the IG data base
            std::iota( tIndexOfSideClusters.begin(), tIndexOfSideClusters.end(), tStartIndex );

            // transform double side cluster indices to the pointers
            std::transform( tIndexOfSideClusters.begin(),
                    tIndexOfSideClusters.end(),
                    aSideSetClusters.begin(),
                    [ this ]( moris_index mClusterIndex ) -> Cluster const * { return &mOutputMesh->mGhostDblSidedSet( mClusterIndex ); } );

            // construct the double side set
            mOutputMesh->mListOfDoubleSideSets( iCounterGlobal ) = new moris::mtk::Double_Side_Set( iSet->get_set_name(),
                    aSideSetClusters,
                    iSet->get_set_colors(),
                    mIGMeshInfo->mSpatialDim );

            // increment the counts by 1
            iCounterGlobal++;
            iCounterLocal++;

        }    // end for: each double sided side set

        // communicate information across all sets
        mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfDoubleSideSets );

    }    // end function: Integration_Mesh_Editor::create_ghost_sets()

    //-----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::set_blockset_topology()
    {
        // loop over the block sets to constrcut a map between names and set toplogy
        for ( const auto& tSet : mOutputMesh->mListOfBlocks )
        {
            mOutputMesh->mNameToCellTopologyMap[ tSet->get_set_name() ] = mInputMesh->get_blockset_topology( tSet->get_set_name() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_vertex_glb_id_to_loc_vertex_ind_map()
    {
        // reserve the space for unordered map
        mOutputMesh->mVertexGlobalIdToLocalIndex.reserve( mOutputMesh->mVertices.size() );

        // create the vertex map used in gen based on the new vertex
        for ( uint iCounter = 0; iCounter < mOutputMesh->mVertices.size(); ++iCounter )
        {
            MORIS_ASSERT( mOutputMesh->mVertices( iCounter ).get_index() == iCounter, "Index alignment issue in vertices" );

            mOutputMesh->mVertexGlobalIdToLocalIndex[ mOutputMesh->mVertices( iCounter ).get_id() ] = iCounter;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::add_vertices( Vector< Vector< moris_index > >& aSideClusterToVertexIndices, Matrix< DDRMat > aVerticesCoords )
    {
        // get number of previous vertices from the mesh
        uint tNumPreviousVertices = mOutputMesh->mVertices.size();

        // determine the total number of the vertices
        // add vertices of each group together
        uint tNumNewVertices = 0;

        uint tNumNewVertices2 = std::accumulate( aSideClusterToVertexIndices.data().begin(),
                aSideClusterToVertexIndices.data().end(),
                tNumNewVertices,
                []( uint aContainter, Vector< moris::moris_index > aVertexGroup ) -> uint { return aContainter + aVertexGroup.size(); } );

        // sum up all the vertices to resize the data
        uint tNumAllVertices = mOutputMesh->mVertices.size() + tNumNewVertices2;
        uint tSpatialDim     = mIGMeshInfo->mSpatialDim;

        // adjust the size of the vertices cel and coordinates
        mOutputMesh->mVertices.reserve( tNumAllVertices );
        mOutputMesh->mVertexCoordinates.resize( tSpatialDim, tNumAllVertices );

        // create new vetices in the mesh with the specfied index coordinates
        for ( uint iVertex = tNumPreviousVertices; iVertex < tNumAllVertices; iVertex++ )
        {
            // create the vertex with the
            mOutputMesh->mVertices.emplace_back( iVertex, mOutputMesh );

            // copy the coordinates
            std::copy( aVerticesCoords.colptr( iVertex - tNumPreviousVertices ),
                    aVerticesCoords.colptr( iVertex - tNumPreviousVertices ) + tSpatialDim,
                    mOutputMesh->mVertexCoordinates.colptr( iVertex ) );
        }

        // create ids, this is created for parallel representation only now
        this->create_parallel_consistent_new_vertex_ids( tNumPreviousVertices );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_parallel_consistent_new_vertex_ids( moris_index aNumPreviousVertices )
    {
        // resize the vertex id list
        mOutputMesh->mVertexIdList.resize( mOutputMesh->mVertices.size() );

        // loop over the newly created vertices to assign ids,  since it is serial id = index +1 ;
        for ( uint iVertex = aNumPreviousVertices; iVertex < mOutputMesh->mVertexIdList.size(); iVertex++ )
        {
            // add the id
            mOutputMesh->mVertexIdList( iVertex ) = iVertex + 1;

            // local to global map for indices
            mOutputMesh->mVertexGlobalIdToLocalIndex[ iVertex + 1 ] = iVertex;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::create_parallel_consistent_cell_ids( moris_index aNumNewCells )
    {
        // resize the id and mesh info list
        mOutputMesh->mCellIdList.resize( mNumPreviousCells + aNumNewCells );
        mOutputMesh->mCellInfoList.resize( mNumPreviousCells + aNumNewCells );

        // cell info  and TRI3 for 2d and TET4 for the 3d case
        mtk::Cell_Info_Factory                   tCellInfoFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo;
        if ( mOutputMesh->get_spatial_dim() == 2 )
        {
            tLinearCellInfo = tCellInfoFactory.create_cell_info_sp( Geometry_Type::TRI, mtk::Interpolation_Order::LINEAR );
        }

        else
        {
            tLinearCellInfo = tCellInfoFactory.create_cell_info_sp( Geometry_Type::TET, mtk::Interpolation_Order::LINEAR );
        }

        // loop over the vertices to assign ids
        // since it is parallel id=index +1
        for ( uint iCell = mNumPreviousCells; iCell < mOutputMesh->mCellIdList.size(); iCell++ )
        {
            mOutputMesh->mCellIdList( iCell )   = iCell + 1;
            mOutputMesh->mCellInfoList( iCell ) = tLinearCellInfo;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::add_cells( Vector< Vector< moris_index > >& aSideClusterToCells,
            Vector< Vector< moris_index > >&                             aCellToVertexIndices )
    {
        // get number of previous cells
        uint tNumPreviousCells = mOutputMesh->mCells.size();

        // get number of new cells, they are stored as a cell of cell of indices
        // sum up all the new cells created
        uint tNumNewCells = 0;
        tNumNewCells      = std::accumulate( aSideClusterToCells.begin(),
                aSideClusterToCells.end(),
                tNumNewCells,
                []( uint aContainter, Vector< moris_index > aVertexGroup ) -> uint { return aContainter + aVertexGroup.size(); } );

        // obtain total number of cells after adding cells
        uint tNumAllCells = tNumPreviousCells + tNumNewCells;

        // reserve enough space in the
        mOutputMesh->mCells.reserve( tNumAllCells );

        // add the cells with the consecutive indices to the mesh
        for ( uint iCell = tNumPreviousCells; iCell < tNumAllCells; iCell++ )
        {
            mOutputMesh->mCells.emplace_back( iCell, mOutputMesh );
        }

        // cell to vertex connectivity of the new cells need to be appended to the previous one
        uint tNumNewCellToVertex = 0;
        tNumNewCellToVertex      = std::accumulate( aCellToVertexIndices.begin(),
                aCellToVertexIndices.end(),
                tNumNewCellToVertex,
                []( uint aContainter, Vector< moris::moris_index > aVertexGroup ) -> uint { return aContainter + aVertexGroup.size(); } );

        // obtain the previous cell to vertex connectivity size
        uint tNumPreviousCellToVertex    = mOutputMesh->mCellToVertices.size();
        uint tPreviousCellToVertexOffSet = mOutputMesh->mCellToVertexOffSet.size();

        // resize the cell to vertex connectivity arrays
        mOutputMesh->mCellToVertexOffSet.resize( tPreviousCellToVertexOffSet + tNumNewCells );
        mOutputMesh->mCellToVertices.resize( tNumNewCellToVertex + tNumPreviousCellToVertex );

        // initialize the offset values for the connectivity cells
        uint tOffSet       = tNumPreviousCellToVertex;
        uint tOffSetOffSet = tPreviousCellToVertexOffSet;

        // loop over the cell groups to establish cell to vertex connectivity
        for ( const auto& iCellToVertexGroup : aCellToVertexIndices )
        {
            std::transform( iCellToVertexGroup.cbegin(),
                    iCellToVertexGroup.cend(),
                    mOutputMesh->mCellToVertices.begin() + tOffSet,
                    [ this ]( moris_index aVertexIndex ) -> moris::mtk::Vertex* { return &mOutputMesh->mVertices( aVertexIndex + mNumPreviousVertices ); } );

            tOffSet += iCellToVertexGroup.size();
        }

        // loop over the cell groups to establish cell to vertex offset based on the connectivity based on the dimension
        if ( mOutputMesh->get_spatial_dim() == 2 )
        {
            for ( const auto& iCellToVertexGroup : aCellToVertexIndices )
            {
                // populate the offset values
                mOutputMesh->mCellToVertexOffSet( tOffSetOffSet ) = mOutputMesh->mCellToVertexOffSet( tOffSetOffSet - 1 ) + iCellToVertexGroup.size();

                tOffSetOffSet++;
            }
        }
        else
        {
            // loop over the data to populate the offsets
            for ( const auto& iCellToVertexGroup : aCellToVertexIndices )
            {
                for ( size_t i = 0; i < iCellToVertexGroup.size() / 4; i++ )
                {
                    // populate the offset values
                    mOutputMesh->mCellToVertexOffSet( tOffSetOffSet ) = mOutputMesh->mCellToVertexOffSet( tOffSetOffSet - 1 ) + 4;

                    tOffSetOffSet++;
                }
            }
        }

        // create cell id and cell info
        this->create_parallel_consistent_cell_ids( tNumNewCells );
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::add_side_clusters( Vector< Vector< moris_index > >& aSideClusterToCells,
            Vector< moris_index >&                                               aSideClusterToIPCell,
            Matrix< DDRMat >&                                                    aVertexParametricCoords,
            Vector< Vector< moris_index > >&                                     aSideClusterToVertexIndices )
    {
        // get number of parametric coordinates
        uint mNumParamCoords = aVertexParametricCoords.n_rows();

        // uint mNumPreviousCells = 0;
        //  cell cluster to ig cell connectivity of the existing mesh
        uint tPerviousSideClusterToIGCell = mOutputMesh->mSideClusterToPrimaryIGCell.size();

        // get total number of the cells connected
        uint tNewSideClusterToIGCell = 0;
        tNewSideClusterToIGCell      = std::accumulate( aSideClusterToCells.begin(),
                aSideClusterToCells.end(),
                tNewSideClusterToIGCell,
                []( uint aContainter, Vector< moris::moris_index > aCellGroup ) -> uint { return aContainter + aCellGroup.size(); } );

        // resize the side cluster to ig cell and their side ordinal connectivity
        mOutputMesh->mSideClusterToPrimaryIGCell.resize( tPerviousSideClusterToIGCell + tNewSideClusterToIGCell );
        mOutputMesh->mSideClusterToPrimaryIGCellSideOrd.resize( tPerviousSideClusterToIGCell + tNewSideClusterToIGCell );

        // offset data for the side cluster to ig cell connectivity
        uint tPerviousSideClusterToIGCellOffset = mOutputMesh->mSideClusterToPrimaryIGCellOffset.size();
        uint tNewSideClusterToIGCellOffset      = aSideClusterToCells.size();

        // resize the side cluster to ig cell connectivity
        mOutputMesh->mSideClusterToPrimaryIGCellOffset.resize( tPerviousSideClusterToIGCellOffset + tNewSideClusterToIGCellOffset );

        // offset info for the connectivity and
        uint tOffSetForOffSet = tPerviousSideClusterToIGCellOffset;
        uint tOffSet          = tPerviousSideClusterToIGCell;

        // loop over the indices and conver them to pointers
        for ( const auto& iCellGroupInSideCluster : aSideClusterToCells )
        {
            std::transform( iCellGroupInSideCluster.cbegin(),
                    iCellGroupInSideCluster.cend(),
                    mOutputMesh->mSideClusterToPrimaryIGCell.begin() + tOffSet,
                    [ this ]( moris_index aCellIndex ) -> moris::mtk::Cell* { return &mOutputMesh->mCells( aCellIndex + mNumPreviousCells ); } );

            // fill out the side ordinals,  by construction they are either zero or 3
            if ( mOutputMesh->get_spatial_dim() == 2 )
            {
                std::fill_n( mOutputMesh->mSideClusterToPrimaryIGCellSideOrd.begin() + tOffSet, iCellGroupInSideCluster.size(), 0 );
            }
            else
            {
                std::fill_n( mOutputMesh->mSideClusterToPrimaryIGCellSideOrd.begin() + tOffSet, iCellGroupInSideCluster.size(), 3 );
            }

            // increment the offset value
            tOffSet += iCellGroupInSideCluster.size();

            // fill out the offset values
            mOutputMesh->mSideClusterToPrimaryIGCellOffset( tOffSetForOffSet ) =
                    mOutputMesh->mSideClusterToPrimaryIGCellOffset( tOffSetForOffSet - 1 ) + iCellGroupInSideCluster.size();

            tOffSetForOffSet++;
        }

        // reserve enough space for side cluster to ip cell
        mOutputMesh->mSideClusterToIPCell.reserve( mOutputMesh->mSideClusterToIPCell.size() + aSideClusterToIPCell.size() );
        // update the side cluter to ip cell connectivity ( it is 1 to 1 ), so there is not offset
        mOutputMesh->mSideClusterToIPCell.append( aSideClusterToIPCell );

        // obtain number of all side clusters in the mesh
        uint tNumAllClusters = mNumPreviousSideCluster + aSideClusterToIPCell.size();

        // update the size of side clusters
        mOutputMesh->mSideClusters.resize( tNumAllClusters );

        // update the side cluster trivial info
        mOutputMesh->mSideClusterIsTrivial.reserve( tNumAllClusters );

        // loop over the new side clusters to constrcut them
        for ( size_t iCluster = mNumPreviousSideCluster; iCluster < tNumAllClusters; iCluster++ )
        {
            // constrcut the side cluster
            mtk::Side_Cluster_DataBase tSideCluster( (moris_index)iCluster, mOutputMesh );

            // save the side cluster
            mOutputMesh->mSideClusters( iCluster ) = tSideCluster;

            // side cluster is not trivial by construction
            mOutputMesh->mSideClusterIsTrivial.push_back( false );
        }

        // vertices to side cluster connectivity
        uint tPreviousSideClusterToVertex = mOutputMesh->mSideClusterToVeretx.size();
        uint tNewSideClusterToVertex      = aVertexParametricCoords.n_rows();

        // resize the side cluster to vertex
        mOutputMesh->mSideClusterToVeretx.resize( tPreviousSideClusterToVertex + tNewSideClusterToVertex );

        // resize the side cluster to vertex offset
        mOutputMesh->mSideClusterToVertexOffSet.resize( mNumPreviousSideCluster + aSideClusterToCells.size() + 1 );

        // initialize the counters in order to find the right index to put into
        tOffSet            = tPreviousSideClusterToVertex;
        uint tOffSetOffSet = mNumPreviousSideCluster + 1;

        // loop over the side cluster vertex groups in order to establish the connectivity and offset data
        for ( const auto& iVertexGroupInSideCluster : aSideClusterToVertexIndices )
        {

            // we connect all the vertices except the top node of the tet and tri
            std::transform( iVertexGroupInSideCluster.cbegin(),
                    iVertexGroupInSideCluster.cend() - 1,
                    mOutputMesh->mSideClusterToVeretx.begin() + tOffSet,
                    [ this ]( moris_index aVertexIndex ) -> moris::mtk::Vertex* { return &mOutputMesh->mVertices( aVertexIndex + mNumPreviousVertices ); } );

            // populate the veretx offset
            mOutputMesh->mSideClusterToVertexOffSet( tOffSetOffSet ) = mOutputMesh->mSideClusterToVertexOffSet( tOffSetOffSet - 1 ) + iVertexGroupInSideCluster.size() - 1;

            // update the offset
            tOffSetOffSet++;
            tOffSet += iVertexGroupInSideCluster.size() - 1;
        }

        // starting row number in order to add it to the map
        uint tRowNumber = mOutputMesh->mCellClusterVertexCoords->n_rows();

        // resize the local vertex coordinates based on the
        mOutputMesh->mCellClusterVertexCoords->resize( mOutputMesh->mCellClusterVertexCoords->n_rows() + mNumParamCoords, mIGMeshInfo->mSpatialDim );

        // combine the local coordinate matrix
        mOutputMesh->mCellClusterVertexCoords->operator()( { mOutputMesh->mCellClusterVertexCoords->n_rows() - mNumParamCoords,
                                                                   mOutputMesh->mCellClusterVertexCoords->n_rows() - 1 },
                { 0, mIGMeshInfo->mSpatialDim - 1 } ) = aVertexParametricCoords.matrix_data();

        // reserve enough space on the map
        mOutputMesh->mSideClusterIndexToRowNumber.reserve( aSideClusterToVertexIndices.size() );

        // side cluster starting index
        uint iCluster = 0;

        // loop over the vertices in side cluster in order to correspond the row number to cluster index
        for ( const auto& iVertexGroup : aSideClusterToVertexIndices )
        {
            mOutputMesh->mSideClusterIndexToRowNumber[ iCluster + mNumPreviousSideCluster ] = tRowNumber;

            // increment  the counts
            tRowNumber += iVertexGroup.size() - 1;
            iCluster++;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::add_double_sided_clusters( uint mNumDblSideCluster, Vector< Vector< moris_index > >& aSideClusterToVertexIndices )
    {
        // resize the dbl sided cluster
        mOutputMesh->mDblSideClusters.reserve( mNumDblSideCluster + mOutputMesh->mDblSideClusters.size() );

        // loop over the double sided clusters
        for ( uint iDblSideCluster = 0; iDblSideCluster < mNumDblSideCluster; iDblSideCluster++ )
        {
            // initialize the vertex pairs of the double sided clusters with the right size
            Vector< mtk::Vertex const * > tVertexPair( aSideClusterToVertexIndices( iDblSideCluster * 2 + 1 ).size() - 1 );

            // transform the vertex pairs indices to the pointers
            std::transform( aSideClusterToVertexIndices( iDblSideCluster * 2 + 1 ).begin(),
                    aSideClusterToVertexIndices( iDblSideCluster * 2 + 1 ).end() - 1,
                    tVertexPair.begin(),
                    [ this ]( moris_index aVertexIndex ) -> Vertex const * { return &mOutputMesh->mVertices( aVertexIndex + mNumPreviousVertices ); } );

            // get the leader side and follower side clusters
            mtk::Cluster* tLeaderCluster   = &mOutputMesh->mSideClusters( 2 * iDblSideCluster + mNumPreviousSideCluster );
            mtk::Cluster* tFollowerCluster = &mOutputMesh->mSideClusters( 2 * iDblSideCluster + 1 + mNumPreviousSideCluster );

            // constrcut the double sided cluster
            mtk::Double_Side_Cluster tDblSideCluster( tLeaderCluster, tFollowerCluster, tVertexPair );

            // save the constructed double sided cluster
            mOutputMesh->mDblSideClusters.push_back( tDblSideCluster );
        }
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::add_double_sided_set( Vector< moris_index >& aDoubleSidedClustersIndex, uint aNumGeometry )
    {
        // get number of geometries from the geometry engine
        Vector< Vector< moris_index > > tDoubleSideSetIndices( aNumGeometry * aNumGeometry );

        // constrcut the phase to phase interaction table
        Matrix< IndexMat > tPhaseInteractionTable( aNumGeometry, aNumGeometry );

        for ( uint iRow = 0; iRow < aNumGeometry; iRow++ )
        {
            for ( uint iCol = 0; iCol < aNumGeometry; iCol++ )
            {
                tPhaseInteractionTable( iRow, iCol ) = iRow * aNumGeometry + iCol;
            }
        }

        // counter for the loop
        uint iCounter = 0;

        // loop over the double sided cluster indices and assign the correct phase combination
        for ( const auto& iIndex : aDoubleSidedClustersIndex )
        {
            tDoubleSideSetIndices( iIndex ).push_back( iCounter );

            iCounter++;
        }

        // obtain number of the double sided sets
        uint tNumNewDoubleSideSet = 0;

        // loop over and if one is not empty then there is a side set
        for ( const auto& iDoubleSideSetIndices : tDoubleSideSetIndices )
        {
            if ( iDoubleSideSetIndices.size() )
            {
                tNumNewDoubleSideSet++;
            }
        }

        // resize number of the double sided set
        mOutputMesh->mListOfDoubleSideSets.resize( mNumPreviousDoubleSideSet + tNumNewDoubleSideSet );

        // counter for the double sided set
        iCounter = 0;

        // loop over the phases to create the double sided set
        for ( uint iFirstPhase = 0; iFirstPhase < aNumGeometry; iFirstPhase++ )
        {
            for ( uint iSecondPhase = 0; iSecondPhase < aNumGeometry; iSecondPhase++ )
            {
                uint tPhaseToPhaseIndex = tPhaseInteractionTable( iFirstPhase, iSecondPhase );

                // check if double sided set is empty
                if ( tDoubleSideSetIndices( tPhaseToPhaseIndex ).size() > 0 )
                {
                    // name the double sided set
                    std::string tDoubleSideSetName = "P" + std::to_string( iFirstPhase ) + std::to_string( iSecondPhase );

                    // colors of the double sided set,  needs to be modified
                    moris::Matrix< moris::IndexMat > tColors = { { 0 } };

                    // create the list of clusters inside the set
                    Vector< Cluster const * > tDoubleSideSetClusters( tDoubleSideSetIndices( tPhaseToPhaseIndex ).size() );

                    // transform with a unary operations indicies to pointers
                    std::transform( tDoubleSideSetIndices( tPhaseToPhaseIndex ).begin(),
                            tDoubleSideSetIndices( tPhaseToPhaseIndex ).end(),
                            tDoubleSideSetClusters.begin(),
                            [ this ]( moris_index mClusterIndex ) -> Cluster const * { return &mOutputMesh->mDblSideClusters( mClusterIndex + mNumPreviousDblSideCluster ); } );

                    // construct the double sided set
                    mOutputMesh->mListOfDoubleSideSets( iCounter + mNumPreviousDoubleSideSet ) =
                            new moris::mtk::Double_Side_Set( tDoubleSideSetName,
                                    tDoubleSideSetClusters,
                                    tColors,
                                    mIGMeshInfo->mSpatialDim );
                    iCounter++;
                }

            }    // end for: each follower phase

        }    // end for: each leader phase candidate

        // communicate information across all sets
        mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfDoubleSideSets );

    }    // end function: Integration_Mesh_Editor::add_double_sided_set()

    //------------------------------------------------------------------------------------------------------------

    // end function: Integration_Mesh_Editor::merge_meshes()

    //------------------------------------------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::recreate_side_sets()
    {
        // call the old side sets
        Vector< moris::mtk::Side_Set* > tSideSets = mOutputMesh->get_side_sets();

        // resize the new side set list
        mOutputMesh->mListOfSideSets.resize( tSideSets.size() );

        // counter for the side set number
        size_t iCounter = 0;

        // loop over the  old sets to create the new sets
        for ( const auto& iSet : tSideSets.data() )
        {
            // create the list of clusters inside the side set
            Vector< Cluster const * > aSideSetClusters;

            // starting index of the side clusters
            uint tStartIndex = mIGMeshInfo->mSideSetToSideClusterOffset( iCounter );

            // resize the cluster list on the set
            aSideSetClusters.resize( mIGMeshInfo->mSideSetToSideClusterOffset( iCounter + 1 ) - mIGMeshInfo->mSideSetToSideClusterOffset( iCounter ) );

            // create a cell with the indices of the side clusters within the set
            Vector< moris_index > tIndexOfSideClusters( mIGMeshInfo->mSideSetToSideClusterOffset( iCounter + 1 ) - mIGMeshInfo->mSideSetToSideClusterOffset( iCounter ) );

            // populate the indices of the side clusters,  they are consecutive by construction in the IG data base
            std::iota( tIndexOfSideClusters.begin(), tIndexOfSideClusters.end(), tStartIndex );

            // convert indicies to pointers
            std::transform( tIndexOfSideClusters.begin(),
                    tIndexOfSideClusters.end(),
                    aSideSetClusters.begin(),
                    [ this ]( moris_index mClusterIndex ) -> Cluster const * { return &mOutputMesh->mSideClusters( mClusterIndex ); } );

            // construct the new side set
            mOutputMesh->mListOfSideSets( iCounter ) =
                    new moris::mtk::Side_Set(
                            iSet->get_set_name(),
                            aSideSetClusters,
                            iSet->get_set_colors(),
                            mIGMeshInfo->mSpatialDim );

            // increment the set count by 1
            iCounter++;

        }    // end for: each side set

        // communicate information across all sets
        mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfSideSets );

        // to delete the old pointers
        for ( const auto& iSet : tSideSets.data() )
        {
            delete iSet;
        }

    }    // end function: Integration_Mesh_Editor::recreate_side_sets()

    //------------------------------------------------------------------------------------------------------------

    void
    Integration_Mesh_Editor::check_input_output_mesh()
    {
        if ( mCheckMesh )
        {
            // checks to confirm that the new mesh is equivalent to the old mesh
            MORIS_ASSERT( this->check_vertices(), "Vertices in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_cells(), "Cells in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_cell_clusters(), "Cell clusters in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_side_clusters(), "side clusters in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_block_sets(), "Block sets in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_side_sets(), "Side Sets in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_double_sided_clusters(), "Double sided clusters in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_ghost_clusters(), "Ghost side clusters and double sided clusters in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_double_sided_sets(), "Double sided sets in the old and new mesh are not identical" );

            MORIS_ASSERT( this->check_maps(), "Maps in the old and new mesh are not identical" );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_vertices()
    {

        bool tOutput = true;

        // check to see if the vertices are stored in consecutive manner and have the same id and indices
        bool tVertexIdAndIndexEqual = std::equal( mOutputMesh->mVertices.begin(),
                mOutputMesh->mVertices.end(),
                mIGMeshInfo->mVertices.begin(),
                []( Vertex_DataBase a, mtk::Vertex const * b ) -> bool { return a.get_id() == b->get_id() && a.get_index() == b->get_index(); } );

        // check if old vertices and new vertices have the same coords
        bool tEqualCoords = std::equal( mOutputMesh->mVertices.begin(),
                mOutputMesh->mVertices.end(),
                mIGMeshInfo->mVertices.begin(),
                []( Vertex_DataBase a, mtk::Vertex const * b ) -> bool {
                    return std::equal( a.get_coords().begin(), a.get_coords().end(), b->get_coords().begin() );
                } );

        // combine the output
        tOutput = tVertexIdAndIndexEqual && tEqualCoords;

        if ( !tOutput )
        {
            // Log n error message
            MORIS_LOG_ERROR( "Vertices are inconsistent, Equal vertex id and index: %d,Equal coordinates %d ",
                    tVertexIdAndIndexEqual,
                    tEqualCoords );

            return tOutput;
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_cells()
    {
        bool tOutput = true;

        // get numer of cells
        uint tNumCells = mInputMesh->get_num_elems();

        // loop over the cells
        for ( size_t iCell = 0; iCell < tNumCells; iCell++ )
        {
            // get the old and new cell
            mtk::Cell const & tCellOldMesh = mInputMesh->get_mtk_cell( iCell );
            mtk::Cell const & tCellNewMesh = mOutputMesh->get_mtk_cell( iCell );

            // check if they have the same ID
            bool tSameCellId = tCellOldMesh.get_id() == tCellNewMesh.get_id();

            // get vertices of the cells
            Vector< Vertex* > tOldVertices = tCellOldMesh.get_vertex_pointers();
            Vector< Vertex* > tNewVertices = tCellNewMesh.get_vertex_pointers();

            // check if the vertices of the cell have the same id and index
            bool tVertexIdAndIndexEqual = std::equal( tOldVertices.begin(),
                    tOldVertices.end(),
                    tNewVertices.begin(),
                    []( Vertex* a, Vertex* b ) -> bool { return a->get_id() == b->get_id() && a->get_index() == b->get_index(); } );

            // overwrite the output with thr
            tOutput = tSameCellId && tVertexIdAndIndexEqual && tOldVertices.size() == tNewVertices.size();

            if ( !tOutput )
            {
                // Log n error message
                MORIS_LOG_ERROR( "Cell number %zu of the mesh, is inconsistent, Equal Cell ID: %d,Equal Vertices %d ",
                        iCell,
                        tSameCellId,
                        tVertexIdAndIndexEqual );
                return tOutput;
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_cell_clusters()
    {
        // initialize the output
        bool tOutput = true;

        // get number of the cell clusters/IP cells
        uint tNumCells = mIPMeshDataBase->get_num_elems();

        // loop over cells ( new and old)
        for ( size_t iCluster = 0; iCluster < tNumCells; iCluster++ )
        {
            // extract the clusters to be examined further
            mtk::Cell_Cluster const & tCellClusterOldMesh = mInputMesh->get_cell_cluster( iCluster );
            mtk::Cell_Cluster const & tCellClusterNewMesh = mOutputMesh->get_cell_cluster( iCluster );

            // extract member data of clusters ( ip cell, primary ig cells, void ig cells, vertices, and local coord matrix)
            moris::mtk::Cell const & tOldIPCell = tCellClusterOldMesh.get_interpolation_cell();
            moris::mtk::Cell const & tNewIPCell = tCellClusterNewMesh.get_interpolation_cell();

            // extract the primary cells of the clusters
            Vector< moris::mtk::Cell const * > const & tOldPrimaryCells = tCellClusterOldMesh.get_primary_cells_in_cluster();
            Vector< moris::mtk::Cell const * > const & tNewPrimaryCells = tCellClusterNewMesh.get_primary_cells_in_cluster();

            // extract the void cells of the clusters
            Vector< moris::mtk::Cell const * > const & tOldVoidCells = tCellClusterOldMesh.get_void_cells_in_cluster();
            Vector< moris::mtk::Cell const * > const & tNewVoidCells = tCellClusterNewMesh.get_void_cells_in_cluster();

            // extract the vertices of the clusters
            Vector< moris::mtk::Vertex const * > tOldVertices = tCellClusterOldMesh.get_vertices_in_cluster();
            Vector< moris::mtk::Vertex const * > tNewVertices = tCellClusterNewMesh.get_vertices_in_cluster();

            // compute the cell cluster measures
            moris::real tMeasureLeaderOld = tCellClusterOldMesh.compute_cluster_cell_measure();
            moris::real tMeasureLeaderNew = tCellClusterNewMesh.compute_cluster_cell_measure();

            // compare the vertex coordinates
            Matrix< DDRMat > tOldLocalCoords = tCellClusterOldMesh.get_vertices_local_coordinates_wrt_interp_cell();
            Matrix< DDRMat > tNewLocalCoords = tCellClusterNewMesh.get_vertices_local_coordinates_wrt_interp_cell();

            // equal measures
            bool tEqualClusterMeasure = tMeasureLeaderNew == tMeasureLeaderOld;

            // compare the primary ig cell id and index
            bool tPrimaryCellIdAndIndexEqual = std::equal( tOldPrimaryCells.cbegin(),
                    tOldPrimaryCells.cend(),
                    tNewPrimaryCells.cbegin(),
                    []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                        return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                    } );

            // compare the void ig cell id and index
            bool tVoidCellIdAndIndexEqual = std::equal( tOldVoidCells.cbegin(),
                    tOldVoidCells.cend(),
                    tNewVoidCells.cbegin(),
                    []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                        return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                    } );

            // compare the vertices id and index
            bool tVertexIdAndIndexEqual = std::equal( tOldVertices.begin(),
                    tOldVertices.end(),
                    tNewVertices.begin(),
                    []( Vertex const * a, Vertex const * b ) -> bool {
                        return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                    } );

            // compare the ip cell
            bool tEqualIPCell = tOldIPCell.get_id() == tNewIPCell.get_id();

            // compare the vertex coordinates
            bool tEqualLocalCoords = std::equal( tOldLocalCoords.begin(), tOldLocalCoords.end(), tNewLocalCoords.begin() );

            // compare the individual coordinate of the cluster
            if ( !tCellClusterNewMesh.is_trivial() )
            {
                for ( uint i = 0; i < tOldVertices.size(); i++ )
                {
                    Matrix< DDRMat > tOldVertLocalCoords = tCellClusterOldMesh.get_vertex_local_coordinate_wrt_interp_cell( tOldVertices( i ) );
                    Matrix< DDRMat > tNewVertLocalCoords = tCellClusterNewMesh.get_vertex_local_coordinate_wrt_interp_cell( tNewVertices( i ) );

                    bool tEqualVertLocalCoords = std::equal( tOldVertLocalCoords.begin(), tOldVertLocalCoords.end(), tNewVertLocalCoords.begin() );

                    MORIS_ERROR( tEqualVertLocalCoords, "Local Corrds of the individual vertices are not the same" );
                }
            }

            // combine the output info
            tOutput = tEqualClusterMeasure
                   && tEqualIPCell
                   && tVertexIdAndIndexEqual
                   && tPrimaryCellIdAndIndexEqual
                   && tVoidCellIdAndIndexEqual
                   && tEqualLocalCoords
                   && tOldVertices.size() == tNewVertices.size();

            if ( !tOutput )
            {
                MORIS_LOG_ERROR( "Cell cluster number %zu of the mesh, is inconsistent, equal cluster measure %d , equal ip cell %d, equal vertices %d, primary cells %d, void cells %d, local coordinates %d",
                        iCluster,
                        tEqualClusterMeasure,
                        tEqualIPCell,
                        tVertexIdAndIndexEqual,
                        tPrimaryCellIdAndIndexEqual,
                        tVoidCellIdAndIndexEqual,
                        tEqualLocalCoords );
                return tOutput;
            }
        }
        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_side_clusters()
    {
        // check the side clusters
        bool tOutput = true;

        // get the side sets
        Vector< moris::mtk::Side_Set* > const & tSideSets = mInputMesh->get_side_sets();

        // counter for clusters
        uint iCounter = 0;

        // loop over sets
        for ( const auto& iSet : tSideSets.data() )
        {

            // if ghost, skip
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) || std::strstr( iSet->get_set_name().c_str(), "Ghost" ) )
            {
                continue;
            }

            // get the cluster of the old mesh
            Vector< Cluster const * > tSideClusters = iSet->get_clusters_on_set();

            // loop over the clusters
            for ( const auto& iCluster : tSideClusters )
            {
                // get old side cluster data
                moris::Matrix< moris::IndexMat >           tSideOrdinalOld  = iCluster->get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellsOld = iCluster->get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexOld       = iCluster->get_vertices_in_cluster();
                uint                                       tIPCellIndexOld  = iCluster->get_interpolation_cell_index();
                Matrix< DDRMat >                           tOldLocalCoords  = iCluster->get_vertices_local_coordinates_wrt_interp_cell();

                // get new side cluster data
                moris::Matrix< moris::IndexMat >           tSideOrdinalNew = mOutputMesh->mSideClusters( iCounter ).get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellNew = mOutputMesh->mSideClusters( iCounter ).get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexNew      = mOutputMesh->mSideClusters( iCounter ).get_vertices_in_cluster();
                uint                                       tIPCellIndexNew = mOutputMesh->mSideClusters( iCounter ).get_interpolation_cell_index();
                Matrix< DDRMat >                           tNewLocalCoords = mOutputMesh->mSideClusters( iCounter ).get_vertices_local_coordinates_wrt_interp_cell();

                // compare the vertices id and index
                bool tVertexIdAndIndexEqual = std::equal( tVertexOld.begin(),
                        tVertexOld.end(),
                        tVertexNew.begin(),
                        []( Vertex const * a, Vertex const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                // compare the local coordinates
                bool tEqualSideOrdinals = std::equal( tSideOrdinalOld.begin(), tSideOrdinalOld.end(), tSideOrdinalNew.begin() );

                // compare the local coordinates
                bool tEqualLocalCoords = std::equal( tOldLocalCoords.begin(), tOldLocalCoords.end(), tNewLocalCoords.begin() );

                // compare the primary ig cell id and index
                bool tPrimaryCellIdAndIndexEqual = std::equal( tPrimaryCellsOld.cbegin(),
                        tPrimaryCellsOld.cend(),
                        tPrimaryCellNew.cbegin(),
                        []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                bool tEqualIPCell = tIPCellIndexOld == tIPCellIndexNew;

                // test the functionality vertex cluster index
                moris_index tVertexOldInd;
                moris_index tVertexNewInd;

                // loop over the vertices
                for ( uint i = 0; i < tVertexOld.size(); i++ )
                {
                    // tVertexNewInd = mOutputMesh->mSideClusters( iCounter ).get_vertex_cluster_index( tVertexNew( i ) );
                    // tVertexOldInd = iCluster->get_vertex_cluster_index( tVertexOld( i ) );
                    tVertexNewInd = 1;
                    tVertexOldInd = 1;

                    MORIS_ERROR( tVertexOldInd == tVertexNewInd, "Vertices do not match!" );
                }

                // individual vertex coordinates
                bool tEqualVertLocalCoords = true;

                // loop over the individual vertices
                for ( uint i = 0; i < tVertexOld.size(); i++ )
                {
                    Matrix< DDRMat > tOldVertLocalCoords = iCluster->get_vertex_local_coordinate_wrt_interp_cell( tVertexOld( i ) );
                    Matrix< DDRMat > tNewVertLocalCoords = mOutputMesh->mSideClusters( iCounter ).get_vertex_local_coordinate_wrt_interp_cell( tVertexNew( i ) );

                    tEqualVertLocalCoords = std::equal( tOldVertLocalCoords.begin(), tOldVertLocalCoords.end(), tNewVertLocalCoords.begin() );
                    MORIS_ERROR( tEqualVertLocalCoords, "Local Corrds are not the same" );
                }

                // combine the output
                tOutput = tEqualIPCell
                       && tVertexIdAndIndexEqual
                       && tPrimaryCellIdAndIndexEqual
                       && tEqualSideOrdinals
                       && tEqualLocalCoords;

                // if output is false return
                if ( !tOutput )
                {
                    MORIS_LOG_ERROR(
                            "Side cluster number %u is inconsistent; "
                            "tEqualIPCell %d \n"
                            "tVertexIdAndIndexEqual %d \n"
                            "tPrimaryCellIdAndIndexEqual %d \n"
                            "tEqualSideOrdinals %d \n"
                            "tEqualLocalCoords %d \n",
                            iCounter,
                            tEqualIPCell,
                            tVertexIdAndIndexEqual,
                            tPrimaryCellIdAndIndexEqual,
                            tEqualSideOrdinals,
                            tEqualLocalCoords );

                    return tOutput;
                }

                // increment the counter by 1
                iCounter++;
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_block_sets()
    {
        // initialize the output
        bool tOutput = true;

        // get the old block sets
        Vector< moris::mtk::Block_Set* > const & tBlocks = mInputMesh->get_block_sets();

        // loop over the blocks
        for ( uint iSet = 0; iSet < mOutputMesh->mListOfBlocks.size(); iSet++ )
        {
            // clusters on the set from the mesh
            Vector< Cluster const * > tClustersOldMesh = tBlocks( iSet )->get_clusters_on_set();
            Vector< Cluster const * > tClustersNewMesh = mOutputMesh->mListOfBlocks( iSet )->get_clusters_on_set();

            // skip check for sets only used for visualization purposes
            if ( !std::strstr( tBlocks( iSet )->get_set_name().c_str(), "Vis" ) )
            {
                // loop over the cluster get their data compare
                for ( uint iCluster = 0; iCluster < tClustersNewMesh.size(); iCluster++ )
                {
                    // member data for the old cluster mesh
                    Vector< moris::mtk::Cell const * > const & tPrimaryCellsOld  = tClustersOldMesh( iCluster )->get_primary_cells_in_cluster();
                    Vector< moris::mtk::Vertex const * >       tVertexOld        = tClustersOldMesh( iCluster )->get_vertices_in_cluster();
                    uint                                       tIPCellIndexOld   = tClustersOldMesh( iCluster )->get_interpolation_cell_index();
                    moris::real                                tMeasureLeaderOld = tClustersOldMesh( iCluster )->compute_cluster_cell_measure();

                    // member data for the new cluster mesh
                    Vector< moris::mtk::Cell const * > const & tPrimaryCellNew   = tClustersNewMesh( iCluster )->get_primary_cells_in_cluster();
                    Vector< moris::mtk::Vertex const * >       tVertexNew        = tClustersNewMesh( iCluster )->get_vertices_in_cluster();
                    uint                                       tIPCellIndexNew   = tClustersNewMesh( iCluster )->get_interpolation_cell_index();
                    moris::real                                tMeasureLeaderNew = tClustersNewMesh( iCluster )->compute_cluster_cell_measure();

                    // equal measures
                    bool tEqualClusterMeasure = tMeasureLeaderNew == tMeasureLeaderOld;

                    // compare the primary ig cell id and index
                    bool tPrimaryCellIdAndIndexEqual = std::equal( tPrimaryCellsOld.cbegin(),
                            tPrimaryCellsOld.cend(),
                            tPrimaryCellNew.cbegin(),
                            []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                                return a->get_index() == b->get_index();
                            } );

                    // compare the vertices id and index
                    bool tVertexIdAndIndexEqual = std::equal( tVertexOld.begin(),
                            tVertexOld.end(),
                            tVertexNew.begin(),
                            []( Vertex const * a, Vertex const * b ) -> bool {
                                return a->get_index() == b->get_index();
                            } );

                    // equal ip cell
                    bool tEqualIPCell = tIPCellIndexOld == tIPCellIndexNew;

                    // combine the outputs
                    tOutput = tEqualClusterMeasure
                           && tEqualIPCell
                           && tVertexIdAndIndexEqual
                           && tPrimaryCellIdAndIndexEqual;

                    if ( !tOutput )
                    {
                        MORIS_LOG_ERROR( "cell cluster number %u of the mesh on the block set, is inconsistent, tEqualClusterMeasure %d,tEqualIPCell %d, Vertices %d, Primary cells %d ",
                                iCluster,
                                tEqualClusterMeasure,
                                tEqualIPCell,
                                tVertexIdAndIndexEqual,
                                tPrimaryCellIdAndIndexEqual );

                        return tOutput;
                    }
                }
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_side_sets()
    {
        // initialize the output
        bool tOutput = true;

        // get the old side sets
        Vector< moris::mtk::Side_Set* > const & tSideSets = mInputMesh->get_side_sets();

        for ( uint iSet = 0; iSet < mOutputMesh->mListOfSideSets.size(); iSet++ )
        {

            // if it is a ghost don't check yet
            if ( std::strstr( mOutputMesh->mListOfSideSets( iSet )->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // obtain old and new clusters
            Vector< Cluster const * > tSideClustersOldMesh = tSideSets( iSet )->get_clusters_on_set();
            Vector< Cluster const * > tSideClustersNewMesh = mOutputMesh->mListOfSideSets( iSet )->get_clusters_on_set();

            for ( uint iCluster = 0; iCluster < tSideClustersNewMesh.size(); iCluster++ )
            {
                // moris::Matrix< moris::IndexMat >              tSideOrdinalOld  = iCluster->get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellsOld = tSideClustersOldMesh( iCluster )->get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexOld       = tSideClustersOldMesh( iCluster )->get_vertices_in_cluster();
                uint                                       tIPCellIndexOld  = tSideClustersOldMesh( iCluster )->get_interpolation_cell_index();

                // moris::Matrix< moris::IndexMat >              tSideOrdinalNew = mSideClusters( iCounter ).get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellNew = tSideClustersNewMesh( iCluster )->get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexNew      = tSideClustersNewMesh( iCluster )->get_vertices_in_cluster();
                uint                                       tIPCellIndexNew = tSideClustersNewMesh( iCluster )->get_interpolation_cell_index();

                // compare the primary ig cell id and index
                bool tPrimaryCellIdAndIndexEqual = std::equal( tPrimaryCellsOld.cbegin(),
                        tPrimaryCellsOld.cend(),
                        tPrimaryCellNew.cbegin(),
                        []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                // compare the vertices id and index
                bool tVertexIdAndIndexEqual = std::equal( tVertexOld.begin(),
                        tVertexOld.end(),
                        tVertexNew.begin(),
                        []( Vertex const * a, Vertex const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                // compare the ip cell
                bool tEqualIPCell = tIPCellIndexOld == tIPCellIndexNew;

                // combine the outputs
                tOutput = tEqualIPCell
                       && tVertexIdAndIndexEqual
                       && tPrimaryCellIdAndIndexEqual;

                if ( !tOutput )
                {
                    MORIS_LOG_ERROR( "Side cluster number %u side sets of the mesh , is inconsistent,tEqualIPCell %d, vertices %d, primary cells %d  ",
                            iCluster,
                            tEqualIPCell,
                            tVertexIdAndIndexEqual,
                            tPrimaryCellIdAndIndexEqual );

                    return tOutput;
                }
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_double_sided_clusters()
    {
        // initialize the output
        bool tOutput = true;

        // get the old side sets
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // counter for the clusters
        uint iCounter = 0;

        // loop over the double sided sets
        for ( auto&& iSet : tDoubleSideSets.data() )
        {
            // if it is a ghost don't check yet
            if ( std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // get clusters on the set
            Vector< Cluster const * > tSideClusters = iSet->get_clusters_on_set();

            for ( const auto& iCluster : tSideClusters )
            {
                // obtain the leader side cluster
                mtk::Cluster const & tLeaderSideCluster = iCluster->get_leader_side_cluster();

                // get the member data of the old mesh
                moris::Matrix< moris::IndexMat >           tSideOrdinalOld  = tLeaderSideCluster.get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellsOld = tLeaderSideCluster.get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexOld       = tLeaderSideCluster.get_vertices_in_cluster();
                uint                                       tIPCellIndexOld  = tLeaderSideCluster.get_interpolation_cell_index();

                // get the member data of the new mesh
                moris::Matrix< moris::IndexMat >           tSideOrdinalNew = mOutputMesh->mDblSideClusters( iCounter ).get_leader_side_cluster().get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellNew = mOutputMesh->mDblSideClusters( iCounter ).get_leader_side_cluster().get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexNew      = mOutputMesh->mDblSideClusters( iCounter ).get_leader_side_cluster().get_vertices_in_cluster();
                uint                                       tIPCellIndexNew = mOutputMesh->mDblSideClusters( iCounter ).get_leader_side_cluster().get_interpolation_cell_index();

                // compare the local coordinates
                bool tEqualLocalCoords = std::equal( tSideOrdinalOld.begin(), tSideOrdinalOld.end(), tSideOrdinalNew.begin() );

                // compare the primary ig cell id and index
                bool tPrimaryCellIdAndIndexEqual = std::equal( tPrimaryCellsOld.cbegin(),
                        tPrimaryCellsOld.cend(),
                        tPrimaryCellNew.cbegin(),
                        []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                // compare the vertices id and index
                bool tVertexIdAndIndexEqual = std::equal( tVertexOld.begin(),
                        tVertexOld.end(),
                        tVertexNew.begin(),
                        []( Vertex const * a, Vertex const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                // compare the IP cell
                bool tEqualIPCell = tIPCellIndexOld == tIPCellIndexNew;

                // combine the outputs
                tOutput = tEqualIPCell
                       && tVertexIdAndIndexEqual
                       && tPrimaryCellIdAndIndexEqual
                       && tEqualLocalCoords;

                if ( !tOutput )
                {
                    MORIS_LOG_ERROR( "double sided cluster number %u of the mesh, is inconsistent, tEqualIPCell: %d , tVertexIdAndIndexEqual: %d ,tPrimaryCellIdAndIndexEqual: %d, tEqualLocalCoords: %d  ",
                            iCounter,
                            tEqualIPCell,
                            tVertexIdAndIndexEqual,
                            tPrimaryCellIdAndIndexEqual,
                            tEqualLocalCoords );

                    MORIS_LOG_ERROR( "double sided cluster %u of the mesh, is inconsistent, tEqualIPCell: %d , tVertexIdAndIndexEqual: %d ,tPrimaryCellIdAndIndexEqual: %d, tEqualLocalCoords: %d  ",
                            iCounter,
                            tEqualIPCell,
                            tVertexIdAndIndexEqual,
                            tPrimaryCellIdAndIndexEqual,
                            tEqualLocalCoords );

                    return tOutput;
                }

                iCounter++;
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_ghost_clusters()
    {
        // initialize the output
        bool tOutput = true;

        // get list of double sided sets
        Vector< moris::mtk::Double_Side_Set* > const & tDoubleSideSets = mInputMesh->get_double_side_sets();

        // counter for number of clusters
        uint iCounter = 0;

        // loop over the sets
        for ( const auto& iSet : tDoubleSideSets.data() )
        {
            // if it is not a ghost don't check
            if ( !std::strstr( iSet->get_set_name().c_str(), "ghost" ) )
            {
                continue;
            }

            // get clusters from the old mesh
            Vector< Cluster const * > tSideClusters = iSet->get_clusters_on_set();

            for ( const auto& iCluster : tSideClusters )
            {
                // get leader side
                mtk::Cluster const & tLeaderSideCluster = iCluster->get_leader_side_cluster();

                // memebr data for the old mesh
                moris::Matrix< moris::IndexMat >           tSideOrdinalOld  = tLeaderSideCluster.get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellsOld = tLeaderSideCluster.get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexOld       = tLeaderSideCluster.get_vertices_in_cluster();
                uint                                       tIPCellIndexOld  = tLeaderSideCluster.get_interpolation_cell_index();

                // member data for the new mesh
                moris::Matrix< moris::IndexMat >           tSideOrdinalNew = mOutputMesh->mGhostDblSidedSet( iCounter ).get_leader_side_cluster().get_cell_side_ordinals();
                Vector< moris::mtk::Cell const * > const & tPrimaryCellNew = mOutputMesh->mGhostDblSidedSet( iCounter ).get_leader_side_cluster().get_primary_cells_in_cluster();
                Vector< moris::mtk::Vertex const * >       tVertexNew      = mOutputMesh->mGhostDblSidedSet( iCounter ).get_leader_side_cluster().get_vertices_in_cluster();
                uint                                       tIPCellIndexNew = mOutputMesh->mGhostDblSidedSet( iCounter ).get_leader_side_cluster().get_interpolation_cell_index();

                // compare the local coordinates
                bool tEqualLocalCoords = std::equal( tSideOrdinalOld.begin(), tSideOrdinalOld.end(), tSideOrdinalNew.begin() );

                // compare the primary ig cell id and index
                bool tPrimaryCellIdAndIndexEqual = std::equal( tPrimaryCellsOld.cbegin(),
                        tPrimaryCellsOld.cend(),
                        tPrimaryCellNew.cbegin(),
                        []( mtk::Cell const * a, mtk::Cell const * b ) -> bool {
                            return a->get_index() == b->get_index();
                        } );

                // compare the vertices id and index
                bool tVertexIdAndIndexEqual = std::equal( tVertexOld.begin(),
                        tVertexOld.end(),
                        tVertexNew.begin(),
                        []( Vertex const * a, Vertex const * b ) -> bool {
                            return a->get_id() == b->get_id() && a->get_index() == b->get_index();
                        } );

                // compare the ip cell
                bool tEqualIPCell = tIPCellIndexOld == tIPCellIndexNew;

                // combine the output
                tOutput = tEqualIPCell
                       && tVertexIdAndIndexEqual
                       && tPrimaryCellIdAndIndexEqual
                       && tEqualLocalCoords;

                if ( !tOutput )
                {
                    MORIS_LOG_ERROR( "ghost cluster number %u of the mesh, is inconsistent,tEqualIPCell %d, tVertices %d, tPrimaryCells %d, tCoords %d ",
                            iCounter,
                            tEqualIPCell,
                            tVertexIdAndIndexEqual,
                            tPrimaryCellIdAndIndexEqual,
                            tEqualLocalCoords );

                    return tOutput;
                }

                // increment the count by 1
                iCounter++;
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Integration_Mesh_Editor::check_double_sided_sets()
    {
        // initialize the output
        bool tOutput = true;

        // get the old and new sets
        Vector< moris::mtk::Double_Side_Set* > const & tSideSetsOld = mInputMesh->get_double_side_sets();
        Vector< moris::mtk::Double_Side_Set* > const & tSideSetsNew = mOutputMesh->get_double_side_sets();

        // loop over the sets
        for ( uint iSet = 0; iSet < tSideSetsOld.size(); iSet++ )
        {
            // compare their names
            std::string tSideSetNameOld = tSideSetsOld( iSet )->get_set_name();
            std::string tSideSetNameNew = tSideSetsNew( iSet )->get_set_name();

            MORIS_ASSERT( tSideSetNameOld == tSideSetNameNew, "SideSet Names are not matching" );

            // get the old and new clusters
            Vector< Cluster const * > tSideClustersOldMesh = tSideSetsOld( iSet )->get_clusters_on_set();
            Vector< Cluster const * > tSideClustersNewMesh = tSideSetsNew( iSet )->get_clusters_on_set();

            // loop over the clusters
            for ( uint iCluster = 0; iCluster < tSideClustersNewMesh.size(); iCluster++ )
            {
                // get the leader and follower side clusters of the old and new mesh
                moris::mtk::Cluster const & tLeaderClusterOld   = tSideClustersOldMesh( iCluster )->get_leader_side_cluster();
                moris::mtk::Cluster const & tFollowerClusterOld = tSideClustersOldMesh( iCluster )->get_follower_side_cluster();

                moris::mtk::Cluster const & tLeaderClusterNew   = tSideClustersNewMesh( iCluster )->get_leader_side_cluster();
                moris::mtk::Cluster const & tFollowerClusterNew = tSideClustersNewMesh( iCluster )->get_follower_side_cluster();

                // check the leader side ordinals for the new and old mesh
                moris::Matrix< moris::IndexMat > tSideOrdLeaderOld = tLeaderClusterOld.get_cell_side_ordinals();
                moris::Matrix< moris::IndexMat > tSideOrdLeaderNew = tLeaderClusterNew.get_cell_side_ordinals();

                bool tEqSideOrdLeader = std::equal( tSideOrdLeaderOld.begin(), tSideOrdLeaderOld.end(), tSideOrdLeaderNew.begin() );

                // check the leader side primary cell ids for the new and old mesh
                moris::Matrix< moris::IndexMat > tCellIndexLeaderOld = tLeaderClusterOld.get_primary_cell_indices_in_cluster();
                moris::Matrix< moris::IndexMat > tCellIndexLeaderNew = tLeaderClusterNew.get_primary_cell_indices_in_cluster();

                bool tEqCellIndexLeader = std::equal( tCellIndexLeaderOld.begin(), tCellIndexLeaderOld.end(), tCellIndexLeaderNew.begin() );

                // check the leader side vertices ids for the new and old mesh
                moris::Matrix< moris::IndexMat > tVertexIndexLeaderOld = tLeaderClusterOld.get_vertex_indices_in_cluster();
                moris::Matrix< moris::IndexMat > tVertexIndexLeaderNew = tLeaderClusterNew.get_vertex_indices_in_cluster();

                bool tEqVertexIndexLeader = std::equal( tVertexIndexLeaderOld.begin(), tVertexIndexLeaderOld.end(), tVertexIndexLeaderNew.begin() );

                // check the leader side ip cell for the new and old mesh
                moris_index tIPCellLeaderOld = tLeaderClusterOld.get_interpolation_cell_index();
                moris_index tIPCellLeaderNew = tLeaderClusterNew.get_interpolation_cell_index();

                bool tEqIPCellLeader = tIPCellLeaderOld == tIPCellLeaderNew;

                // check the leader side trivial for the new and old mesh
                bool tTrivialLeaderOld = tLeaderClusterOld.is_trivial();
                bool tTrivialLeaderNew = tLeaderClusterNew.is_trivial();

                bool tEqTrivialLeader = tTrivialLeaderOld == tTrivialLeaderNew;

                // check the leader side local coordinates for the new and old mesh
                moris::Matrix< moris::DDRMat > tCoordsLeaderOld = tLeaderClusterOld.get_vertices_local_coordinates_wrt_interp_cell();
                moris::Matrix< moris::DDRMat > tCoordsLeaderNew = tLeaderClusterNew.get_vertices_local_coordinates_wrt_interp_cell();

                bool tEqCoordsLeader = std::equal( tCoordsLeaderOld.begin(), tCoordsLeaderOld.end(), tCoordsLeaderNew.begin() );

                // check the cluster measure for non-ghost side clusters
                bool tEqMeasureLeader = true;
                if ( !std::strstr( tSideSetNameOld.c_str(), "ghost" ) )
                {
                    moris::real tMeasureLeaderOld = tLeaderClusterOld.compute_cluster_cell_measure();
                    moris::real tMeasureLeaderNew = tLeaderClusterNew.compute_cluster_cell_measure();

                    tEqMeasureLeader = tMeasureLeaderNew == tMeasureLeaderOld;
                }

                // combine the leader side output
                bool tOutputLeader = tEqSideOrdLeader
                                  && tEqCellIndexLeader
                                  && tEqVertexIndexLeader
                                  && tEqIPCellLeader
                                  && tEqTrivialLeader
                                  && tEqCoordsLeader
                                  && tEqMeasureLeader;

                // Do the same for the Follower Side
                moris::Matrix< moris::IndexMat > tSideOrdFollowerOld = tFollowerClusterOld.get_cell_side_ordinals();
                moris::Matrix< moris::IndexMat > tSideOrdFollowerNew = tFollowerClusterNew.get_cell_side_ordinals();

                bool tEqSideOrdFollower = std::equal( tSideOrdFollowerOld.begin(), tSideOrdFollowerOld.end(), tSideOrdFollowerNew.begin() );

                moris::Matrix< moris::IndexMat > tCellIndexFollowerOld = tFollowerClusterOld.get_primary_cell_indices_in_cluster();
                moris::Matrix< moris::IndexMat > tCellIndexFollowerNew = tFollowerClusterNew.get_primary_cell_indices_in_cluster();

                bool tEqCellIndexFollower = std::equal( tCellIndexFollowerOld.begin(), tCellIndexFollowerOld.end(), tCellIndexFollowerNew.begin() );

                moris::Matrix< moris::IndexMat > tVertexIndexFollowerOld = tFollowerClusterOld.get_vertex_indices_in_cluster();
                moris::Matrix< moris::IndexMat > tVertexIndexFollowerNew = tFollowerClusterNew.get_vertex_indices_in_cluster();

                bool tEqVertexIndexFollower = std::equal( tVertexIndexFollowerOld.begin(), tVertexIndexFollowerOld.end(), tVertexIndexFollowerNew.begin() );

                moris_index tIPCellFollowerOld = tFollowerClusterOld.get_interpolation_cell_index();
                moris_index tIPCellFollowerNew = tFollowerClusterNew.get_interpolation_cell_index();

                bool tEqIPCellFollower = tIPCellFollowerOld == tIPCellFollowerNew;

                bool tTrivialFollowerOld = tFollowerClusterOld.is_trivial();
                bool tTrivialFollowerNew = tFollowerClusterNew.is_trivial();

                bool tEqTrivialFollower = tTrivialFollowerOld == tTrivialFollowerNew;

                moris::Matrix< moris::DDRMat > tCoordsFollowerOld = tFollowerClusterOld.get_vertices_local_coordinates_wrt_interp_cell();
                moris::Matrix< moris::DDRMat > tCoordsFollowerNew = tFollowerClusterNew.get_vertices_local_coordinates_wrt_interp_cell();

                bool tEqCoordsFollower = std::equal( tCoordsFollowerOld.begin(), tCoordsFollowerOld.end(), tCoordsFollowerNew.begin() );

                bool tEqMeasureFollower = true;
                if ( !std::strstr( tSideSetNameOld.c_str(), "ghost" ) )
                {
                    moris::real tMeasureFollowerOld = tFollowerClusterOld.compute_cluster_cell_measure();
                    moris::real tMeasureFollowerNew = tFollowerClusterNew.compute_cluster_cell_measure();

                    tEqMeasureFollower = tMeasureFollowerOld == tMeasureFollowerNew;
                }

                // combine the output for the follower side
                bool tOutputFollower = tEqSideOrdFollower
                                    && tEqCellIndexFollower
                                    && tEqVertexIndexFollower
                                    && tEqIPCellFollower
                                    && tEqTrivialFollower
                                    && tEqCoordsFollower
                                    && tEqMeasureFollower;

                tOutput = tOutputFollower && tOutputLeader;

                if ( !tOutput )
                {
                    // follower info
                    MORIS_LOG_ERROR( "Cluster number %u of the mesh, is inconsistent,tEqSideOrdFollower %d ,tEqCellIndexFollower%d ,tEqVertexIndexFollower %d ,tEqIPCellFollower %d ,tEqTrivialFollower %d , tEqCoordsFollower %d ,tEqMeasureFollower%d   ",
                            iCluster,
                            tEqSideOrdFollower,
                            tEqCellIndexFollower,
                            tEqVertexIndexFollower,
                            tEqIPCellFollower,
                            tEqTrivialFollower,
                            tEqCoordsFollower,
                            tEqMeasureFollower );

                    // leader info
                    MORIS_LOG_ERROR( "Cluster number %u of the mesh, is inconsistent,tEqSideOrdFollower %d ,tEqCellIndexFollower%d ,tEqVertexIndexFollower %d ,tEqIPCellFollower %d ,tEqTrivialFollower %d , tEqCoordsFollower %d ,tEqMeasureFollower%d   ",
                            iCluster,
                            tEqSideOrdLeader,
                            tEqCellIndexLeader,
                            tEqVertexIndexLeader,
                            tEqIPCellLeader,
                            tEqTrivialLeader,
                            tEqCoordsLeader,
                            tEqMeasureLeader );

                    return tOutput;
                }
            }
        }

        return tOutput;
    }

    //--------------------------------------------------------------------------------------------------------------
    bool
    Integration_Mesh_Editor::check_maps()
    {
        // obtain vertex map from the mesh
        std::unordered_map< moris_id, moris_index > tVertexGlobalIdToLocalIndex = mInputMesh->get_vertex_glb_id_to_loc_vertex_ind_map();

        // compare the vertex map
        bool tVertexMapEqual = mOutputMesh->mVertexGlobalIdToLocalIndex == tVertexGlobalIdToLocalIndex;

        // if they are not equal put a message
        if ( !tVertexMapEqual )
        {
            MORIS_LOG_ERROR( "Vertex Global To Local Map does not match!" );
        }

        return tVertexMapEqual;
    }

}    // namespace moris::mtk
