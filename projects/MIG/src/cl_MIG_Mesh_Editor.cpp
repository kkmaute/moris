/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MIG_Mesh_Editor.cpp
 *
 */

#include "cl_MTK_Set_Communicator.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Side_Cluster_DataBase.hpp"
#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_MTK_Mesh_DataBase_IP.hpp"
#include "cl_MTK_Integration_Mesh_Editor.hpp"
#include "cl_Matrix.hpp"
#include <numeric>

#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_MIG_Mesh_Editor.hpp"
#include "cl_MIG_Periodic_2D.hpp"
#include "cl_MIG_Periodic_3D.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_Tracer.hpp"

void moris::mig::Periodic_Mesh_Editor::merge_meshes()
{
    this->recreate_side_sets();

    // // double Sided Clusters

    // get number of cells and reserve enough space for cells
    uint tNumDblSideClusters = mIGMeshInfo->mDoubleSidedClusterToNewSideClusterIndex.size();

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
                [ this ]( moris_index aVertexIndex ) -> mtk::Vertex const * { return &mOutputMesh->mVertices( aVertexIndex ); } );

        // constrcut the double sided cluster
        mtk::Double_Side_Cluster tDblSideCluster( tLeaderCluster, tFollowerCluster, tVertexPair );

        // save the constructed double sided cluster
        mOutputMesh->mDblSideClusters( iDblCluster ) = tDblSideCluster;
    }

    // // Double Sided Set Data

    // get the old double sided sets
    Vector< mtk::Double_Side_Set* > tDoubleSideSets = mOutputMesh->get_double_side_sets();

    // counter for the DoubleSidedSets
    size_t iCounter = 0;

    // loop over the old double sided sets
    for ( const auto& iSet : tDoubleSideSets.data() )
    {
        if ( iCounter == ( mNumPreviousDoubleSideSet ) ) break;

        if ( strstr( iSet->get_set_name().c_str(), "P" ) )
        {
            break;
        }

        // if the double sided cluster is ghost skip it
        if ( strstr( iSet->get_set_name().c_str(), "ghost" ) )
        {
            break;
        }

        // create the list of clusters inside the set
        Vector< mtk::Cluster const * > aSideSetClusters;

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
                [ this ]( moris_index mClusterIndex ) -> mtk::Cluster const * { return &mOutputMesh->mDblSideClusters( mClusterIndex ); } );

        // construct the double sided set
        mOutputMesh->mListOfDoubleSideSets( iCounter ) =
                new mtk::Double_Side_Set( iSet->get_set_name(),
                        aSideSetClusters,
                        iSet->get_set_colors(),
                        mIGMeshInfo->mSpatialDim );

        // increment the set count by 1
        iCounter++;

    }    // end for: each double sided side set

    // communicate information across all sets
    mtk::Set_Communicator tSetCommunicator( mOutputMesh->mListOfDoubleSideSets );

    // delete the old data to prevent memory leaks
    for ( int iSet = iCounter - 1; iSet >= 0; iSet-- )
    {
        delete tDoubleSideSets( iSet );
    }

    // collet the sets and make collect them
    mOutputMesh->collect_all_sets();
}
//------------------------------------------------------------------------------------------------------------
void moris::mig::Periodic_Mesh_Editor::reconstruct_connectivity()
{
    uint iCounter = 0;
    // populate the cell to vertex connectivity list
    for ( const int& iVertex : mIGMeshInfo->mCellToVertexIndicies )
    {
        mOutputMesh->mCellToVertices( iCounter ) = &mOutputMesh->mVertices( iVertex );
        iCounter++;
    }

    iCounter = 0;
    // populate the cell cluster to primary cell connectivity list
    for ( const int& iPrimaryCell : mIGMeshInfo->mCellClusterToPrimaryIGCellIndices )
    {
        mOutputMesh->mCellClusterToPrimaryIGCell( iCounter ) = &mOutputMesh->mCells( iPrimaryCell );
        iCounter++;
    }

    iCounter = 0;
    // populate the cell cluster to void cell connectivity list
    for ( const int& iVoidCell : mIGMeshInfo->mCellClusterToVoidIGCellIndices )
    {
        mOutputMesh->mCellClusterToVoidIGCell( iCounter ) = &mOutputMesh->mCells( iVoidCell );
        iCounter++;
    }

    iCounter = 0;
    // populate the cell cluster to vertices connectivity list
    for ( const int& iVertex : mIGMeshInfo->mCellClusterToVertexIndices )
    {
        mOutputMesh->mCellClusterToVeretx( iCounter ) = &mOutputMesh->mVertices( iVertex );
        iCounter++;
    }

    iCounter = 0;
    // populate the cell cluster to vertices connectivity list
    for ( const int& iPrimaryCell : mIGMeshInfo->mSideClusterToPrimaryIGCellIndices )
    {
        mOutputMesh->mSideClusterToPrimaryIGCell( iCounter ) = &mOutputMesh->mCells( iPrimaryCell );
        iCounter++;
    }

    iCounter = 0;
    // populate the cell cluster to vertices connectivity list
    for ( const int& iVertex : mIGMeshInfo->mSideClusterToVertexIndices )
    {
        mOutputMesh->mSideClusterToVeretx( iCounter ) = &mOutputMesh->mVertices( iVertex );
        iCounter++;
    }


    iCounter = 0;
    // populate the cluster to vertices,  convert indices to pointers
    for ( uint iVertex = 0; iVertex < mIGMeshInfo->mGhostLeaderFollowerVertexIndices.size(); iVertex++ )
    {
        mOutputMesh->mGhostLeaderFollowerToVertex( iCounter ) = &mOutputMesh->mVertices( mIGMeshInfo->mGhostLeaderFollowerVertexIndices( iVertex ) );
        iCounter++;
    }

    // get number of the double sided clusters
    uint tNumGhostClusters = mIGMeshInfo->mGhostLeaderToIPCellIndex.size();

    iCounter = 0;
    // loop over the number of the clusters to add leader and follower in order
    for ( uint iCell = 0; iCell < tNumGhostClusters; iCell++ )
    {
        // add the ig cell pointer of the leader and then salve
        mOutputMesh->mGhostLeaderFollowerIGCellList( 2 * iCounter )     = &mOutputMesh->get_mtk_cell( mIGMeshInfo->mGhostLeaderToIGCellIndex( iCell ) );
        mOutputMesh->mGhostLeaderFollowerIGCellList( 2 * iCounter + 1 ) = &mOutputMesh->get_mtk_cell( mIGMeshInfo->mGhostFollowerToIGCellIndex( iCell ) );

        iCounter++;
    }

    // populate the outward data that is being returned
    for ( auto& iCluster : mOutputMesh->mCellClusters )
    {
        iCluster.set_outward_data();
    }

    // populate the outward data that is being returned
    for ( auto& iCluster : mOutputMesh->mSideClusters )
    {
        iCluster.set_outward_data();
    }


    // get number of the double sided clusters
    uint tNumGhostClusters2 = mIGMeshInfo->mGhostLeaderToIPCellIndex.size();

    // get the starting index of the ghost clusters
    moris_index tGhostOverallClusterIndex = mOutputMesh->mSideClusters.size();

    // loop over number of double sided clusters to create the new ghost clusters
    for ( uint iCluster = 0; iCluster < tNumGhostClusters2; iCluster++ )
    {
        mOutputMesh->mGhostLeader( iCluster ).update_cluster_index( tGhostOverallClusterIndex++ );
        mOutputMesh->mGhostLeader( iCluster ).set_outward_data();

        mOutputMesh->mGhostFollower( iCluster ).update_cluster_index( tGhostOverallClusterIndex++ );
        mOutputMesh->mGhostFollower( iCluster ).set_outward_data();

        // initialize the vertexpair list for the double sided cluster
        Vector< mtk::Vertex const * > tVertexPair( mIGMeshInfo->mGhostToVertexOffset( iCluster + 1 ) - mIGMeshInfo->mGhostToVertexOffset( iCluster ) );

        // transform indices to pointers
        std::transform( mIGMeshInfo->mGhostToVertexPairIndices.begin() + mIGMeshInfo->mGhostToVertexOffset( iCluster ),
                mIGMeshInfo->mGhostToVertexPairIndices.begin() + mIGMeshInfo->mGhostToVertexOffset( iCluster + 1 ),
                tVertexPair.begin(),
                [ this ]( moris_index aVertexIndex ) -> mtk::Vertex const * { return &mOutputMesh->mVertices( aVertexIndex ); } );

        // constrcut the double sided clusters
        mtk::Double_Side_Cluster tDblSideCluster( &mOutputMesh->mGhostLeader( iCluster ), &mOutputMesh->mGhostFollower( iCluster ), tVertexPair );

        // save the double sided clusters
        mOutputMesh->mGhostDblSidedSet( iCluster ) = tDblSideCluster;
    }
}
namespace moris::mig
{
    //------------------------------------------------------------------------------------------------------------
    Periodic_Mesh_Editor::Periodic_Mesh_Editor( mtk::Integration_Mesh_DataBase_IG* aIGMesh, mig::Periodic_2D* aPeriodicData2D )
            : mPeriodicData2D( aPeriodicData2D )
    {
        mOutputMesh = aIGMesh;
    }

    Periodic_Mesh_Editor::Periodic_Mesh_Editor( mtk::Integration_Mesh_DataBase_IG* aIGMesh, mig::Periodic_3D* aPeriodicData3D )
            : mPeriodicData3D( aPeriodicData3D )
    {
        mOutputMesh = aIGMesh;
    }
    //------------------------------------------------------------------------------------------------------------

    Periodic_Mesh_Editor::~Periodic_Mesh_Editor()
    {
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_Mesh_Editor::perform()
    {
        Tracer tTracer( "MIG", "Editor", "Merging Meshes" );

        // link the newly created nodes to the geomtry engine
        this->link_nodes_to_geomtry_engine();

        // depending on the dimension add data to the data base
        if ( mPeriodicData3D == nullptr )
        {
            this->construct_periodic_data_base(
                    mPeriodicData2D->mSideClusterToVertexIndices,
                    mPeriodicData2D->mVerticesCoords,
                    mPeriodicData2D->mSideClusterToCells,
                    mPeriodicData2D->mCellToVertexIndices,
                    mPeriodicData2D->mSideClusterToIPCell,
                    mPeriodicData2D->mVertexParametricCoords,
                    mPeriodicData2D->mDoubleSidedClustersIndex,
                    mPeriodicData2D->mNumDblSideCluster,
                    mGeometryEngine->get_num_phases() );
        }

        else
        {
            this->construct_periodic_data_base(
                    mPeriodicData3D->mSideClusterToVertexIndices,
                    mPeriodicData3D->mVerticesCoords,
                    mPeriodicData3D->mSideClusterToCells,
                    mPeriodicData3D->mCellToVertexIndices,
                    mPeriodicData3D->mSideClusterToIPCell,
                    mPeriodicData3D->mVertexParametricCoords,
                    mPeriodicData3D->mDoubleSidedClustersIndex,
                    mPeriodicData3D->mNumDblSideCluster,
                    mGeometryEngine->get_num_phases() );
        }
    }

    //------------------------------------------------------------------------------------------------------------
    void
    Periodic_Mesh_Editor::set_geometry_engine( moris::gen::Geometry_Engine* aGeometryEngine )
    {
        mGeometryEngine = aGeometryEngine;
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_Mesh_Editor::link_nodes_to_geomtry_engine()
    {
        // get number of current nodes
        uint tNumPreviousVertices = mOutputMesh->get_num_nodes();

        // update underlying ids and owners of interpolation nodes in GE
        if ( mOutputMesh->get_spatial_dim() == 2 )
        {
            for ( uint iVertex = 0; iVertex < mPeriodicData2D->mNumVertices; iVertex++ )
            {
                moris_index tNodeIndex = iVertex + tNumPreviousVertices;

                // FIXME - it is for serial
                moris_id    tNodeId    = tNodeIndex + 1;
                moris_index tNodeOwner = 0;

                // add nodes to the ge
                mGeometryEngine->update_intersection_node( tNodeIndex, tNodeId, tNodeOwner );
            }
        }
        else
        {
            for ( uint iVertex = 0; iVertex < mPeriodicData3D->mNumVertices; iVertex++ )
            {
                moris_index tNodeIndex = iVertex + tNumPreviousVertices;

                // FIXME - it is for serial
                moris_id    tNodeId    = tNodeIndex + 1;
                moris_index tNodeOwner = 0;

                // add nodes to the ge
                mGeometryEngine->update_intersection_node( tNodeIndex, tNodeId, tNodeOwner );
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------

    void
    Periodic_Mesh_Editor::construct_periodic_data_base(
            Vector< Vector< moris_index > >& aSideClusterToVertexIndices,
            Matrix< DDRMat >                 aVerticesCoords,
            Vector< Vector< moris_index > >& aSideClusterToCells,
            Vector< Vector< moris_index > >& aCellToVertexIndices,
            Vector< moris_index >&           aSideClusterToIPCell,
            Matrix< DDRMat >&                aVertexParametricCoords,
            Vector< moris_index >&           aDoubleSidedClustersIndex,
            uint                             mNumDblSideCluster,
            uint                             aNumGeometry )
    {
        // get the previous sizes to
        mNumPreviousVertices       = mOutputMesh->mVertices.size();
        mNumPreviousCells          = mOutputMesh->mCells.size();
        mNumPreviousSideCluster    = mOutputMesh->mSideClusters.size();
        mNumPreviousDblSideCluster = mOutputMesh->mDblSideClusters.size();
        mNumPreviousDoubleSideSet  = mOutputMesh->mListOfDoubleSideSets.size();

        // overwrite the mesh info
        mIGMeshInfo = mOutputMesh->mIGMeshInfo;

        // add the data to the database
        this->add_vertices( aSideClusterToVertexIndices, aVerticesCoords );
        this->add_cells( aSideClusterToCells, aCellToVertexIndices );
        this->add_side_clusters( aSideClusterToCells,
                aSideClusterToIPCell,
                aVertexParametricCoords,
                aSideClusterToVertexIndices );
        this->add_double_sided_clusters( mNumDblSideCluster, aSideClusterToVertexIndices );
        this->add_double_sided_set( aDoubleSidedClustersIndex, aNumGeometry );

        // reconstruct connectivity that is based in pointers
        this->reconstruct_connectivity();

        // reconstrcut double sided sets
        this->merge_meshes();
    }
}    // namespace moris::mig
