/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Ghost_Stabilization.cpp
 *
 */
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Enums.hpp"
#include "fn_norm.hpp"

#include "fn_stringify_matrix.hpp"
#include "fn_determine_cell_topology.hpp"

namespace xtk
{
    Ghost_Stabilization::Ghost_Stabilization()
            : mXTKModel( nullptr )
    {
    }

    // ----------------------------------------------------------------------------------

    Ghost_Stabilization::Ghost_Stabilization( Model* aXTKModel )
            : mXTKModel( aXTKModel )
    {
        mMinMeshIndex = mXTKModel->get_enriched_interp_mesh( 0 ).mMeshIndices.min();
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::setup_ghost_stabilization()
    {
        // get access to the B-spline mesh information
        mMeshIndices = mXTKModel->mEnrichment->get_mesh_indices();

        // initialize object carrying generated data for ghost facets
        Ghost_Setup_Data tGhostSetupData;

        // construct trivial subphase interpolation cells
        this->construct_ip_ig_cells_for_ghost_side_clusters( tGhostSetupData );

        // setup the side sets
        this->declare_ghost_double_side_sets_in_mesh( tGhostSetupData );

        // Construct Ghost Double Side Clusters and Sets
        this->construct_ghost_double_side_sets_in_mesh( tGhostSetupData );

        // Look through the vertices used in ghost stabilization
        // and identify which ones do not have their t-matrix.
        // Retrieve these t-matrices from their owner. This ensures
        // the solver has all appropriate information downstream.
        this->identify_and_setup_aura_vertices_in_ghost( tGhostSetupData );
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::setup_ghost_stabilization_new()
    {
        // get access to the B-spline mesh information
        mMeshIndices      = mXTKModel->mEnrichment->get_mesh_indices();
        mBsplineMeshInfos = mXTKModel->mCutIntegrationMesh->get_bspline_mesh_info();

        // initialize object carrying generated data for ghost facets
        Ghost_Setup_Data tGhostSetupData;

        // check whether a linear Lagrange mesh is used
        tGhostSetupData.mLinearBackgroundMesh = this->is_linear_ip_mesh();

        // setup the side sets
        this->declare_ghost_double_side_sets_in_mesh_new( tGhostSetupData );

        // Construct Ghost Double Side Clusters and Sets
        this->construct_ghost_double_side_sets_in_mesh_new( tGhostSetupData );

        // Look through the vertices used in ghost stabilization
        // and identify which ones do not have their t-matrix.
        // Retrieve these t-matrices from their owner. This ensures
        // the solver has all appropriate information downstream.
        this->identify_and_setup_aura_vertices_in_ghost( tGhostSetupData );
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::visualize_ghost_on_mesh( moris_index const & aBulkPhase )
    {
        Tracer tTracer( "XTK", "Ghost", "visualize_ghost_on_mesh" );

        // get the enriched integration mesh
        Enriched_Integration_Mesh& tEnrIgMesh = mXTKModel->get_enriched_integ_mesh( 0 );

        // get the set name
        std::string tGhostName = this->get_ghost_dbl_side_set_name( aBulkPhase );

        // get the double side set index in the mesh
        moris_index tDSSIndexInMesh = tEnrIgMesh.get_double_sided_set_index( tGhostName );

        // get the cell topology for the IP cells for visualization
        mtk::CellTopology tFacetTopo = determine_cell_topology( mXTKModel->get_spatial_dim(), mtk::Interpolation_Order::LINEAR, mtk::CellShape::RECTANGULAR );

        moris_index tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set( tDSSIndexInMesh, "ghost_ss_" + std::to_string( aBulkPhase ) );
        tEnrIgMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_" + std::to_string( aBulkPhase ), tFacetTopo );

        tEnrIgMesh.setup_color_to_set();
        tEnrIgMesh.collect_all_sets( false );
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::visualize_ghost_on_mesh_new(
            moris_index const & aBsplineMeshListIndex,
            moris_index const & aBulkPhaseIndex )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Ghost", "visualize ghost on mesh (new approach)" );

        // get the enriched integration mesh
        Enriched_Integration_Mesh& tEnrIgMesh = mXTKModel->get_enriched_integ_mesh( 0 );

        // get the B-spline mesh index
        moris_index tBspMeshIndex = mMeshIndices( aBsplineMeshListIndex );

        // get the set name
        std::string tGhostName = this->get_ghost_dbl_side_set_name( aBulkPhaseIndex, tBspMeshIndex );

        // get the double side set index in the mesh
        moris_index tDSSIndexInMesh = tEnrIgMesh.get_double_sided_set_index( tGhostName );

        // create strings for Ghost Side set and block set names
        std::string tGhostSsName = "Vis_Ghost_SS_B" + std::to_string( tBspMeshIndex ) + "_p" + std::to_string( aBulkPhaseIndex );
        std::string tGhostBsName = "Vis_Ghost_BS_B" + std::to_string( tBspMeshIndex ) + "_p" + std::to_string( aBulkPhaseIndex );

        // get the index of the side set
        moris_index tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set( tDSSIndexInMesh, tGhostSsName );

        // get the cell topology for the IP cells for visualization
        mtk::CellTopology tFacetTopo = determine_cell_topology( mXTKModel->get_spatial_dim(), mtk::Interpolation_Order::LINEAR, mtk::CellShape::RECTANGULAR );

        // create a block set for
        tEnrIgMesh.create_block_set_from_cells_of_side_set( tSSIndex, tGhostBsName, tFacetTopo );

        // finalize sets
        tEnrIgMesh.setup_color_to_set();
        tEnrIgMesh.collect_all_sets( false );
    }

    // ----------------------------------------------------------------------------------

    std::string
    Ghost_Stabilization::get_ghost_dbl_side_set_name( moris_index const & aBulkPhase )
    {
        MORIS_ASSERT( aBulkPhase < (moris_index)mXTKModel->get_geom_engine()->get_num_bulk_phase(),
                "Ghost_Stabilization::get_ghost_dbl_side_set_name() - Bulk Phase index out of bounds." );

        return "ghost_p" + std::to_string( aBulkPhase );
    }

    // ----------------------------------------------------------------------------------

    std::string
    Ghost_Stabilization::get_ghost_dbl_side_set_name(
            moris_index const & aBulkPhase,
            moris_index const & aBspMeshIndex )
    {
        MORIS_ASSERT( aBulkPhase < (moris_index)mXTKModel->get_geom_engine()->get_num_bulk_phase(),
                "Ghost_Stabilization::get_ghost_dbl_side_set_name_new() - Bulk Phase index out of bounds." );

        return "ghost_B" + std::to_string( aBspMeshIndex ) + "_p" + std::to_string( aBulkPhase );
    }

    // ----------------------------------------------------------------------------------

    Memory_Map
    Ghost_Stabilization::get_memory_usage()
    {
        // Ghost is an algorithm, all data created by
        // ghost is placed in the meshes
        Memory_Map tMM;
        tMM.mMemoryMapData[ "mXTKModel ptr" ] = sizeof( mXTKModel );
        return tMM;
    }

    // ----------------------------------------------------------------------------------

    // FIXME: Keenan - code below needs to be split up in smaller units
    void
    Ghost_Stabilization::construct_ip_ig_cells_for_ghost_side_clusters( Ghost_Setup_Data& aGhostSetupData )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Ghost", "construct_ip_ig_cells_for_ghost_side_clusters" );

        // access enriched integration mesh and enriched interp mesh
        Enriched_Integration_Mesh&   tEnrIgMesh = mXTKModel->get_enriched_integ_mesh( 0 );
        Enriched_Interpolation_Mesh& tEnrIpMesh = mXTKModel->get_enriched_interp_mesh( 0 );

        // check whether a linear Lagrange mesh is used
        aGhostSetupData.mLinearBackgroundMesh = this->is_linear_ip_mesh();

        // get the groupings of interpolation cells
        Cell< Interpolation_Cell_Unzipped* >         tOwnedInterpCells;
        Cell< Cell< Interpolation_Cell_Unzipped* > > tNotOwnedInterpCells;
        Cell< uint >                                 tProcRanks;
        tEnrIpMesh.get_owned_and_not_owned_enriched_interpolation_cells( tOwnedInterpCells, tNotOwnedInterpCells, tProcRanks );

        // get new index for new ghost interpolation cell, start where the enriched Ip mesh stops
        uint tCurrentNewInterpCellIndex = tEnrIpMesh.get_num_entities( mtk::EntityRank::ELEMENT );

        // allocate data in ghost setup data
        Cut_Integration_Mesh*                                 tCutIgMesh            = mXTKModel->get_cut_integration_mesh();
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseNeighborhood = tCutIgMesh->get_subphase_neighborhood();
        aGhostSetupData.mSubphaseIndexToInterpolationCellIndex.resize( tSubphaseNeighborhood->mSubphaseToSubPhase.size(), MORIS_INDEX_MAX );

        // figure out how many of the interpolation cells are not trivial (this value is the number of new interpolation cells)
        // iterate through owned interpolation cells and give them an id if they are not trivial clusters
        // Also, we're only interested in the nontrivial cell clusters so keep track of those.
        uint tNumNewInterpCellsNotOwned = 0;
        uint tNumNewInterpCellsOwned    = 0;

        Cell< Interpolation_Cell_Unzipped* > tNonTrivialOwnedInterpCells;

        // go over all enriched Ip Cells on processor
        for ( moris::size_t iEnrIpCell = 0; iEnrIpCell < tOwnedInterpCells.size(); iEnrIpCell++ )
        {
            // get the cluster and the subphase index associated with it
            xtk::Cell_Cluster const & tCluster  = tEnrIgMesh.get_xtk_cell_cluster( *tOwnedInterpCells( iEnrIpCell ) );
            moris_index               tSubphase = tOwnedInterpCells( iEnrIpCell )->get_subphase_index();

            // if the cluster is non-trivial ...
            if ( !tCluster.is_trivial() )
            {
                // ... make a 'note' that the enriched IP cell underlying the cluster needs to be used for creating a Ghost Facet
                tNumNewInterpCellsOwned++;
                tNonTrivialOwnedInterpCells.push_back( tOwnedInterpCells( iEnrIpCell ) );
                tCurrentNewInterpCellIndex++;
            }

            // store the enriched IP cell index associated with the current subphase in Ghost setup data
            aGhostSetupData.mSubphaseIndexToInterpolationCellIndex( tSubphase ) = tOwnedInterpCells( iEnrIpCell )->get_index();
        }

        // allocate new interpolation cell ids
        moris_id tCurrentId = tEnrIpMesh.allocate_entity_ids( tNumNewInterpCellsOwned, mtk::EntityRank::ELEMENT, false );

        // initialize maps that give each new non-trivial owned enriched Interp cell its global IDs
        Cell< moris_id >                         tNewNonTrivialOwnedInterpCellsIds( tNumNewInterpCellsOwned );
        std::unordered_map< moris_id, moris_id > tBaseEnrIdToIndexInNonTrivialOwned;

        // assign new non trivial owned enriched Interp cell Ids
        for ( moris::size_t iNtIpCell = 0; iNtIpCell < tNonTrivialOwnedInterpCells.size(); iNtIpCell++ )
        {
            // list IDs for new non trivial owned enriched Interp cell
            tNewNonTrivialOwnedInterpCellsIds( iNtIpCell ) = tCurrentId++;

            // map giving the list index for a new cell's ID
            tBaseEnrIdToIndexInNonTrivialOwned[ tNonTrivialOwnedInterpCells( iNtIpCell )->get_id() ] = iNtIpCell;
        }

        Cell< Cell< Interpolation_Cell_Unzipped* > > tNonTrivialNotOwnedInterpCells( tNotOwnedInterpCells.size() );

        // setup requests for not owned not trivial
        for ( moris::size_t iProc = 0; iProc < tNotOwnedInterpCells.size(); iProc++ )
        {
            for ( moris::size_t iIpCell = 0; iIpCell < tNotOwnedInterpCells( iProc ).size(); iIpCell++ )
            {
                // get the cluster and the subphase index associated with it
                xtk::Cell_Cluster const & tCluster  = tEnrIgMesh.get_xtk_cell_cluster( *tNotOwnedInterpCells( iProc )( iIpCell ) );
                moris_index               tSubphase = tNotOwnedInterpCells( iProc )( iIpCell )->get_subphase_index();

                // if the cluster is non-trivial ...
                if ( !tCluster.is_trivial() )
                {
                    // ... make a 'note' that the enriched IP cell underlying the cluster needs to be used for creating a Ghost Facet
                    tNumNewInterpCellsNotOwned++;
                    tNonTrivialNotOwnedInterpCells( iProc ).push_back( tNotOwnedInterpCells( iProc )( iIpCell ) );
                    tCurrentNewInterpCellIndex++;
                }

                // store the enriched IP cell index associated with the current subphase in Ghost setup data
                aGhostSetupData.mSubphaseIndexToInterpolationCellIndex( tSubphase ) = tNotOwnedInterpCells( iProc )( iIpCell )->get_index();
            }
        }
    }
    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_ip_cell_id_answers(
            Cell< Matrix< IndexMat > >&               aReceivedEnrCellIds,
            Cell< moris_id >&                         aNewInterpCellIds,
            Cell< Matrix< IndexMat > >&               aEnrCellIds,
            std::unordered_map< moris_id, moris_id >& aBaseEnrIdToIndexInNonTrivialOwned )
    {
        aEnrCellIds.resize( aReceivedEnrCellIds.size() );

        for ( uint iP = 0; iP < aReceivedEnrCellIds.size(); iP++ )
        {
            aEnrCellIds( iP ).resize( 1, aReceivedEnrCellIds( iP ).numel() );

            if ( aReceivedEnrCellIds( iP )( 0, 0 ) != MORIS_INDEX_MAX )
            {

                for ( uint iC = 0; iC < aReceivedEnrCellIds( iP ).numel(); iC++ )
                {
                    auto tIter = aBaseEnrIdToIndexInNonTrivialOwned.find( aReceivedEnrCellIds( iP )( iC ) );

                    if ( tIter == aBaseEnrIdToIndexInNonTrivialOwned.end() )
                    {
                        MORIS_ASSERT( tIter != aBaseEnrIdToIndexInNonTrivialOwned.end(), "Enriched cell id not in map" );
                    }
                    aEnrCellIds( iP )( iC ) = aNewInterpCellIds( tIter->second );
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::identify_and_setup_aura_vertices_in_ghost( Ghost_Setup_Data& aGhostSetupData )
    {
        // log this function when verbose output is requested
        Tracer tTracer( "XTK", "Ghost", "identify and setup aura vertices in ghost" );

        // access the enriched ip mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Get the communication table and map */

        // get the communication table
        Matrix< IdMat > tCommTable = mXTKModel->get_communication_table();

        /* ---------------------------------------------------------------------------------------- */
        /* Communicate the T-matrices wrt. every B-spline mesh*/

        for ( uint iMeshIndex = 0; iMeshIndex < tEnrInterpMesh.mMeshIndices.numel(); iMeshIndex++ )
        {
            // get the current B-spline mesh index
            moris_index tMeshIndex = tEnrInterpMesh.mMeshIndices( iMeshIndex );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 0: Figure out which T-matrices need to be communicated with whom */

            // unzipped vertices in ghost that do not have interpolation
            Cell< mtk::Vertex* > tGhostVerticesWithoutInterpolation;

            // an interpolation cell that the ghost vertex without interpolation is connected to
            // this is needed for communication routine
            Cell< mtk::Cell const * > tGhostIpCellConnectedToVertex;

            // find those UIPVs without a T-matrix wrt the current B-spline mesh
            this->collect_unzipped_IP_vertices_without_interpolation(
                    aGhostSetupData,
                    tGhostVerticesWithoutInterpolation,
                    tGhostIpCellConnectedToVertex,
                    tMeshIndex );

            /* ---------------------------------------------------------------------------------------- */
            /* The following steps are only necessary if code runs in parallel */

            if ( par_size() == 1 )    // serial
            {
                // check that all entities are owned in serial
                MORIS_ASSERT( tGhostVerticesWithoutInterpolation.size() == 0,
                        "Ghost_Stabilization::identify_and_setup_aura_vertices_in_ghost() - "
                        "Code running in serial, but there are nodes with missing T-matrices." );
            }
            else    // parallel
            {
                /* ---------------------------------------------------------------------------------------- */
                /* Prepare requests for non-owned entities */

                // initialize lists of information that identifies entities (on other procs)
                Cell< Cell< moris_index > > tIPVertIndsToProcs;    // UIPV indices for communication (local to current proc, just used for construction of arrays)
                Cell< Matrix< IdMat > >     tBaseVertexIds;        // base IP vertex IDs the UIPVs are constructed from
                Cell< Matrix< IdMat > >     tUnzippedIpCellIds;    // UIPC IDs these vertices belong to

                // initialize a map relating position in communication array to position in list of all requested
                Cell< Cell< moris_index > > tNotOwnedIPVertIndsInNotOwnedList;

                // fill identifying information
                this->prepare_requests_for_T_matrices_without_interpolation(
                        tGhostVerticesWithoutInterpolation,
                        tGhostIpCellConnectedToVertex,
                        tNotOwnedIPVertIndsInNotOwnedList,
                        tIPVertIndsToProcs,
                        tBaseVertexIds,
                        tUnzippedIpCellIds );

                /* ---------------------------------------------------------------------------------------- */
                /* Step 4: Send and Receive requests about non-owned entities to and from other procs */

                // initialize arrays for receiving
                Cell< Matrix< IdMat > > tReceivedBaseVertexIds;
                Cell< Matrix< IdMat > > tReceivedUnzippedIpCellIds;

                // communicate information
                moris::communicate_mats( tCommTable, tBaseVertexIds, tReceivedBaseVertexIds );
                moris::communicate_mats( tCommTable, tUnzippedIpCellIds, tReceivedUnzippedIpCellIds );

                // clear memory not needed anymore
                tBaseVertexIds.clear();
                tUnzippedIpCellIds.clear();

                /* ---------------------------------------------------------------------------------------- */
                /* Step 5: Find answers to the requests */

                // initialize lists of ID answers to other procs
                Cell< Matrix< DDRMat > >   tTMatrixWeights;
                Cell< Matrix< IdMat > >    tTMatrixIds;
                Cell< Matrix< IdMat > >    tTMatrixOwners;
                Cell< Matrix< IndexMat > > tTMatrixOffsets;

                this->prepare_answers_for_T_matrices(
                        tReceivedBaseVertexIds,
                        tReceivedUnzippedIpCellIds,
                        tTMatrixWeights,
                        tTMatrixIds,
                        tTMatrixOwners,
                        tTMatrixOffsets,
                        tMeshIndex );

                // clear memory from requests (the answers to which have been found)
                tReceivedBaseVertexIds.clear();
                tReceivedUnzippedIpCellIds.clear();

                /* ---------------------------------------------------------------------------------------- */
                /* Step 6: Send and receive answers to and from other procs */

                // initialize arrays for receiving
                Cell< Matrix< DDRMat > >   tReceivedTMatrixWeights;
                Cell< Matrix< IdMat > >    tReceivedTMatrixIds;
                Cell< Matrix< IdMat > >    tReceivedTMatrixOwners;
                Cell< Matrix< IndexMat > > tReceivedTMatrixOffsets;

                // communicate answers
                moris::communicate_mats( tCommTable, tTMatrixWeights, tReceivedTMatrixWeights );
                moris::communicate_mats( tCommTable, tTMatrixIds, tReceivedTMatrixIds );
                moris::communicate_mats( tCommTable, tTMatrixOwners, tReceivedTMatrixOwners );
                moris::communicate_mats( tCommTable, tTMatrixOffsets, tReceivedTMatrixOffsets );

                // clear unused memory
                tTMatrixWeights.clear();
                tTMatrixIds.clear();
                tTMatrixOwners.clear();
                tTMatrixOffsets.clear();

                /* ---------------------------------------------------------------------------------------- */
                /* Step 7: Store answers for not owned entities */

                this->handle_requested_T_matrix_answers(
                        tGhostIpCellConnectedToVertex,
                        tNotOwnedIPVertIndsInNotOwnedList,
                        tIPVertIndsToProcs,
                        tReceivedTMatrixWeights,
                        tReceivedTMatrixIds,
                        tReceivedTMatrixOwners,
                        tReceivedTMatrixOffsets,
                        tMeshIndex );

            }    // end if: parallel

        }        // end for: every B-spline mesh index

    }            // end function: Ghost_Stabilization::identify_and_setup_aura_vertices_in_ghost

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::collect_unzipped_IP_vertices_without_interpolation(
            Ghost_Setup_Data&          aGhostSetupData,
            Cell< mtk::Vertex* >&      aGhostVerticesWithoutInterpolation,
            Cell< mtk::Cell const * >& aGhostIpCellConnectedToVertex,
            const moris_index          aMeshIndex )
    {
        // get the number of ghost double side sets
        uint tNumGhostSets = aGhostSetupData.mDblSideSetIndexInMesh.size();

        // map listing vertices that are already marked for communication
        std::unordered_map< moris_index, bool > tGhostVerticesWithOutInterpolationMap;

        // iterate through sets and then clusters and gather the T-matrix information from those needed
        for ( uint iGhostSet = 0; iGhostSet < tNumGhostSets; iGhostSet++ )
        {
            // get the ghost clusters on the current set
            Cell< mtk::Cluster const * > tDblSideSetClusters =
                    mXTKModel->get_enriched_integ_mesh( 0 ).get_double_side_set_cluster( aGhostSetupData.mDblSideSetIndexInMesh( iGhostSet ) );

            // iterate through ghost clusters
            for ( uint iGhostCluster = 0; iGhostCluster < tDblSideSetClusters.size(); iGhostCluster++ )
            {
                // get the leader and follower UIPCs
                mtk::Cell const & tLeaderIpCell   = tDblSideSetClusters( iGhostCluster )->get_interpolation_cell( mtk::Leader_Follower::LEADER );
                mtk::Cell const & tFollowerIpCell = tDblSideSetClusters( iGhostCluster )->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER );

                // get the vertices attached to leader/follower UIPCs
                Cell< mtk::Vertex* > tLeaderVertices   = tLeaderIpCell.get_vertex_pointers();
                Cell< mtk::Vertex* > tFollowerVertices = tFollowerIpCell.get_vertex_pointers();

                // iterate through leader vertices and place them in the correct list
                for ( uint iVertOnLeader = 0; iVertOnLeader < tLeaderVertices.size(); iVertOnLeader++ )
                {
                    // get the ID of the current unzipped vertex
                    moris_index tVertexId = tLeaderVertices( iVertOnLeader )->get_id();

                    // find out whether a T-matrix exists wrt. the current B-spline mesh
                    bool tHasInterpolation = tLeaderVertices( iVertOnLeader )->has_interpolation( aMeshIndex );

                    // add to list of vertices without interpolation if it hasn't been already
                    if ( !tHasInterpolation )
                    {
                        // check that it hasn't been added already to avoid unnecessary communication
                        if ( tGhostVerticesWithOutInterpolationMap.find( tVertexId ) == tGhostVerticesWithOutInterpolationMap.end() )
                        {
                            tGhostVerticesWithOutInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithoutInterpolation.push_back( tLeaderVertices( iVertOnLeader ) );
                            aGhostIpCellConnectedToVertex.push_back( &tLeaderIpCell );
                        }
                    }
                }    // end for: leader vertices

                // iterate through follower vertices and place them in the correct list
                for ( uint iVertOnFollower = 0; iVertOnFollower < tFollowerVertices.size(); iVertOnFollower++ )
                {
                    // get the ID of the current unzipped vertex
                    moris_index tVertexId = tFollowerVertices( iVertOnFollower )->get_id();

                    // find out whether a T-matrix exists wrt to the current B-spline mesh
                    bool tHasInterpolation = tFollowerVertices( iVertOnFollower )->has_interpolation( aMeshIndex );

                    // add to list of vertices without interpolation if it hasn't been already
                    if ( !tHasInterpolation )
                    {
                        // check that it hasn't been added already to avoid unnecessary communication
                        if ( tGhostVerticesWithOutInterpolationMap.find( tVertexId ) == tGhostVerticesWithOutInterpolationMap.end() )
                        {
                            tGhostVerticesWithOutInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithoutInterpolation.push_back( tFollowerVertices( iVertOnFollower ) );
                            aGhostIpCellConnectedToVertex.push_back( &tFollowerIpCell );
                        }
                    }
                }    // end for: follower vertices
            }        // end for: ghost clusters
        }            // end for: ghost sets

        // size out unused memory
        aGhostVerticesWithoutInterpolation.shrink_to_fit();
        aGhostIpCellConnectedToVertex.shrink_to_fit();
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_requests_for_T_matrices_without_interpolation(
            Cell< mtk::Vertex* > const &      aGhostVerticesWithoutInterpolation,
            Cell< mtk::Cell const * > const & aGhostIpCellConnectedToVertex,
            Cell< Cell< moris_index > >&      aNotOwnedIPVertIndsInNotOwnedList,
            Cell< Cell< moris_index > >&      aIPVertIndsToProcs,
            Cell< Matrix< IdMat > >&          aBaseVertexIds,
            Cell< Matrix< IdMat > >&          aUnzippedIpCellIds )
    {
        // get the communication table and map
        Matrix< IdMat >                   tCommTable              = mXTKModel->get_communication_table();
        uint                              tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = mXTKModel->get_communication_map();

        // initialize lists of identifying information
        aIPVertIndsToProcs.resize( tCommTableSize );
        aBaseVertexIds.resize( tCommTableSize );
        aUnzippedIpCellIds.resize( tCommTableSize );

        // get the number of UIPVs without a T-matrix on the executing processor
        uint tNumVertsWithoutTmat = aGhostVerticesWithoutInterpolation.size();

        // initialize map relating position of vertex in communication array back to position in not owned list
        aNotOwnedIPVertIndsInNotOwnedList.resize( tCommTableSize );

        // go through Vertices that executing proc knows about, but doesn't have a T-matrix for ...
        for ( uint iVertWithoutTmat = 0; iVertWithoutTmat < tNumVertsWithoutTmat; iVertWithoutTmat++ )
        {
            // ... get the owner of the UIPC the UIPV is attached to ...
            moris_index tOwnerProc = aGhostIpCellConnectedToVertex( iVertWithoutTmat )->get_owner();
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Ghost_Stabilization::prepare_requests_for_T_matrices_without_interpolation() - "
                    "Entity owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // ... get index of the UIPV ...
            moris_index tVertIndex = aGhostVerticesWithoutInterpolation( iVertWithoutTmat )->get_index();

            // ... and finally add the UIPV without a T-matrix in the list of entities to be requested from that owning proc
            aIPVertIndsToProcs( tProcDataIndex ).push_back( tVertIndex );

            // relate back position in comm arrays to position in the not owned list
            aNotOwnedIPVertIndsInNotOwnedList( tProcDataIndex ).push_back( iVertWithoutTmat );
        }

        // size out unused memory
        aIPVertIndsToProcs.shrink_to_fit();
        aNotOwnedIPVertIndsInNotOwnedList.shrink_to_fit();

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned entities to be sent to each processor processor
            uint tNumNotOwnedEntitiesOnProc = aIPVertIndsToProcs( iProc ).size();

            // allocate matrices
            aBaseVertexIds( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );
            aUnzippedIpCellIds( iProc ).resize( tNumNotOwnedEntitiesOnProc, 1 );

            // go through the entities for which IDs will be requested by the other processor
            for ( uint iVert = 0; iVert < tNumNotOwnedEntitiesOnProc; iVert++ )
            {
                // get the position of the information in the not-owned arrays
                moris_index tPosInNotOwnedArray = aNotOwnedIPVertIndsInNotOwnedList( iProc )( iVert );

                // get the UIPV on the executing proc and its base vertex ID
                mtk::Vertex const * tUIPV       = aGhostVerticesWithoutInterpolation( tPosInNotOwnedArray );
                moris_id            tBaseVertId = tUIPV->get_base_vertex()->get_id();

                // get the index and ID of the UIPC the current UIPV is attached to
                mtk::Cell const * tUIPC   = aGhostIpCellConnectedToVertex( tPosInNotOwnedArray );
                moris_id          tUipcId = tUIPC->get_id();

                // store the identifying information in the output arrays
                aBaseVertexIds( iProc )( iVert )     = tBaseVertId;
                aUnzippedIpCellIds( iProc )( iVert ) = tUipcId;
            }
        }    // end for: processors in comm-table
    }        // end function: Ghost_Stabilization::prepare_requests_for_T_matrices_without_interpolation()

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_answers_for_T_matrices(
            Cell< Matrix< IdMat > > const & aReceivedBaseVertexIds,
            Cell< Matrix< IdMat > > const & aReceivedUnzippedIpCellIds,
            Cell< Matrix< DDRMat > >&       aTMatrixWeights,
            Cell< Matrix< IdMat > >&        aTMatrixIds,
            Cell< Matrix< IdMat > >&        aTMatrixOwners,
            Cell< Matrix< IndexMat > >&     aTMatrixOffsets,
            const moris_index               aMeshIndex )
    {
        // access the enriched ip mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // get size of the communication table
        uint tCommTableSize = mXTKModel->get_communication_table().numel();

        // initialize communication arrays with correct size
        aTMatrixWeights.resize( tCommTableSize );
        aTMatrixIds.resize( tCommTableSize );
        aTMatrixOwners.resize( tCommTableSize );
        aTMatrixOffsets.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedBaseVertexIds.size() == tCommTableSize && aReceivedUnzippedIpCellIds.size() == tCommTableSize,
                "Ghost_Stabilization::prepare_answers_for_T_matrices() - Received information incomplete." );

        // initialize array storing pointers to the collected vertex interpolations
        Cell< Cell< Vertex_Enrichment const * > > tVertexInterpolations( tCommTableSize );

        // stores for each processor how big the linear arrays storing the T-matrix information will be
        Cell< moris_index > tDataSize( tCommTableSize, 0 );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of T-matrices requested from the current proc position
            uint tNumReceivedReqs = aReceivedBaseVertexIds( iProc ).numel();

            // initialize temporary array storing the enriched T-matrices found
            tVertexInterpolations( iProc ).resize( tNumReceivedReqs );

            // initialize information storing which index ranges in the communicated linear array belong to which T-matrices
            aTMatrixOffsets( iProc ).set_size( tNumReceivedReqs + 1, 1, -1 );
            aTMatrixOffsets( iProc )( 0 ) = 0;

            // iterate through the entities for which the IDs are requested
            for ( uint iVert = 0; iVert < tNumReceivedReqs; iVert++ )
            {
                // get the IDs from received information describing the entity
                moris_id tBaseVertexId = aReceivedBaseVertexIds( iProc )( iVert );
                moris_id tUipcId       = aReceivedUnzippedIpCellIds( iProc )( iVert );

                // get the UIPV
                moris_index                           tVertexIndex = this->get_enriched_interpolation_vertex( tBaseVertexId, tUipcId );
                Interpolation_Vertex_Unzipped const * tUIPV        = tEnrInterpMesh.get_unzipped_vertex_pointer( tVertexIndex );

                // get the enriched T-matrix
                Vertex_Enrichment const * tTMatrix = tUIPV->get_xtk_interpolation( aMeshIndex );

                // check that the enriched T-matrix is actually known
                MORIS_ASSERT( tTMatrix->get_base_vertex_interpolation() != nullptr,
                        "Ghost_Stabilization::prepare_answers_for_T_matrices() - "
                        "Owning (this) proc has a nullptr for the vertex interpolation on unzipped vertex #%i (ID: %i).",
                        tVertexIndex,
                        tUIPV->get_id() );

                // store pointer to the enriched T-matrix for later copy
                tVertexInterpolations( iProc )( iVert ) = tTMatrix;

                // check the size of the nodal T-matrix
                moris_index tNumBasisIdsOnTMat = tTMatrix->get_basis_indices().numel();

                // store the offsets in the linear data array
                aTMatrixOffsets( iProc )( iVert ) = tDataSize( iProc );

                // set the offset for the linear arrays using the size of the nodal T-matrices
                tDataSize( iProc ) = tDataSize( iProc ) + tNumBasisIdsOnTMat;

            }    // end for: each vertex on the current proc communicated with

            // store the size of the total array as the last offset
            aTMatrixOffsets( iProc )( tNumReceivedReqs ) = tDataSize( iProc );

        }    // end for: each proc communicated with

        // initialize return data size
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            aTMatrixWeights( iProc ).resize( tDataSize( iProc ), 1 );
            aTMatrixIds( iProc ).resize( tDataSize( iProc ), 1 );
            aTMatrixOwners( iProc ).resize( tDataSize( iProc ), 1 );
        }

        // populate the T-matrix data
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of T-matrices requested from the current proc position
            uint tNumReceivedReqs = aReceivedBaseVertexIds( iProc ).numel();

            // answer all of them
            for ( uint iVert = 0; iVert < tNumReceivedReqs; iVert++ )
            {
                // access the enriched T-matrix
                Vertex_Enrichment const * tEnrTMat = tVertexInterpolations( iProc )( iVert );

                // get the range of where the T-matrix information sits in the linear array
                moris_index tStart = aTMatrixOffsets( iProc )( iVert );
                moris_index tEnd   = aTMatrixOffsets( iProc )( iVert + 1 ) - 1;

                // fill the communication arrays
                aTMatrixWeights( iProc )( { tStart, tEnd }, { 0, 0 } ) = tEnrTMat->get_basis_weights().matrix_data();
                aTMatrixIds( iProc )( { tStart, tEnd }, { 0, 0 } )     = tEnrTMat->get_basis_ids().matrix_data();
                aTMatrixOwners( iProc )( { tStart, tEnd }, { 0, 0 } )  = tEnrTMat->get_owners().matrix_data();
            }
        }

        // size out unused memory
        aTMatrixWeights.shrink_to_fit();
        aTMatrixIds.shrink_to_fit();
        aTMatrixOwners.shrink_to_fit();
        aTMatrixOffsets.shrink_to_fit();

    }    // end function: Ghost_Stabilization::prepare_answers_for_T_matrices()

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::handle_requested_T_matrix_answers(
            Cell< mtk::Cell const * > const &   aGhostIpCellConnectedToVertex,
            Cell< Cell< moris_index > > const & aNotOwnedIPVertIndsInNotOwnedList,
            Cell< Cell< moris_index > > const & aIPVertIndsToProcs,
            Cell< Matrix< DDRMat > > const &    aReceivedTMatrixWeights,
            Cell< Matrix< IdMat > > const &     aReceivedTMatrixIds,
            Cell< Matrix< IdMat > > const &     aReceivedTMatrixOwners,
            Cell< Matrix< IndexMat > > const &  aReceivedTMatrixOffsets,
            const moris_index                   aMeshIndex )
    {
        // access the enriched ip mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // get the communication map
        std::map< moris_id, moris_index > tProcIdToCommTableIndex = mXTKModel->mCutIntegrationMesh->get_communication_map();

        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedTMatrixOffsets.size(); iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedEntityIds = aReceivedTMatrixOffsets( iProc ).numel() - 1;

            // make sure everything has been answered
            MORIS_ASSERT( tNumReceivedEntityIds == aIPVertIndsToProcs( iProc ).size(),
                    "Ghost_Stabilization::handle_requested_T_matrix_answers() - "
                    "Request arrays sent to and answers received from proc #%i have different size.",
                    mXTKModel->get_communication_table()( iProc ) );

            // assign IDs to each communicated entity
            for ( uint iVert = 0; iVert < tNumReceivedEntityIds; iVert++ )
            {
                // get the vertex
                moris_index                    tVertexIndex = aIPVertIndsToProcs( iProc )( iVert );
                Interpolation_Vertex_Unzipped& tVertex      = tEnrInterpMesh.get_xtk_interp_vertex( tVertexIndex );

                // get access to the enriched T-matrix
                Vertex_Enrichment* tEnrTMat = tVertex.get_xtk_interpolation( aMeshIndex );

                // get the index ranges and size of the T-matrices to received
                moris_index tStart = aReceivedTMatrixOffsets( iProc )( iVert );
                moris_index tEnd   = aReceivedTMatrixOffsets( iProc )( iVert + 1 ) - 1;
                uint        tSize  = (uint)( tEnd - tStart + 1 );

                // extract the information from the linear communication arrays
                Matrix< IdMat >  tBasisIds     = aReceivedTMatrixIds( iProc )( { tStart, tEnd }, { 0, 0 } );
                Matrix< IdMat >  tBasisOwners  = aReceivedTMatrixOwners( iProc )( { tStart, tEnd }, { 0, 0 } );
                Matrix< DDRMat > tBasisWeights = aReceivedTMatrixWeights( iProc )( { tStart, tEnd }, { 0, 0 } );

                // collect the indices for the basis IDs received
                Matrix< IndexMat > tBasisIndices( tSize, 1, -1 );
                for ( uint iBF = 0; iBF < tSize; iBF++ )
                {
                    // check that this basis is indeed owned by another proc
                    moris_id tBasisOwner = tBasisOwners( iBF );

                    // get the current BF's ID
                    moris_id tBasisId = tBasisIds( iBF );

                    // add this basis to the mesh if it does not exists on the current partition
                    if ( !tEnrInterpMesh.basis_exists_on_partition( aMeshIndex, tBasisId ) )
                    {
                        // get the bulk-phase the basis interpolates into
                        moris_index                  tUipcIndexInNotOwnedData = aNotOwnedIPVertIndsInNotOwnedList( iProc )( iVert );
                        moris_index                  tUipcIndex               = aGhostIpCellConnectedToVertex( tUipcIndexInNotOwnedData )->get_index();
                        Interpolation_Cell_Unzipped* tEnrIpCell               = tEnrInterpMesh.get_enriched_interpolation_cells()( tUipcIndex );
                        moris_index                  tBulkPhase               = tEnrIpCell->get_bulkphase_index();

                        // add basis ID to partition
                        tEnrInterpMesh.add_basis_function( aMeshIndex, tBasisId, tBasisOwner, tBulkPhase );
                    }

                    // find and store the basis index local to the executing processor
                    tBasisIndices( iBF ) = tEnrInterpMesh.get_enr_basis_index_from_enr_basis_id( aMeshIndex, tBasisId );

                    // if the basis has an owning proc that is not in the comm table, add it to the comm table
                    if ( tProcIdToCommTableIndex.find( tBasisOwner ) == tProcIdToCommTableIndex.end() && tBasisOwner != par_rank() )
                    {
                        tEnrInterpMesh.add_proc_to_comm_table( tBasisOwner );
                        tProcIdToCommTableIndex = mXTKModel->mCutIntegrationMesh->get_communication_map();
                    }
                }

                // Setup the basis index to ID map for the T-matrix
                IndexMap& tVertEnrichMap = tEnrTMat->get_basis_map();
                for ( uint iBF = 0; iBF < tSize; iBF++ )
                {
                    moris_index tBasisIndex       = tBasisIndices( iBF );
                    tVertEnrichMap[ tBasisIndex ] = iBF;
                }

                // get the basis indices from the basis ids
                tEnrTMat->add_basis_information( tBasisIndices, tBasisIds );
                tEnrTMat->add_basis_weights( tBasisIndices, tBasisWeights );
                tEnrTMat->add_basis_owners( tBasisIndices, tBasisOwners );
                tEnrTMat->add_base_vertex_interpolation( nullptr );
                // base vertex interpolation does not exists (other  proc)

            }    // end for: each vertex communicated with current proc

        }        // end for: each proc answers are received from

    }            // end function: Ghost_Stabilization::handle_requested_T_matrix_answers()

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::get_ip_vertices_in_ghost_sets(
            Ghost_Setup_Data&                 aGhostSetupData,
            moris::Cell< mtk::Vertex* >&      aGhostVerticesWithInterpolation,
            moris::Cell< mtk::Vertex* >&      aGhostVerticesWithoutInterpolation,
            moris::Cell< mtk::Cell const * >& aGhostIpCellConnectedToVertex )
    {
        uint tNumGhostSets = aGhostSetupData.mDblSideSetIndexInMesh.size();

        std::unordered_map< moris_index, bool > tGhostVerticesWithInterpolationMap;
        std::unordered_map< moris_index, bool > tGhostVerticesWithOutInterpolationMap;


        // iterate through vertices and gather the
        for ( uint iS = 0; iS < tNumGhostSets; iS++ )
        {
            // get the clusters
            moris::Cell< mtk::Cluster const * > tDblSideSetClusters =
                    mXTKModel->get_enriched_integ_mesh( 0 ).get_double_side_set_cluster( aGhostSetupData.mDblSideSetIndexInMesh( iS ) );

            // iterate through clusters
            for ( uint iC = 0; iC < tDblSideSetClusters.size(); iC++ )
            {
                // get the leader ip cell
                moris::mtk::Cell const & tLeaderIpCell = tDblSideSetClusters( iC )->get_interpolation_cell( mtk::Leader_Follower::LEADER );

                // get the follower ip cell
                moris::mtk::Cell const & tFollowerIpCell = tDblSideSetClusters( iC )->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER );

                // get the vertices attached to leader/follower cells
                moris::Cell< mtk::Vertex* > tLeaderVertices   = tLeaderIpCell.get_vertex_pointers();
                moris::Cell< mtk::Vertex* > tFollowerVertices = tFollowerIpCell.get_vertex_pointers();

                // iterate through leader vertices and place them in the correct list
                for ( uint iV = 0; iV < tLeaderVertices.size(); iV++ )
                {
                    // get the ID of the current unzipped vertex
                    moris_index tVertexId = tLeaderVertices( iV )->get_id();

                    // find out whether a T-matrix exists wrt. each B-spline mesh
                    bool tHasInterpolation = true;
                    for ( uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
                    {
                        tHasInterpolation = ( tHasInterpolation && tLeaderVertices( iV )->has_interpolation( mMeshIndices( iBspMesh ) ) );
                    }

                    // add to vertices without interpolation
                    if ( !tHasInterpolation )
                    {
                        if ( tGhostVerticesWithOutInterpolationMap.find( tVertexId ) == tGhostVerticesWithOutInterpolationMap.end() )
                        {
                            tGhostVerticesWithOutInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithoutInterpolation.push_back( tLeaderVertices( iV ) );
                            aGhostIpCellConnectedToVertex.push_back( &tLeaderIpCell );
                        }
                    }

                    // add to vertices with interpolation
                    else if ( tHasInterpolation )
                    {
                        if ( tGhostVerticesWithInterpolationMap.find( tVertexId ) == tGhostVerticesWithInterpolationMap.end() )
                        {
                            tGhostVerticesWithInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithInterpolation.push_back( tLeaderVertices( iV ) );
                        }
                    }
                }

                // iterate through follower vertices and place them in the correct list
                for ( uint iV = 0; iV < tFollowerVertices.size(); iV++ )
                {
                    // get the ID of the current unzipped vertex
                    moris_index tVertexId = tFollowerVertices( iV )->get_id();

                    // find out whether a T-matrix exists wrt. each B-spline mesh
                    bool tHasInterpolation = true;
                    for ( uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
                    {
                        tHasInterpolation = ( tHasInterpolation && tFollowerVertices( iV )->has_interpolation( mMeshIndices( iBspMesh ) ) );
                    }

                    // add to vertices without interpolation
                    if ( !tHasInterpolation )
                    {
                        if ( tGhostVerticesWithOutInterpolationMap.find( tVertexId ) == tGhostVerticesWithOutInterpolationMap.end() )
                        {
                            tGhostVerticesWithOutInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithoutInterpolation.push_back( tFollowerVertices( iV ) );
                            aGhostIpCellConnectedToVertex.push_back( &tFollowerIpCell );
                        }
                    }

                    // add to vertices with interpolation
                    else if ( tHasInterpolation )
                    {
                        if ( tGhostVerticesWithInterpolationMap.find( tVertexId ) == tGhostVerticesWithInterpolationMap.end() )
                        {
                            tGhostVerticesWithInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithInterpolation.push_back( tFollowerVertices( iV ) );
                        }
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Ghost_Stabilization::get_enriched_interpolation_vertex(
            moris_index const & aBGVertId,
            moris_index const & aEnrichedIpCellId )
    {
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // get the interpolation cell index using the id
        moris_index tCellIndex = tEnrInterpMesh.get_loc_entity_ind_from_entity_glb_id( aEnrichedIpCellId, mtk::EntityRank::ELEMENT, 0 );

        // get the cell
        Interpolation_Cell_Unzipped* tEnrIpCell = tEnrInterpMesh.get_enriched_interpolation_cells()( tCellIndex );

        // get the vertices
        moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVertexPointers = tEnrIpCell->get_xtk_interpolation_vertices();

        moris_index tVertexPointerInd = 0;
        uint        tCount            = 0;

        for ( uint i = 0; i < tVertexPointers.size(); i++ )
        {
            if ( tVertexPointers( i )->get_base_vertex()->get_id() == aBGVertId )
            {
                tVertexPointerInd = tVertexPointers( i )->get_index();
                tCount++;
            }
        }

        MORIS_ERROR( tCount == 1, "Enriched interpolation vertex not found or found more than once" );

        return tVertexPointerInd;
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::declare_ghost_double_side_sets_in_mesh( Ghost_Setup_Data& aGhostSetupData )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Ghost", "declare_ghost_double_side_sets_in_mesh" );

        // get number of Bulk phases on mesh
        uint tNumBulkPhases = mXTKModel->get_geom_engine()->get_num_bulk_phase();

        // initialize container with Ghost Set Names
        Cell< std::string > tGhostDoubleSideNames( tNumBulkPhases );

        // for every bulk phase ...
        for ( moris::moris_index iBulkPhase = 0; iBulkPhase < (moris_index)tNumBulkPhases; iBulkPhase++ )
        {
            // .. create name for ghost set associated with current bulk phase and save it in container
            std::string tGhostSideSetName       = this->get_ghost_dbl_side_set_name( iBulkPhase );
            tGhostDoubleSideNames( iBulkPhase ) = tGhostSideSetName;
        }

        // register and get the dbl. side set indices for the ghost side sets in the overall mesh data
        aGhostSetupData.mDblSideSetIndexInMesh = mXTKModel->get_enriched_integ_mesh( 0 ).register_double_side_set_names( tGhostDoubleSideNames );
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::declare_ghost_double_side_sets_in_mesh_new( Ghost_Setup_Data& aGhostSetupData )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Ghost", "declare_ghost_double_side_sets_in_mesh" );

        // get number of Bulk phases on mesh
        uint tNumBulkPhases = mXTKModel->get_geom_engine()->get_num_bulk_phase();

        // get number of B-spline meshes
        uint tNumBspMeshes = mMeshIndices.numel();

        // initialize container with Ghost Set Names
        Cell< std::string > tGhostDoubleSideNames( tNumBulkPhases * tNumBspMeshes );

        // initialize counter for Ghost side sets
        uint tNumGhostSideSets = 0;

        // for every B-spline mesh index ...
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the B-spline mesh index
            moris_index tBspMeshIndex = mMeshIndices( iBspMesh );

            // ... and every bulk phase ...
            for ( moris::moris_index iBulkPhase = 0; iBulkPhase < (moris_index)tNumBulkPhases; iBulkPhase++ )
            {
                // .. create name for ghost set associated with current bulk phase and B-spline mesh and save it in container
                std::string tGhostSideSetName              = this->get_ghost_dbl_side_set_name( iBulkPhase, tBspMeshIndex );
                tGhostDoubleSideNames( tNumGhostSideSets ) = tGhostSideSetName;

                // count number of Ghost side sets
                tNumGhostSideSets++;
            }
        }

        // register and get the dbl. side set indices for the ghost side sets in the overall mesh data
        aGhostSetupData.mDblSideSetIndexInMesh = mXTKModel->get_enriched_integ_mesh( 0 ).register_double_side_set_names( tGhostDoubleSideNames );
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh( Ghost_Setup_Data& aGhostSetupData )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Ghost", "construct_ghost_double_side_sets_in_mesh" );

        // get the enriched IP and IG meshes
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();
        Enriched_Integration_Mesh&   tEnrIntegMesh  = mXTKModel->get_enriched_integ_mesh();

        // get list of all enr. (i.e. unzipped) IP cells
        Cell< Interpolation_Cell_Unzipped* >& tEnrIpCells = tEnrInterpMesh.get_enriched_interpolation_cells();
        Cut_Integration_Mesh*                 tCutIgMesh  = mXTKModel->get_cut_integration_mesh();

        // access subphase neighborhood information
        std::shared_ptr< Subphase_Neighborhood_Connectivity >                tSubphaseNeighborhood               = tCutIgMesh->get_subphase_neighborhood();
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphaseToSubphase                 = tSubphaseNeighborhood->mSubphaseToSubPhase;
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphaseToSubphaseMySideOrds       = tSubphaseNeighborhood->mSubphaseToSubPhaseMySideOrds;
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSubphaseToSubphaseNeighborSideOrds = tSubphaseNeighborhood->mSubphaseToSubPhaseNeighborSideOrds;
        moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tTransitionLocation                 = tSubphaseNeighborhood->mTransitionNeighborCellLocation;

        // get number of bulk phases in the mesh
        uint tNumBulkPhases = mXTKModel->get_geom_engine()->get_num_bulk_phase();

        // number of subphases in mesh
        uint tNumSubphases = tSubphaseToSubphase.size();

        // reserve space in data
        const uint tReserveDim = 10;

        // estimate: number of ghost facets is in the same order of magnitude as the number of subphases
        uint tReserveSize = tNumBulkPhases * std::min( tReserveDim, tNumSubphases );

        // reserve memory in ghost setup data according to estimate
        aGhostSetupData.mLeaderSideIpCells.reserve( tReserveSize );
        aGhostSetupData.mFollowerSideIpCells.reserve( tReserveSize );
        aGhostSetupData.mLeaderSideIgCellSideOrds.reserve( tReserveSize );
        aGhostSetupData.mFollowerSideIgCellSideOrds.reserve( tReserveSize );
        aGhostSetupData.mNonTrivialFlag.reserve( tReserveSize );
        aGhostSetupData.mTransitionLocation.reserve( tReserveSize );

        aGhostSetupData.mLeaderSideIpCells.resize( tNumBulkPhases );
        aGhostSetupData.mFollowerSideIpCells.resize( tNumBulkPhases );
        aGhostSetupData.mLeaderSideIgCellSideOrds.resize( tNumBulkPhases );
        aGhostSetupData.mFollowerSideIgCellSideOrds.resize( tNumBulkPhases );
        aGhostSetupData.mNonTrivialFlag.resize( tNumBulkPhases );
        aGhostSetupData.mTransitionLocation.resize( tNumBulkPhases );

        // for debug
        // mXTKModel->print_subphase_neighborhood();

        // initialize flag indicating whether its trivial or not
        moris_index tTrivial         = 0;
        moris_index tNonTrivialCount = 0;

        // ----------------------------------------------------------------------------------
        // loop collecting ghost side sets to be constructed
        // iterate through subphases
        for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
        {
            // iterate through subphase neighbors
            for ( uint jSpNeighbor = 0; jSpNeighbor < tSubphaseToSubphase( iSP )->size(); jSpNeighbor++ )
            {
                // if I am the one constructing this subphase then add it to ghost setup data
                if ( this->create_ghost( aGhostSetupData, (moris_index)iSP, ( *tSubphaseToSubphase( iSP ) )( jSpNeighbor ), tTrivial ) )
                {
                    moris_index tFirstInterpCellIndex =
                            aGhostSetupData.mSubphaseIndexToInterpolationCellIndex( iSP );

                    moris_index tSecondInterpCellIndex =
                            aGhostSetupData.mSubphaseIndexToInterpolationCellIndex( ( *tSubphaseToSubphase( iSP ) )( jSpNeighbor ) );

                    Interpolation_Cell_Unzipped* tFirstInterpCell  = tEnrIpCells( tFirstInterpCellIndex );
                    Interpolation_Cell_Unzipped* tSecondInterpCell = tEnrIpCells( tSecondInterpCellIndex );

                    // get the bulk phase
                    moris_index tBulkPhase = tFirstInterpCell->get_bulkphase_index();

                    MORIS_ASSERT( tBulkPhase == tSecondInterpCell->get_bulkphase_index(),
                            "Bulk phase between neighboring subphases does not match" );

                    // setup ip cell indices in ghost setup data
                    aGhostSetupData.mLeaderSideIpCells( tBulkPhase ).push_back( tFirstInterpCell->get_index() );
                    aGhostSetupData.mFollowerSideIpCells( tBulkPhase ).push_back( tSecondInterpCell->get_index() );

                    // setup ig cells in ghost set up data
                    aGhostSetupData.mLeaderSideIgCellSideOrds( tBulkPhase ).push_back( ( *tSubphaseToSubphaseMySideOrds( iSP ) )( jSpNeighbor ) );
                    aGhostSetupData.mFollowerSideIgCellSideOrds( tBulkPhase ).push_back( ( *tSubphaseToSubphaseNeighborSideOrds( iSP ) )( jSpNeighbor ) );

                    // mark as trivial or not in ghost setup data
                    aGhostSetupData.mNonTrivialFlag( tBulkPhase ).push_back( tTrivial );

                    // mark the transition location
                    aGhostSetupData.mTransitionLocation( tBulkPhase ).push_back( ( *tTransitionLocation( iSP ) )( jSpNeighbor ) );

                    // increment count
                    if ( tTrivial > 0 )
                    {
                        tNonTrivialCount++;
                    }
                }
            }
        }

        // check that reserved size was appropriate
        MORIS_ASSERT( aGhostSetupData.mLeaderSideIpCells.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mLeaderSideIpCells too small, increase by %zu\n",
                aGhostSetupData.mLeaderSideIpCells.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mFollowerSideIpCells.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mFollowerSideIpCells too small, increase by %zu\n",
                aGhostSetupData.mFollowerSideIpCells.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mLeaderSideIgCellSideOrds.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mLeaderSideIgCellSideOrds too small, increase by %zu\n",
                aGhostSetupData.mLeaderSideIgCellSideOrds.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mFollowerSideIgCellSideOrds.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mFollowerSideIgCellSideOrds too small, increase by %zu\n",
                aGhostSetupData.mFollowerSideIgCellSideOrds.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mNonTrivialFlag.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mTrivialFlag too small, increase by %zu\n",
                aGhostSetupData.mNonTrivialFlag.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mTransitionLocation.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mTransitionLocation too small, increase by %zu\n",
                aGhostSetupData.mTransitionLocation.size() / tReserveSize );

        // free up unused memory
        aGhostSetupData.mLeaderSideIpCells.shrink_to_fit();
        aGhostSetupData.mFollowerSideIpCells.shrink_to_fit();
        aGhostSetupData.mLeaderSideIgCellSideOrds.shrink_to_fit();
        aGhostSetupData.mFollowerSideIgCellSideOrds.shrink_to_fit();
        aGhostSetupData.mNonTrivialFlag.shrink_to_fit();
        aGhostSetupData.mTransitionLocation.shrink_to_fit();

        // allocate ids for non-trivial integration cells
        moris_id tCurrentId = 0;

        // reserve global entity ids for additional elements that need to be created to get dbl. sided facets ...
        if ( aGhostSetupData.mLinearBackgroundMesh )
        {
            // ... for non-trivial element transitions only in the linear case
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tNonTrivialCount, mtk::EntityRank::ELEMENT );
        }
        else
        {
            // ... // TODO: I don't understand this bit
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tEnrIntegMesh.get_num_entities( mtk::EntityRank::ELEMENT ), mtk::EntityRank::ELEMENT );
        }

        moris_id tCurrentIndex = tEnrIntegMesh.get_num_entities( mtk::EntityRank::ELEMENT );

        // get total number of ghost facets for each bulk phase on current processor
        Matrix< DDUMat > tLocalNumberOfGhostFacets( aGhostSetupData.mLeaderSideIpCells.size(), 1 );
        for ( uint iBulkPhase = 0; iBulkPhase < aGhostSetupData.mLeaderSideIpCells.size(); iBulkPhase++ )
        {
            tLocalNumberOfGhostFacets( iBulkPhase ) = aGhostSetupData.mLeaderSideIpCells( iBulkPhase ).size();
        }

        // get global (i.e. across all procs) number of ghost facets for each bulk phase
        Matrix< DDUMat > tTotalNumberOfGhostFacets = sum_all_matrix( tLocalNumberOfGhostFacets );

        // build list of double side set indices
        moris::Cell< moris_index > tDoubleSideSetIndexList;
        tDoubleSideSetIndexList.reserve( aGhostSetupData.mLeaderSideIpCells.size() );

        // ----------------------------------------------------------------------------------
        // loop constructing the actual ghost facets
        // iterate through bulk phases
        for ( uint iBulkPhase = 0; iBulkPhase < aGhostSetupData.mLeaderSideIpCells.size(); iBulkPhase++ )
        {
            // allocate space in the integration mesh double side sets
            tEnrIntegMesh.mDoubleSideSets( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) ).resize( aGhostSetupData.mLeaderSideIpCells( iBulkPhase ).size() );

            MORIS_LOG_SPEC( "Total Ghost Facets for Bulk Phase " + std::to_string( iBulkPhase ), tTotalNumberOfGhostFacets( iBulkPhase ) );

            // iterate through double sides in this bulk phase
            for ( uint iGhostFacet = 0; iGhostFacet < aGhostSetupData.mLeaderSideIpCells( iBulkPhase ).size(); iGhostFacet++ )
            {
                // create a new side cluster for each of the pairs
                std::shared_ptr< Side_Cluster > tFollowerSideCluster =
                        this->create_follower_side_cluster( aGhostSetupData, tEnrIpCells, iBulkPhase, iGhostFacet, tCurrentIndex, tCurrentId );

                std::shared_ptr< Side_Cluster > tLeaderSideCluster =
                        this->create_leader_side_cluster( aGhostSetupData, tEnrIpCells, iBulkPhase, iGhostFacet, tFollowerSideCluster.get(), tCurrentIndex, tCurrentId );

                // verify the subphase cluster
                MORIS_ASSERT( tFollowerSideCluster->mInterpolationCell->get_bulkphase_index() == (moris_index)iBulkPhase,
                        "Bulk phase mismatch on follower side of double side set cluster" );

                MORIS_ASSERT( tLeaderSideCluster->mInterpolationCell->get_bulkphase_index() == (moris_index)iBulkPhase,
                        "Bulk phase mismatch on leader side of double side set cluster" );

                // add to side clusters the integration mesh
                tEnrIntegMesh.mDoubleSideSetsLeaderIndex( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() );

                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tLeaderSideCluster );

                tEnrIntegMesh.mDoubleSideSetsFollowerIndex( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() );

                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tFollowerSideCluster );

                // create double side cluster
                std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster = std::make_shared< mtk::Double_Side_Cluster >(
                        tLeaderSideCluster.get(),
                        tFollowerSideCluster.get(),
                        tLeaderSideCluster->mVerticesInCluster );

                // add to integration mesh
                tEnrIntegMesh.mDoubleSideClusters.push_back( tDblSideCluster );

                // add to the integration mesh
                tEnrIntegMesh.mDoubleSideSets( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) )( iGhostFacet ) = tDblSideCluster;
            }

            tDoubleSideSetIndexList.push_back( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) );

            tEnrIntegMesh.set_double_side_set_colors(
                    aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ),
                    { { (moris_index)iBulkPhase } },
                    { { (moris_index)iBulkPhase } } );
        }

        // commit and communicate the double side sets
        tEnrIntegMesh.commit_double_side_set( tDoubleSideSetIndexList );
        tEnrIntegMesh.communicate_sets_of_type( mtk::SetType::DOUBLE_SIDED_SIDESET );

        tEnrIntegMesh.collect_all_sets();
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh_new( Ghost_Setup_Data& aGhostSetupData )
    {
        // log/trace this function
        Tracer tTracer( "XTK", "Ghost", "Construct Ghost double side sets (new approach)" );

        // get the enriched IP and IG meshes
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();
        Enriched_Integration_Mesh&   tEnrIntegMesh  = mXTKModel->get_enriched_integ_mesh();

        // get the background mesh
        moris::mtk::Mesh* tLagrangeBackgroundMesh = &mXTKModel->get_background_mesh();

        // get list of all enr. (i.e. unzipped) IP cells
        Cell< Interpolation_Cell_Unzipped* >& tEnrIpCells = tEnrInterpMesh.get_enriched_interpolation_cells();
        Cut_Integration_Mesh*                 tCutIgMesh  = mXTKModel->get_cut_integration_mesh();

        // get number of bulk phases in the mesh
        uint tNumBulkPhases = mXTKModel->get_geom_engine()->get_num_bulk_phase();

        // estimate: number of ghost facets is in the same order of magnitude as the number of subphases
        // FIXME: estimate base on SPs, should be based on SPGs for better accuracy
        uint tReserveSize = tNumBulkPhases * std::min( (uint)10, tCutIgMesh->get_num_subphases() );

        // get number of B-spline meshes
        uint tNumBspMeshes = mMeshIndices.numel();

        // size lists in Ghost setup data
        aGhostSetupData.mLeaderSideIpCellsNew.resize( tNumBspMeshes );
        aGhostSetupData.mFollowerSideIpCellsNew.resize( tNumBspMeshes );
        aGhostSetupData.mLeaderSideIgCellSideOrdsNew.resize( tNumBspMeshes );
        aGhostSetupData.mFollowerSideIgCellSideOrdsNew.resize( tNumBspMeshes );
        aGhostSetupData.mNonTrivialFlagNew.resize( tNumBspMeshes );
        aGhostSetupData.mTransitionLocationNew.resize( tNumBspMeshes );

        // initialize counter keeping track of the number of non-trivial Lagrange element transitions
        uint tNonTrivialCount = 0;

        // ----------------------------------------------------------------------------------
        // loops collecting ghost facets to be constructed
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // get the mesh index
            moris_index tBspMeshIndex = mMeshIndices( iBspMesh );

            // access subphase group neighborhood information
            std::shared_ptr< Subphase_Neighborhood_Connectivity >                tSpgNeighborhood          = tCutIgMesh->get_subphase_group_neighborhood( iBspMesh );
            moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSpgToSpg                 = tSpgNeighborhood->mSubphaseToSubPhase;
            moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSpgToSpgMySideOrds       = tSpgNeighborhood->mSubphaseToSubPhaseMySideOrds;
            moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const & tSpgToSpgNeighborSideOrds = tSpgNeighborhood->mSubphaseToSubPhaseNeighborSideOrds;
            // moris::Cell< std::shared_ptr< moris::Cell< moris_index > > > const& tTransitionLocation       = tSpgNeighborhood->mTransitionNeighborCellLocation;

            // reserve memory in ghost setup data according to estimate
            aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mFollowerSideIpCellsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mLeaderSideIgCellSideOrdsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mFollowerSideIgCellSideOrdsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mNonTrivialFlagNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mTransitionLocationNew( iBspMesh ).reserve( tReserveSize );

            aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mFollowerSideIpCellsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mLeaderSideIgCellSideOrdsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mFollowerSideIgCellSideOrdsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mNonTrivialFlagNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mTransitionLocationNew( iBspMesh ).resize( tNumBulkPhases );

            // quick access to SPG list for current Bspline mesh
            moris::Cell< Subphase_Group* > const & tSubphaseGroups = mBsplineMeshInfos( iBspMesh )->mSubphaseGroups;
            uint                                   tNumSPGs        = tSubphaseGroups.size();

            // sanity check
            MORIS_ASSERT( tSpgToSpg.size() == tNumSPGs,
                    "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh_new() - "
                    "Number of SPGs in SPG neighborhood and in B-spline mesh info don't match." );

            // get the relationship between the unzipped IP cells, the base IP cells and SPGs from the enrichment
            moris::Cell< moris::Cell< moris_index > > const & tBaseIpCellAndSpgToUnzipping =
                    mXTKModel->mEnrichment->get_unzipped_IP_cell_to_base_IP_cell_and_SPG( iBspMesh );

            // loop collecting ghost side sets to be constructed for the current mesh
            // go through all SPG - neighbor pairs
            for ( uint iSPG = 0; iSPG < tNumSPGs; iSPG++ )
            {
                // get the Bulk-phase and B-spline element indices for the current SPG
                moris_index tBulkPhaseIndex       = tSubphaseGroups( iSPG )->get_bulk_phase();
                moris_index tLeaderBspElemIndex   = tSubphaseGroups( iSPG )->get_bspline_cell_index();
                moris_index tLeaderLocalSpgIndex  = tSubphaseGroups( iSPG )->get_local_index();
                moris_index tLeaderRepIpCellIndex = mBsplineMeshInfos( iBspMesh )->mExtractionCellsIndicesInBsplineCells( tLeaderBspElemIndex )( 0 );

                // access SPGs that are connected to the current SPG
                moris::Cell< moris_index > const & tNeighborSPGs     = *tSpgToSpg( iSPG );
                moris::Cell< moris_index > const & tMySideOrds       = *tSpgToSpgMySideOrds( iSPG );
                moris::Cell< moris_index > const & tNeighborSideOrds = *tSpgToSpgNeighborSideOrds( iSPG );

                uint tNumNeighborSPGs = tNeighborSPGs.size();

                for ( uint iNeighborSPG = 0; iNeighborSPG < tNumNeighborSPGs; iNeighborSPG++ )
                {
                    // get neighbor SPG's index
                    moris_index tNeighborSpgIndex = tNeighborSPGs( iNeighborSPG );

                    // get the side ordinal of the leader cell for the transition
                    moris_index tMySideOrdinal = tMySideOrds( iNeighborSPG );

                    // sanity check that both SPGs have the same bulk-phase indices
                    MORIS_ASSERT( tBulkPhaseIndex == tSubphaseGroups( tNeighborSpgIndex )->get_bulk_phase(),
                            "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh_new() - "
                            "Bulk phase between neighboring subphase groups does not match" );

                    // check whether ghost facts need to be created for this SPG pair
                    moris_index tIsNonTrivial      = 0;
                    bool        tCreateGhostFacets = this->create_ghost_new( iBspMesh, iSPG, tNeighborSpgIndex, tIsNonTrivial );

                    // mark facets for Ghost construction
                    if ( tCreateGhostFacets )
                    {
                        // initialize lists of base Ip cells on the leader B-spline element
                        moris::Cell< mtk::Cell* > tBaseIpCellsOnMyOrdinal;

                        // get the base IP elements on the first B-spline elements' facet
                        // Note: HMR uses 1-based indices for the side ordinals, hence the "+1"
                        tLagrangeBackgroundMesh->get_elements_in_interpolation_cluster_and_side_ordinal(
                                tLeaderRepIpCellIndex, tBspMeshIndex, tMySideOrdinal + 1, tBaseIpCellsOnMyOrdinal );

                        // get SPG index local to the neighbor B-spline element
                        moris_index tSecondLocalSpgIndex = tSubphaseGroups( tNeighborSpgIndex )->get_local_index();

                        // get the side ordinal from the neighbor's point of view
                        moris_index tNeighborSideOrdinal = tNeighborSideOrds( iNeighborSPG );

                        // go over all IP cells on the side of the leader B-spline element
                        for ( uint iMyIpCells = 0; iMyIpCells < tBaseIpCellsOnMyOrdinal.size(); iMyIpCells++ )
                        {
                            // current Base IP cell index currently treated
                            moris_index tBaseIpCellIndex = tBaseIpCellsOnMyOrdinal( iMyIpCells )->get_index();

                            // get the UIPC corresponding to the current base IP cell on the given bulk phase and SPG
                            moris_index tUipcIndex = tBaseIpCellAndSpgToUnzipping( tBaseIpCellIndex )( tLeaderLocalSpgIndex );

                            // refinement level of the current base IP cell (to be filled)
                            moris_index tMyRefineLevel = -1;

                            // list of information describing the transition
                            moris::Cell< moris_index > tNeighborElements;
                            moris::Cell< moris_index > tNeighborSideOrdinals;
                            moris::Cell< moris_index > tTransitionLocations;
                            moris::Cell< moris_index > tNeighborRefinementLevels;

                            // get the IP cells the current cell is connected to and whether it's a non-trivial element transition
                            bool tNonTrivialFlag = tLagrangeBackgroundMesh->get_elements_connected_to_element_through_face_ord(
                                    tBaseIpCellIndex,
                                    tMySideOrdinal,
                                    tMyRefineLevel,
                                    tNeighborElements,
                                    tNeighborSideOrdinals,
                                    tTransitionLocations,
                                    tNeighborRefinementLevels );

                            // get the number of Lagrange elements the current Lagrange cell is connected to
                            uint tNumConnectedIpCells = tNeighborElements.size();

                            // check if the current element is the small one,
                            // if so, assume the bigger neighbor to be the leader
                            if ( tMyRefineLevel > tNeighborRefinementLevels( 0 ) )
                            {
                                // check that there's only one neighbor in this case
                                MORIS_ASSERT( tNeighborElements.size() == 1,
                                        "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh_new() - "
                                        "Transition from small to big with more than one neighbor. This is not possible." );

                                // get the neighbor base IP cell index
                                moris_index tNeighborBaseIpCellIndex = tNeighborElements( 0 );

                                // get the neighbor UIPC index
                                moris_index tNeighborUipcIndex = tBaseIpCellAndSpgToUnzipping( tNeighborBaseIpCellIndex )( tSecondLocalSpgIndex );

                                // store in Ghost setup data the Follower UIPCs used and their side ordinals
                                aGhostSetupData.mFollowerSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tUipcIndex );
                                aGhostSetupData.mFollowerSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tMySideOrdinal );

                                // store in Ghost setup data the Leader UIPCs used and their side ordinals
                                aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborUipcIndex );
                                aGhostSetupData.mLeaderSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborSideOrdinal );

                                // store the transition location and that the transition is non-trivial
                                aGhostSetupData.mTransitionLocationNew( iBspMesh )( tBulkPhaseIndex ).push_back( tTransitionLocations( 0 ) );
                                aGhostSetupData.mNonTrivialFlagNew( iBspMesh )( tBulkPhaseIndex ).push_back( (moris_index)tNonTrivialFlag );

                                // count number of non-trivial element transitions
                                tNonTrivialCount++;
                            }

                            // the neighbor is of same level, or finer
                            else
                            {
                                // go over the Lagrange elements the current Lagrange element is connected to and store element pairs in the Ghost setup data
                                for ( uint iNeighborIpCell = 0; iNeighborIpCell < tNumConnectedIpCells; iNeighborIpCell++ )
                                {
                                    // get the neighbor base IP cell index
                                    moris_index tNeighborBaseIpCellIndex = tNeighborElements( iNeighborIpCell );

                                    // get the neighbor UIPC index
                                    moris_index tNeighborUipcIndex = tBaseIpCellAndSpgToUnzipping( tNeighborBaseIpCellIndex )( tSecondLocalSpgIndex );

                                    // store in Ghost setup data the Follower UIPCs used and their side ordinals
                                    aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tUipcIndex );
                                    aGhostSetupData.mLeaderSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tMySideOrdinal );

                                    // store in Ghost setup data the Leader UIPCs used and their side ordinals
                                    aGhostSetupData.mFollowerSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborUipcIndex );
                                    aGhostSetupData.mFollowerSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborSideOrdinal );

                                    // store the transition location and that the transition is non-trivial
                                    aGhostSetupData.mTransitionLocationNew( iBspMesh )( tBulkPhaseIndex ).push_back( tTransitionLocations( iNeighborIpCell ) );
                                    aGhostSetupData.mNonTrivialFlagNew( iBspMesh )( tBulkPhaseIndex ).push_back( (moris_index)tNonTrivialFlag );

                                    // count number of non-trivial element transitions
                                    if ( tNonTrivialFlag )
                                    {
                                        tNonTrivialCount++;
                                    }
                                }
                            }
                        }    // end for: loop over IP cells connected to leader B-spline cell facet

                    }        // end if: mark elements for ghost side construction
                }            // end for: loop over neighboring SPGs
            }                // end for: loop over all SPGs

            // check that reserved size was appropriate
            MORIS_ASSERT( aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh ).size() < tReserveSize,
                    "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of aGhostSetupData too small, increased by %zu\n",
                    aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh ).size() / tReserveSize );

            // free up unused memory
            aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mFollowerSideIpCellsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mLeaderSideIgCellSideOrdsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mFollowerSideIgCellSideOrdsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mNonTrivialFlagNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mTransitionLocationNew( iBspMesh ).shrink_to_fit();
        }

        // allocate ids for non-trivial integration cells
        moris_id tCurrentId = 0;

        // reserve global entity ids for additional elements that need to be created to get dbl. sided facets ...
        if ( aGhostSetupData.mLinearBackgroundMesh )
        {
            // ... for non-trivial element transitions only in the linear case
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tNonTrivialCount, mtk::EntityRank::ELEMENT );
        }
        else
        {
            // ... ?
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tEnrIntegMesh.get_num_entities( mtk::EntityRank::ELEMENT ), mtk::EntityRank::ELEMENT );
        }

        // Get next free index for new cells/elements
        moris_id tCurrentIndex = tEnrIntegMesh.get_num_entities( mtk::EntityRank::ELEMENT );

        // get total number of ghost facets for each B-spline mesh and bulk phase on current processor
        Matrix< DDUMat > tLocalNumberOfGhostFacets( tNumBspMeshes, tNumBulkPhases, 1 );
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            for ( uint iBulkPhase = 0; iBulkPhase < tNumBulkPhases; iBulkPhase++ )
            {
                tLocalNumberOfGhostFacets( iBspMesh, iBulkPhase ) = aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh )( iBulkPhase ).size();
            }
        }

        // get global (i.e. across all procs) number of ghost facets for each bulk phase
        Matrix< DDUMat > tTotalNumberOfGhostFacets = sum_all_matrix( tLocalNumberOfGhostFacets );

        // initialize list of double side set indices
        moris::Cell< moris_index > tDoubleSideSetIndexList;
        tDoubleSideSetIndexList.reserve( tNumBspMeshes * tNumBulkPhases );

        // initialize Ghost-Dbl.-SS. index
        moris_index tGhostDblSsIndex = 0;

        // ----------------------------------------------------------------------------------
        // loops constructing the actual ghost facets
        // iterate through B-spline meshes ...
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            // ... and bulk phases
            for ( uint iBulkPhase = 0; iBulkPhase < tNumBulkPhases; iBulkPhase++ )
            {
                // get the index of the current Double-Side-Set in the enriched integration mesh
                moris_index tCurrentDblSsIndexInEnrIgMesh = aGhostSetupData.mDblSideSetIndexInMesh( tGhostDblSsIndex );

                // store index in list of all Ghost Dbl.-SS.'s
                tDoubleSideSetIndexList.push_back( tCurrentDblSsIndexInEnrIgMesh );

                // get number of Ghost facets to be created
                uint tNumGhostFacetsInDblSS = aGhostSetupData.mLeaderSideIpCellsNew( iBspMesh )( iBulkPhase ).size();

                // resize list of Double sided clusters for current Dbl.-SS. in enr. IG mesh
                tEnrIntegMesh.mDoubleSideSets( tCurrentDblSsIndexInEnrIgMesh ).resize( tNumGhostFacetsInDblSS );

                // print number of Ghost facets to console
                MORIS_LOG_SPEC(
                        "Total Ghost Facets for B-spline mesh " + std::to_string( mMeshIndices( iBspMesh ) ) + " Bulk Phase " + std::to_string( iBulkPhase ),
                        tLocalNumberOfGhostFacets( iBspMesh, iBulkPhase ) );

                // iterate through double sides in this bulk phase
                for ( uint iGhostFacet = 0; iGhostFacet < tNumGhostFacetsInDblSS; iGhostFacet++ )
                {
                    // create a new single side cluster - follower side
                    std::shared_ptr< Side_Cluster > tFollowerSideCluster =
                            this->create_follower_side_cluster_new(
                                    aGhostSetupData,
                                    tEnrIpCells,
                                    iBspMesh,
                                    iBulkPhase,
                                    iGhostFacet,
                                    tCurrentIndex,
                                    tCurrentId );

                    // create a new single side cluster - leader side
                    std::shared_ptr< Side_Cluster > tLeaderSideCluster =
                            this->create_leader_side_cluster_new(
                                    aGhostSetupData,
                                    tEnrIpCells,
                                    iBspMesh,
                                    iBulkPhase,
                                    iGhostFacet,
                                    tFollowerSideCluster.get(),
                                    tCurrentIndex,
                                    tCurrentId );

                    // add single side cluster for leader side to the integration mesh
                    tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tLeaderSideCluster );

                    // store the index of the single side clusters associated with the double side clusters
                    tEnrIntegMesh.mDoubleSideSetsLeaderIndex( tCurrentDblSsIndexInEnrIgMesh ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() - 1 );

                    // add single side cluster for follower side to the integration mesh
                    tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tFollowerSideCluster );

                    // store the index of the single side clusters associated with the double side clusters
                    tEnrIntegMesh.mDoubleSideSetsFollowerIndex( tCurrentDblSsIndexInEnrIgMesh ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() - 1 );

                    // create double side cluster
                    std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster = std::make_shared< mtk::Double_Side_Cluster >(
                            tLeaderSideCluster.get(),
                            tFollowerSideCluster.get(),
                            tLeaderSideCluster->mVerticesInCluster );

                    // add double side cluster to list of all double side clusters in integration mesh
                    tEnrIntegMesh.mDoubleSideClusters.push_back( tDblSideCluster );

                    // add double side cluster to list of double side clusters on current set
                    tEnrIntegMesh.mDoubleSideSets( tCurrentDblSsIndexInEnrIgMesh )( iGhostFacet ) = tDblSideCluster;

                }    // end for: loop over Ghost facets in Dbl.-SS.

                // set color of the Double Side set for visualization
                tEnrIntegMesh.set_double_side_set_colors(
                        tCurrentDblSsIndexInEnrIgMesh,
                        { { tGhostDblSsIndex } },
                        { { tGhostDblSsIndex } } );

                // update index for the next Ghost Dbl.-SS.
                tGhostDblSsIndex++;

            }    // end for: loop over bulk phases
        }        // end for: loop over B-spline meshes

        // commit Ghost Dbl.-SS.'s to the enriched integration mesh
        tEnrIntegMesh.commit_double_side_set( tDoubleSideSetIndexList );
        tEnrIntegMesh.communicate_sets_of_type( mtk::SetType::DOUBLE_SIDED_SIDESET );

        // finalize enriched integration mesh
        tEnrIntegMesh.collect_all_sets();
    }

    // ----------------------------------------------------------------------------------

    bool
    Ghost_Stabilization::create_ghost(
            Ghost_Setup_Data&   aGhostSetupData,
            moris_index const & aFirstSubphase,
            moris_index const & aSecondSubphase,
            moris_index&        aTrivialFlag )
    {
        // Rules:
        // 1. Only create ghost facets between a subphase created inside an intersected
        //    cell and its neighbors.
        // 2. The owning processor of the leader (first) subphase constructs the ghost facet.
        // 3. Construct from coarse to fine in HMR

        // make sure flag is set to true, this is only turned to false on transition from coarse to fine cells
        aTrivialFlag = 0;

        moris_index tFirstSubphaseId  = mXTKModel->get_subphase_id( aFirstSubphase );
        moris_index tSecondSubphaseId = mXTKModel->get_subphase_id( aSecondSubphase );

        MORIS_ASSERT( tFirstSubphaseId != tSecondSubphaseId,
                "Subphase neighbor relation inconsistent\n" );

        // interpolation cell for this subphase
        moris_index tFirstInterpCell  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex( aFirstSubphase );
        moris_index tSecondInterpCell = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex( aSecondSubphase );

        // interpolation mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // MTK Cells
        mtk::Cell const & tFirstCell  = tEnrInterpMesh.get_mtk_cell( tFirstInterpCell );
        mtk::Cell const & tSecondCell = tEnrInterpMesh.get_mtk_cell( tSecondInterpCell );

        // Levels
        uint tFirstLevel  = tFirstCell.get_level();
        uint tSecondLevel = tSecondCell.get_level();

        // owners of interpolation cells
        moris_index tFirstOwnerIndex = tFirstCell.get_owner();
        // moris_index tSecondOwnerIndex = tSecondCell.get_owner();

        // proc rank
        moris_index tProcRank = par_rank();

        // if both subphases are not in child meshes return false
        moris::mtk::Cell* tSubphaseParentCell0 = mXTKModel->get_cut_integration_mesh()->get_subphase_parent_cell( aFirstSubphase );
        moris::mtk::Cell* tSubphaseParentCell1 = mXTKModel->get_cut_integration_mesh()->get_subphase_parent_cell( aSecondSubphase );
        if ( !mXTKModel->get_cut_integration_mesh()->parent_cell_has_children( tSubphaseParentCell0->get_index() )
                && !mXTKModel->get_cut_integration_mesh()->parent_cell_has_children( tSubphaseParentCell1->get_index() ) )
        {
            return false;
        }

        // Check based on refinement level of subphases

        // if the first subphase is finer than the second subphase,
        // do not construct ghost
        if ( tFirstLevel > tSecondLevel )
        {
            return false;
        }

        // if the first subphase is coarser than the second subphase
        // do construct ghost if proc owns first subphase
        if ( tFirstLevel < tSecondLevel )
        {
            if ( tFirstOwnerIndex == tProcRank )
            {
                aTrivialFlag = 1;
                return true;
            }

            return false;
        }

        // first set of checks passed - subphases are on same refinement level

        // do not construct if first subphase ID is smaller than second subphase ID
        if ( tFirstSubphaseId < tSecondSubphaseId )
        {
            return false;
        }

        // second set of checks passed - first subphase ID is larger than second one

        // do not construct if I do not own first subphase
        if ( tFirstOwnerIndex != tProcRank )
        {
            return false;
        }

        // if(tSecondOwnerIndex != tProcRank )
        // {
        //     return false;
        // }

        // third set of checks passed - first subphase ID is larger than second one and it is owned by current proc

        return true;
    }

    // ----------------------------------------------------------------------------------

    bool
    Ghost_Stabilization::create_ghost_new(
            moris_index const & aBspMeshListIndex,
            moris_index const & aFirstSpgIndex,
            moris_index const & aSecondSpgIndex,
            moris_index&        aNonTrivialFlag )
    {
        // make sure flag is set to true, this is only turned to false on transition from coarse to fine cells
        aNonTrivialFlag = 0;

        // quick access to the B-spline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mBsplineMeshInfos( aBspMeshListIndex );

        // get the B-spline element indices
        moris_index tFirstBspElemIndex  = tBsplineMeshInfo->mSubphaseGroups( aFirstSpgIndex )->get_bspline_cell_index();
        moris_index tSecondBspElemIndex = tBsplineMeshInfo->mSubphaseGroups( aSecondSpgIndex )->get_bspline_cell_index();

        // check that the two SPGs are on different elements
        MORIS_ASSERT( tFirstBspElemIndex != tSecondBspElemIndex,
                "Ghost_Stabilization::create_ghost_new() - Subphase group neighborhood relation inconsistent: SPGs on same B-spline element" );

        // check whether B-spline elements are intersected
        bool tFirstBspElemIsCut  = tBsplineMeshInfo->mSpgIndicesInBsplineCells( tFirstBspElemIndex ).size() > 1;
        bool tSecondBspElemIsCut = tBsplineMeshInfo->mSpgIndicesInBsplineCells( tSecondBspElemIndex ).size() > 1;

        // if neither is cut, don't construct a Ghost facet
        if ( !tFirstBspElemIsCut && !tSecondBspElemIsCut )
        {
            return false;
        }

        // get the owner of IP elements inside the first B-spline element
        mtk::Cell const * tFirstBaseIpCell = tBsplineMeshInfo->mExtractionCellsInBsplineCells( tFirstBspElemIndex )( 0 );
        moris_id          tFirstOwner      = tFirstBaseIpCell->get_owner();

        // get current proc rank
        moris_id tProcRank = par_rank();

        // if the IP Cells within the first B-spline element are not owned by current proc, don't construct a Ghost facet
        if ( tFirstOwner != tProcRank )
        {
            return false;
        }

        // get the refinement levels of the B-spline elements
        uint tFirstBspLevel  = tBsplineMeshInfo->mBsplineCellLevels( tFirstBspElemIndex );
        uint tSecondBspLevel = tBsplineMeshInfo->mBsplineCellLevels( tSecondBspElemIndex );

        // make sure ghost always gets constructed starting from the finer B-spline element
        if ( tFirstBspLevel < tSecondBspLevel )    // if current Bsp-element is coarser (i.e. on a lower refinement level) than the neighbor
        {
            return false;
        }

        // construct non-trivial ghost facet if current B-spline element is coarser
        if ( tFirstBspLevel > tSecondBspLevel )    // if current Bsp-element is coarser (on a lower refinement level) than the neighbor
        {
            aNonTrivialFlag = 1;
            return true;
        }

        // else: both B-spline elements are on the same refinement level

        // get the IDs of the Subphase groups
        moris_index tFirstSpgId  = tBsplineMeshInfo->mSubphaseGroups( aFirstSpgIndex )->get_id();
        moris_index tSecondSpgId = tBsplineMeshInfo->mSubphaseGroups( aSecondSpgIndex )->get_id();

        // check that the two SPG IDs are valid
        MORIS_ASSERT( ( tFirstSpgId != MORIS_ID_MAX ) && ( tSecondSpgId != MORIS_ID_MAX ),
                "Ghost_Stabilization::create_ghost_new() - At least one SPG ID within the pair is invalid." );

        // check that the two SPGs are different
        MORIS_ASSERT( tFirstSpgId != tSecondSpgId,
                "Ghost_Stabilization::create_ghost_new() - Subphase group neighborhood relation inconsistent: SPG IDs equal" );

        // ghost facet only gets constructed by element with SPG with bigger ID
        // to avoid duplicate Ghost-facets when B-spline elements are on same proc and same refinement level
        if ( tFirstSpgId > tSecondSpgId )
        {
            return true;
        }
        else
        {
            return false;
        }

        // check that all cases are covered by the logic; decision on ghost facet construction should exist at this point
        MORIS_ERROR( false, "Ghost_Stabilization::create_ghost_new() - logic broken; unknown if a ghost facet needs to be constructed or not." );
        return false;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_follower_side_cluster(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBulkIndex,
            uint const &                          aCellIndex,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // cut integration mesh access
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create a new side cluster for the follower
        std::shared_ptr< Side_Cluster > tFollowerSideCluster = std::make_shared< Side_Cluster >();

        // give the cluster the enriched interpolation cell
        tFollowerSideCluster->mInterpolationCell = aEnrIpCells( aGhostSetupData.mFollowerSideIpCells( aBulkIndex )( aCellIndex ) );

        // follower cluster is always trivial because the small facet is always the follower in the case of HMR hanging nodes
        tFollowerSideCluster->mTrivial = true;

        // create a linear IG cell over the whole Lagrange interpolation element and store this as being part of the new side cluster
        tFollowerSideCluster->mIntegrationCells = {
            this->get_linear_ig_cell(
                    aGhostSetupData,
                    aEnrIpCells( aGhostSetupData.mFollowerSideIpCells( aBulkIndex )( aCellIndex ) ),
                    aCurrentIndex,
                    aCurrentId )
        };

        // allocate space in integration cell side ordinals
        tFollowerSideCluster->mIntegrationCellSideOrdinals = { { aGhostSetupData.mFollowerSideIgCellSideOrds( aBulkIndex )( aCellIndex ) } };

        // add geometric vertices to the cluster // note: why is this
        // tFollowerSideCluster->mVerticesInCluster =
        //         tFollowerSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( tFollowerSideCluster->mIntegrationCellSideOrdinals( 0 ) );

        // add the corresponding cell (i.e. bulk) cluster to the new side cluster
        tFollowerSideCluster->mAssociatedCellCluster =
                &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tFollowerSideCluster->mInterpolationCell );

        // find and add the group of vertices to the new side cluster
        std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tFollowerSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );
        tFollowerSideCluster->set_ig_vertex_group( tVertexGroup );

        // return the shared pointer to the newly generated ghost follower side cluster
        return tFollowerSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_follower_side_cluster_new(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBsplineMeshListIndex,
            uint const &                          aBulkPhaseIndex,
            uint const &                          aGhostFacetIndexInSet,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // cut integration mesh access
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create a new side cluster for the follower
        std::shared_ptr< Side_Cluster > tFollowerSideCluster = std::make_shared< Side_Cluster >();

        // get the index to the Enr. IP cell the side cluster is attached to
        moris_index tFollowerEnrIpCellIndex = aGhostSetupData.mFollowerSideIpCellsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );

        // give the cluster the enriched interpolation cell
        tFollowerSideCluster->mInterpolationCell = aEnrIpCells( tFollowerEnrIpCellIndex );

        // follower cluster is always trivial because the small facet is always the follower in the case of HMR hanging nodes
        tFollowerSideCluster->mTrivial = true;

        // add integration cell
        tFollowerSideCluster->mIntegrationCells = {
            this->get_linear_ig_cell(
                    aGhostSetupData,
                    aEnrIpCells( tFollowerEnrIpCellIndex ),
                    aCurrentIndex,
                    aCurrentId )
        };

        // allocate space in integration cell side ordinals
        moris_index tSideOrdinal                           = aGhostSetupData.mFollowerSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );
        tFollowerSideCluster->mIntegrationCellSideOrdinals = { { tSideOrdinal } };

        // assign bulk cluster corresponding to new side cluster
        tFollowerSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tFollowerSideCluster->mInterpolationCell );

        // leader vertex group
        std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tFollowerSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

        tFollowerSideCluster->set_ig_vertex_group( tVertexGroup );

        return tFollowerSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_leader_side_cluster(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBulkIndex,
            uint const &                          aCellIndex,
            Side_Cluster*                         aFollowerSideCluster,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // get access to the cut integration mesh
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create a leader side cluster
        std::shared_ptr< Side_Cluster > tLeaderSideCluster = std::make_shared< Side_Cluster >();

        // assign the UIPC to the new side cluster
        tLeaderSideCluster->mInterpolationCell = aEnrIpCells( aGhostSetupData.mLeaderSideIpCells( aBulkIndex )( aCellIndex ) );

        // Case: non-trivial Lagrange element transition
        if ( aGhostSetupData.mNonTrivialFlag( aBulkIndex )( aCellIndex ) > 0 )
        {
            // flag the leader side as non-trivial
            tLeaderSideCluster->mTrivial = false;

            // create new integration cell using the vertices on the follower facet and the adjacent vertices of the base interpolation cell
            std::shared_ptr< xtk::Cell_XTK_No_CM > tNewIgCell = this->create_non_trivial_leader_ig_cell(
                    aGhostSetupData,
                    aBulkIndex,
                    aCellIndex,
                    tLeaderSideCluster.get(),
                    aFollowerSideCluster,
                    aCurrentIndex,
                    aCurrentId );

            // get the local coordinates of the IP cell vertices including the hanging vertex
            Cell< Matrix< DDRMat > > tLocCoords;
            this->get_local_coords_on_transition_side(
                    aGhostSetupData.mLeaderSideIgCellSideOrds( aBulkIndex )( aCellIndex ),
                    aGhostSetupData.mTransitionLocation( aBulkIndex )( aCellIndex ),
                    tLocCoords );

            // add integration cell
            tLeaderSideCluster->mIntegrationCells.push_back( tNewIgCell.get() );

            // add side ordinal relative to the integration cell
            tLeaderSideCluster->mIntegrationCellSideOrdinals = { { this->get_side_ordinals_for_non_trivial_leader() } };

            // get the vertex group
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tLeaderSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

            // Add the follower vertex to the vertex group
            moris::Cell< moris::mtk::Vertex const * > tFollowerVertices =
                    aFollowerSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aFollowerSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // iterate through vertices
            for ( uint iV = 0; iV < tFollowerVertices.size(); iV++ )
            {
                if ( !tVertexGroup->vertex_is_in_group( tFollowerVertices( iV )->get_index() ) )
                {
                    tVertexGroup->add_vertex( tFollowerVertices( iV ), std::make_shared< Matrix< DDRMat > >( tLocCoords( iV ) ) );
                }

                MORIS_ASSERT(
                        moris::norm( ( *tVertexGroup->get_vertex_local_coords( tFollowerVertices( iV )->get_index() ) - tLocCoords( iV ) ) ) < 1e-8,
                        "Ghost_Stabilization::create_leader_side_cluster() - Local coord issue" );
            }

            // set the vertex group
            tLeaderSideCluster->set_ig_vertex_group( tVertexGroup );

            // // get the vertices on the side ordinal
            // tLeaderSideCluster->mVerticesInCluster = tFollowerVertices;

            // add the local coordinates
            tLeaderSideCluster->mVertexLocalCoords = tLocCoords;

            // associated cell cluster
            tLeaderSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tLeaderSideCluster->mInterpolationCell );

            // verify new cluster
            mtk::Mesh_Checker tCheck;
            MORIS_ASSERT(
                    tCheck.verify_side_cluster( tLeaderSideCluster.get(), mtk::Leader_Follower::LEADER ),
                    "Ghost_Stabilization::create_leader_side_cluster() - "
                    "Invalid Side Cluster Created Check the local coordinates" );

            // place the new ig cell in the background mesh
            mXTKModel->get_cut_integration_mesh()->add_integration_cell( tNewIgCell->get_index(), tNewIgCell );
        }

        // Case: trivial lagrange element transition
        else
        {
            // flag the leader side as trivial
            tLeaderSideCluster->mTrivial = true;

            // create a linear IG cell over the whole Lagrange interpolation element and store this as being part of the new side cluster
            tLeaderSideCluster->mIntegrationCells = {
                this->get_linear_ig_cell(
                        aGhostSetupData,
                        aEnrIpCells( aGhostSetupData.mLeaderSideIpCells( aBulkIndex )( aCellIndex ) ),
                        aCurrentIndex,
                        aCurrentId )
            };

            // add side ordinal on which the ghost facet sits relative to the leader cluster
            moris_index tSideOrdinal                         = aGhostSetupData.mLeaderSideIgCellSideOrds( aBulkIndex )( aCellIndex );
            tLeaderSideCluster->mIntegrationCellSideOrdinals = { { tSideOrdinal } };

            // get and add the vertices on that side ordinal
            tLeaderSideCluster->mVerticesInCluster =
                    tLeaderSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( tLeaderSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // get and set the associated Cell (i.e. bulk) cluster
            tLeaderSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tLeaderSideCluster->mInterpolationCell );

            // get and set the vertex group from the IP cell
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tLeaderSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );
            tLeaderSideCluster->set_ig_vertex_group( tVertexGroup );
        }

        // return the shared pointer to the leader side cluster constructed
        return tLeaderSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_leader_side_cluster_new(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBsplineMeshListIndex,
            uint const &                          aBulkPhaseIndex,
            uint const &                          aGhostFacetIndexInSet,
            Side_Cluster*                         aFollowerSideCluster,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // get access to the cut integration mesh
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create the leader side cluster
        std::shared_ptr< Side_Cluster > tLeaderSideCluster = std::make_shared< Side_Cluster >();

        // get the index to the Enr. IP cell the side cluster is attached to
        moris_index tLeaderEnrIpCellIndex = aGhostSetupData.mLeaderSideIpCellsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );

        // set the enr. IP cell the side cluster is attached to
        tLeaderSideCluster->mInterpolationCell = aEnrIpCells( tLeaderEnrIpCellIndex );

        // Setup the leader side cluster
        if ( aGhostSetupData.mNonTrivialFlagNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet ) > 0 )
        {
            // flag the leader side as non-trivial
            tLeaderSideCluster->mTrivial = false;

            // create new integration cell using the vertices on the follower facet and the adjacent vertices of the base interpolation cell
            std::shared_ptr< xtk::Cell_XTK_No_CM > tNewIgCell = this->create_non_trivial_leader_ig_cell_new(
                    aGhostSetupData,
                    aBsplineMeshListIndex,
                    aBulkPhaseIndex,
                    aGhostFacetIndexInSet,
                    tLeaderSideCluster.get(),
                    aFollowerSideCluster,
                    aCurrentIndex,
                    aCurrentId );

            // get the local coordinates of the IP cell vertices including the hanging vertex
            Cell< Matrix< DDRMat > > tLocCoords;
            this->get_local_coords_on_transition_side(
                    aGhostSetupData.mLeaderSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet ),
                    aGhostSetupData.mTransitionLocationNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet ),
                    tLocCoords );

            // add integration cell
            tLeaderSideCluster->mIntegrationCells.push_back( tNewIgCell.get() );

            // add side ordinal relative to the integration cell
            tLeaderSideCluster->mIntegrationCellSideOrdinals = { { this->get_side_ordinals_for_non_trivial_leader() } };

            // get the vertex group
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tLeaderSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

            // Add the follower vertex to the vertex group
            moris::Cell< moris::mtk::Vertex const * > tFollowerVertices =
                    aFollowerSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aFollowerSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // iterate through vertices
            for ( uint iV = 0; iV < tFollowerVertices.size(); iV++ )
            {
                if ( !tVertexGroup->vertex_is_in_group( tFollowerVertices( iV )->get_index() ) )
                {
                    tVertexGroup->add_vertex( tFollowerVertices( iV ), std::make_shared< Matrix< DDRMat > >( tLocCoords( iV ) ) );
                }

                MORIS_ASSERT(
                        moris::norm( ( *tVertexGroup->get_vertex_local_coords( tFollowerVertices( iV )->get_index() ) - tLocCoords( iV ) ) ) < 1e-8,
                        "Ghost_Stabilization::create_leader_side_cluster() - Local coord issue" );
            }

            // set the vertex group
            tLeaderSideCluster->set_ig_vertex_group( tVertexGroup );

            // add the local coordinates
            tLeaderSideCluster->mVertexLocalCoords = tLocCoords;

            // associated cell cluster
            tLeaderSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tLeaderSideCluster->mInterpolationCell );

            // verify new cluster
            mtk::Mesh_Checker tCheck;
            MORIS_ASSERT(
                    tCheck.verify_side_cluster( tLeaderSideCluster.get(), mtk::Leader_Follower::LEADER ),
                    "Ghost_Stabilization::create_leader_side_cluster_new() - "
                    "Invalid Side Cluster Created Check the local coordinates" );

            // place the new ig cell in the background mesh
            mXTKModel->get_cut_integration_mesh()->add_integration_cell( tNewIgCell->get_index(), tNewIgCell );
        }

        // Case: trivial Lagrange element transition
        else
        {
            // flag the leader side as trivial
            tLeaderSideCluster->mTrivial = true;

            // add integration cell
            tLeaderSideCluster->mIntegrationCells = {
                this->get_linear_ig_cell(
                        aGhostSetupData,
                        aEnrIpCells( tLeaderEnrIpCellIndex ),
                        aCurrentIndex,
                        aCurrentId )
            };

            // add side ordinal relative to the integration cell
            moris_index tSideOrdinal                         = aGhostSetupData.mLeaderSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );
            tLeaderSideCluster->mIntegrationCellSideOrdinals = { { tSideOrdinal } };

            // add the vertices on the side ordinal
            tLeaderSideCluster->mVerticesInCluster =
                    tLeaderSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( tLeaderSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // get and set the cell cluster
            tLeaderSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tLeaderSideCluster->mInterpolationCell );

            // get the vertex group from the IP cell
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tLeaderSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

            // set Side cluster to use vertex group from IP cell
            tLeaderSideCluster->set_ig_vertex_group( tVertexGroup );
        }

        // return side cluster constructed
        return tLeaderSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< xtk::Cell_XTK_No_CM >
    Ghost_Stabilization::create_non_trivial_leader_ig_cell(
            Ghost_Setup_Data& aGhostSetupData,
            uint const &      aBulkIndex,
            uint const &      aCellIndex,
            Side_Cluster*     aLeaderSideCluster,
            Side_Cluster*     aFollowerSideCluster,
            moris_index&      aCurrentIndex,
            moris_index&      aCurrentId )
    {
        // make sure the corresponding follower side cluster is valid
        MORIS_ASSERT( aFollowerSideCluster->get_cells_in_side_cluster().size() == 1,
                "Ghost_Stabilization::create_non_trivial_leader_ig_cell() - "
                "Follower side cluster should have exactly one integration cell." );

        // get the vertices on the side for the follower side cluster
        moris::Cell< moris::mtk::Vertex const * > tFollowerVertices =
                aFollowerSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aFollowerSideCluster->mIntegrationCellSideOrdinals( 0 ) );

        // get the leader UIPC
        Interpolation_Cell_Unzipped const * tLeaderIpCell = aLeaderSideCluster->mInterpolationCell;

        // base leader IP cell
        moris::mtk::Cell const * tBaseLeaderCell = tLeaderIpCell->get_base_cell();

        // get the connectivity information from the cell
        std::shared_ptr< moris::mtk::Cell_Info > tCellInfo = tLeaderIpCell->get_cell_info_sp();

        // adjacent side ordinal on leader
        uint tAdjFacetOrd = tCellInfo->get_adjacent_side_ordinal(
                aGhostSetupData.mLeaderSideIgCellSideOrds( aBulkIndex )( aCellIndex ) );

        // setup the vertices and local coordinates of the vertices relative to the cell
        moris::Cell< moris::mtk::Vertex const * > tAdjVertices =
                tBaseLeaderCell->get_geometric_vertices_on_side_ordinal( tAdjFacetOrd );

        // line up the ordering of the leader and follower side vertices such that they match
        moris::Cell< moris::mtk::Vertex const * > tPermutedFollowerVertices;
        moris::Cell< moris::mtk::Vertex const * > tPermutedAdjVertices;
        this->permute_follower_vertices(
                tFollowerVertices,
                tAdjVertices,
                tPermutedFollowerVertices,
                tPermutedAdjVertices );

        // collect all vertices on the ghost facet
        moris::Cell< moris::mtk::Vertex* > tCellVertices( tPermutedAdjVertices.size() + tPermutedFollowerVertices.size() );

        uint tCount = 0;
        for ( auto iV : tPermutedAdjVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }
        for ( auto iV : tPermutedFollowerVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }

        // linear cell info
        mtk::Cell_Info_Factory                   tCellInfoFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo =
                tCellInfoFactory.create_cell_info_sp( tLeaderIpCell->get_geometry_type(), mtk::Interpolation_Order::LINEAR );

        // create a new integration cell that does not have a child mesh association
        std::shared_ptr< xtk::Cell_XTK_No_CM > tIgCell = std::make_shared< xtk::Cell_XTK_No_CM >(
                aCurrentId,
                aCurrentIndex,
                tLeaderIpCell->get_owner(),
                tLinearCellInfo,
                tCellVertices );

        // increment current id and index
        aCurrentId++;
        aCurrentIndex++;

        return tIgCell;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< xtk::Cell_XTK_No_CM >
    Ghost_Stabilization::create_non_trivial_leader_ig_cell_new(
            Ghost_Setup_Data& aGhostSetupData,
            uint const &      aBsplineMeshListIndex,
            uint const &      aBulkIndex,
            uint const &      aCellIndex,
            Side_Cluster*     aLeaderSideCluster,
            Side_Cluster*     aFollowerSideCluster,
            moris_index&      aCurrentIndex,
            moris_index&      aCurrentId )
    {
        // make sure the corresponding follower side cluster is valid
        MORIS_ASSERT( aFollowerSideCluster->get_cells_in_side_cluster().size() == 1,
                "Ghost_Stabilization::create_non_trivial_leader_ig_cell() - "
                "Follower side cluster should have exactly one integration cell." );

        // get the vertices on the side for the follower side cluster
        moris::Cell< moris::mtk::Vertex const * > tFollowerVertices =
                aFollowerSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aFollowerSideCluster->mIntegrationCellSideOrdinals( 0 ) );

        // get the leader UIPC
        Interpolation_Cell_Unzipped const * tLeaderIpCell = aLeaderSideCluster->mInterpolationCell;

        // base leader IP cell
        moris::mtk::Cell const * tBaseLeaderCell = tLeaderIpCell->get_base_cell();

        // get the connectivity information from the cell
        std::shared_ptr< moris::mtk::Cell_Info > tCellInfo = tLeaderIpCell->get_cell_info_sp();

        // adjacent side ordinal on leader
        uint tAdjFacetOrd = tCellInfo->get_adjacent_side_ordinal(
                aGhostSetupData.mLeaderSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkIndex )( aCellIndex ) );

        // setup the vertices and local coordinates of the vertices relative to the cell
        moris::Cell< moris::mtk::Vertex const * > tAdjVertices =
                tBaseLeaderCell->get_geometric_vertices_on_side_ordinal( tAdjFacetOrd );

        // line up the ordering of the leader and follower side vertices such that they match
        moris::Cell< moris::mtk::Vertex const * > tPermutedFollowerVertices;
        moris::Cell< moris::mtk::Vertex const * > tPermutedAdjVertices;
        this->permute_follower_vertices(
                tFollowerVertices,
                tAdjVertices,
                tPermutedFollowerVertices,
                tPermutedAdjVertices );

        // collect all vertices on the ghost facet
        moris::Cell< moris::mtk::Vertex* > tCellVertices( tPermutedAdjVertices.size() + tPermutedFollowerVertices.size() );

        uint tCount = 0;
        for ( auto iV : tPermutedAdjVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }
        for ( auto iV : tPermutedFollowerVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }

        // linear cell info
        mtk::Cell_Info_Factory                   tCellInfoFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo =
                tCellInfoFactory.create_cell_info_sp( tLeaderIpCell->get_geometry_type(), mtk::Interpolation_Order::LINEAR );

        // create a new integration cell that does not have a child mesh association
        std::shared_ptr< xtk::Cell_XTK_No_CM > tIgCell = std::make_shared< xtk::Cell_XTK_No_CM >(
                aCurrentId,
                aCurrentIndex,
                tLeaderIpCell->get_owner(),
                tLinearCellInfo,
                tCellVertices );

        // increment current id and index
        aCurrentId++;
        aCurrentIndex++;

        return tIgCell;
    }

    // ----------------------------------------------------------------------------------

    mtk::Cell*
    Ghost_Stabilization::get_linear_ig_cell(
            Ghost_Setup_Data&            aGhostSetupData,
            Interpolation_Cell_Unzipped* aInterpCell,
            moris_index&                 aCurrentIndex,
            moris_index&                 aCurrentId )
    {
        mtk::Cell* tCell = nullptr;

        auto tIter = aGhostSetupData.mLinearIgCellIndex.find( aInterpCell->get_id() );
        if ( aGhostSetupData.mLinearBackgroundMesh == true )
        {
            tCell = aInterpCell->get_base_cell();
        }
        else if ( tIter == aGhostSetupData.mLinearIgCellIndex.end() )
        {
            tCell = this->create_linear_ig_cell( aGhostSetupData, aInterpCell, aCurrentIndex, aCurrentId );
        }
        else
        {
            moris_index tIndex = tIter->second;
            tCell              = aGhostSetupData.mLinearIgCells( tIndex );
        }

        return tCell;
    }
    // ----------------------------------------------------------------------------------
    mtk::Cell*
    Ghost_Stabilization::create_linear_ig_cell( Ghost_Setup_Data& aGhostSetupData,
            Interpolation_Cell_Unzipped const *                   aInterpCell,
            moris_index&                                          aCurrentIndex,
            moris_index&                                          aCurrentId )
    {
        MORIS_ASSERT( aGhostSetupData.mLinearIgCellIndex.find( aInterpCell->get_id() ) == aGhostSetupData.mLinearIgCellIndex.end(), "Trying to create linear ig cell twice" );

        // get the base cell
        moris::mtk::Cell const * tBaseCell = aInterpCell->get_base_cell();

        // get the geometric vertices
        moris::Cell< moris::mtk::Vertex* > tAllVertices = tBaseCell->get_vertex_pointers();

        moris::Cell< moris::mtk::Vertex* > tLinearVertices;

        uint tSpatialDim = mXTKModel->get_spatial_dim();

        MORIS_ERROR( tBaseCell->get_geometry_type() == mtk::Geometry_Type::HEX || tBaseCell->get_geometry_type() == mtk::Geometry_Type::QUAD,
                "Ghost only tested on quads and hex" );

        if ( tSpatialDim == 2 )
        {
            tLinearVertices.resize( 4 );
            tLinearVertices( 0 ) = tAllVertices( 0 );
            tLinearVertices( 1 ) = tAllVertices( 1 );
            tLinearVertices( 2 ) = tAllVertices( 2 );
            tLinearVertices( 3 ) = tAllVertices( 3 );
        }
        else if ( tSpatialDim == 3 )
        {
            tLinearVertices.resize( 8 );

            tLinearVertices( 0 ) = tAllVertices( 0 );
            tLinearVertices( 1 ) = tAllVertices( 1 );
            tLinearVertices( 2 ) = tAllVertices( 2 );
            tLinearVertices( 3 ) = tAllVertices( 3 );
            tLinearVertices( 4 ) = tAllVertices( 4 );
            tLinearVertices( 5 ) = tAllVertices( 5 );
            tLinearVertices( 6 ) = tAllVertices( 6 );
            tLinearVertices( 7 ) = tAllVertices( 7 );
        }
        else
        {
            MORIS_ERROR( 0, "Only 2/3d" );
        }

        // linear cell info
        mtk::Cell_Info_Factory                   tCellInfoFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo = tCellInfoFactory.create_cell_info_sp( aInterpCell->get_geometry_type(), mtk::Interpolation_Order::LINEAR );


        // create a new integration cell that does not have a child mesh association
        std::shared_ptr< xtk::Cell_XTK_No_CM > tIgCell = std::make_shared< xtk::Cell_XTK_No_CM >(
                aCurrentId,
                aCurrentIndex,
                aInterpCell->get_owner(),
                tLinearCellInfo,
                tLinearVertices );


        // add to map
        aGhostSetupData.mLinearIgCellIndex[ aInterpCell->get_id() ] = (moris_index)aGhostSetupData.mLinearIgCells.size();

        // add to data
        aGhostSetupData.mLinearIgCells.push_back( tIgCell.get() );

        // add to background mesh
        mXTKModel->get_cut_integration_mesh()->add_integration_cell( tIgCell->get_index(), tIgCell );

        aCurrentIndex++;
        aCurrentId++;

        return tIgCell.get();
    }


    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::permute_follower_vertices(
            moris::Cell< moris::mtk::Vertex const * > const & aFollowerVertices,
            moris::Cell< moris::mtk::Vertex const * > const & aLeaderVertices,
            moris::Cell< moris::mtk::Vertex const * >&        aPermutedFollowerVertices,
            moris::Cell< moris::mtk::Vertex const * >&        aPermutedAdjMastVertices )
    {
        uint tSpatialDim = mXTKModel->get_spatial_dim();

        switch ( tSpatialDim )
        {
            case 2:
            {
                aPermutedFollowerVertices = { aFollowerVertices( 1 ), aFollowerVertices( 0 ) };
                aPermutedAdjMastVertices  = { aLeaderVertices( 0 ), aLeaderVertices( 1 ) };
                break;
            }
            default:
            {
                aPermutedFollowerVertices = { aFollowerVertices( 0 ), aFollowerVertices( 3 ), aFollowerVertices( 2 ), aFollowerVertices( 1 ) };
                aPermutedAdjMastVertices  = { aLeaderVertices( 0 ), aLeaderVertices( 3 ), aLeaderVertices( 2 ), aLeaderVertices( 1 ) };
            }
        }
    }
    // ----------------------------------------------------------------------------------
    void
    Ghost_Stabilization::get_local_coords_on_transition_side(
            moris_index const &       aMySideOrdinal,
            moris_index const &       aTransitionLoc,
            Cell< Matrix< DDRMat > >& aLocCoord )
    {
        uint tTag = 100 * mXTKModel->get_spatial_dim() + 10 * aMySideOrdinal + aTransitionLoc;

        switch ( tTag )
        {
            case 200:
                aLocCoord = { { { 0, -1 } }, { { -1, -1 } } };
                break;
            case 201:
                aLocCoord = { { { 1, -1 } }, { { 0, -1 } } };
                break;
            case 211:
                aLocCoord = { { { 1, 0 } }, { { 1, -1 } } };
                break;
            case 213:
                aLocCoord = { { { 1, 1 } }, { { 1, 0 } } };
                break;
            case 223:
                aLocCoord = { { { -1, 1 } }, { { 0, 1 } } };
                break;
            case 222:
                aLocCoord = { { { 0, 1 } }, { { 1, 1 } } };
                break;
            case 230:
                aLocCoord = { { { -1, 0 } }, { { -1, 1 } } };
                break;
            case 232:
                aLocCoord = { { { -1, -1 } }, { { -1, 0 } } };
                break;
            case 300:
                aLocCoord = { { { -1, -1, -1 } }, { { -1, -1, 0 } }, { { 0, -1, 0 } }, { { 0, -1, -1 } } };
                break;
            case 301:
                aLocCoord = { { { 0, -1, -1 } }, { { 0, -1, 0 } }, { { 1, -1, 0 } }, { { 1, -1, -1 } } };
                break;
            case 304:
                aLocCoord = { { { 0, -1, 0 } }, { { 0, -1, 1 } }, { { 1, -1, 1 } }, { { 1, -1, 0 } } };
                break;
            case 305:
                aLocCoord = { { { -1, -1, 0 } }, { { -1, -1, 1 } }, { { 0, -1, 1 } }, { { 0, -1, 0 } } };
                break;
            case 311:
                aLocCoord = { { { 1, -1, -1 } }, { { 1, -1, 0 } }, { { 1, 0, 0 } }, { { 1, 0, -1 } } };
                break;
            case 313:
                aLocCoord = { { { 1, 0, -1 } }, { { 1, 0, 0 } }, { { 1, 1, 0 } }, { { 1, 1, -1 } } };
                break;
            case 315:
                aLocCoord = { { { 1, 0, 0 } }, { { 1, 0, 1 } }, { { 1, 1, 1 } }, { { 1, 1, 0 } } };
                break;
            case 317:
                aLocCoord = { { { 1, -1, 0 } }, { { 1, -1, 1 } }, { { 1, 0, 1 } }, { { 1, 0, 0 } } };
                break;
            case 322:
                aLocCoord = { { { 1, 1, -1 } }, { { 1, 1, 0 } }, { { 0, 1, 0 } }, { { 0, 1, -1 } } };
                break;
            case 323:
                aLocCoord = { { { 0, 1, -1 } }, { { 0, 1, 0 } }, { { -1, 1, 0 } }, { { -1, 1, -1 } } };
                break;
            case 326:
                aLocCoord = { { { 0, 1, 0 } }, { { 0, 1, 1 } }, { { -1, 1, 1 } }, { { -1, 1, 0 } } };
                break;
            case 327:
                aLocCoord = { { { 1, 1, 0 } }, { { 1, 1, 1 } }, { { 0, 1, 1 } }, { { 0, 1, 0 } } };
                break;
            case 330:
                aLocCoord = { { { -1, 0, -1 } }, { { -1, 1, -1 } }, { { -1, 1, 0 } }, { { -1, 0, 0 } } };
                break;
            case 332:
                aLocCoord = { { { -1, -1, -1 } }, { { -1, 0, -1 } }, { { -1, 0, 0 } }, { { -1, -1, 0 } } };
                break;
            case 334:
                aLocCoord = { { { -1, -1, 0 } }, { { -1, 0, 0 } }, { { -1, 0, 1 } }, { { -1, -1, 1 } } };
                break;
            case 336:
                aLocCoord = { { { -1, 0, 0 } }, { { -1, 1, 0 } }, { { -1, 1, 1 } }, { { -1, 0, 1 } } };
                break;
            case 340:
                aLocCoord = { { { 0, 0, -1 } }, { { 1, 0, -1 } }, { { 1, 1, -1 } }, { { 0, 1, -1 } } };
                break;
            case 341:
                aLocCoord = { { { -1, 0, -1 } }, { { 0, 0, -1 } }, { { 0, 1, -1 } }, { { -1, 1, -1 } } };
                break;
            case 342:
                aLocCoord = { { { -1, -1, -1 } }, { { 0, -1, -1 } }, { { 0, 0, -1 } }, { { -1, 0, -1 } } };
                break;
            case 343:
                aLocCoord = { { { 0, -1, -1 } }, { { 1, -1, -1 } }, { { 1, 0, -1 } }, { { 0, 0, -1 } } };
                break;
            case 354:
                aLocCoord = { { { -1, -1, 1 } }, { { -1, 0, 1 } }, { { 0, 0, 1 } }, { { 0, -1, 1 } } };
                break;
            case 355:
                aLocCoord = { { { 0, -1, 1 } }, { { 0, 0, 1 } }, { { 1, 0, 1 } }, { { 1, -1, 1 } } };
                break;
            case 356:
                aLocCoord = { { { 0, 0, 1 } }, { { 0, 1, 1 } }, { { 1, 1, 1 } }, { { 1, 0, 1 } } };
                break;
            case 357:
                aLocCoord = { { { -1, 0, 1 } }, { { -1, 1, 1 } }, { { 0, 1, 1 } }, { { 0, 0, 1 } } };
                break;
            default:
                MORIS_ERROR( 0, "Invalid tag (100*spatial dim + 10 * side ord + transition location)" );
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Ghost_Stabilization::get_side_ordinals_for_non_trivial_leader()
    {
        if ( mXTKModel->get_spatial_dim() == 2 )
        {
            return 2;
        }
        else if ( mXTKModel->get_spatial_dim() == 3 )
        {
            return 5;
        }
        else
        {
            MORIS_ERROR( false,
                    "Ghost_Stabilization::get_side_ordinals_for_non_trivial_leader() - "
                    "Invalid spatial dimension for ghost" );

            return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Ghost_Stabilization::is_linear_ip_mesh()
    {
        Enriched_Interpolation_Mesh& tEnrIpMesh = mXTKModel->get_enriched_interp_mesh( 0 );

        MORIS_ERROR( tEnrIpMesh.get_num_entities( mtk::EntityRank::ELEMENT ) > 0, "Cannot deduce type on empty mesh." );

        // get the first cell
        mtk::Cell const & tCell0 = tEnrIpMesh.get_mtk_cell( 0 );

        if ( tCell0.get_interpolation_order() == mtk::Interpolation_Order::LINEAR )
        {
            return true;
        }

        else
        {
            return false;
        }
    }

    // ----------------------------------------------------------------------------------

}    // namespace xtk
