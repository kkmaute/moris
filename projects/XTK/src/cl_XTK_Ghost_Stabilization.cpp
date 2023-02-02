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

        // initialize object carrying generateddata for ghost facets
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
        enum CellTopology tFacetTopo = determine_cell_topology( mXTKModel->get_spatial_dim(), mtk::Interpolation_Order::LINEAR, CellShape::RECTANGULAR );

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
        enum CellTopology tFacetTopo = determine_cell_topology( mXTKModel->get_spatial_dim(), mtk::Interpolation_Order::LINEAR, CellShape::RECTANGULAR );

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
        uint tCurrentNewInterpCellIndex = tEnrIpMesh.get_num_entities( EntityRank::ELEMENT );

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
        moris_id tCurrentId = tEnrIpMesh.allocate_entity_ids( tNumNewInterpCellsOwned, EntityRank::ELEMENT, false );

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
        Tracer tTracer( "XTK", "Ghost", "identify_and_setup_aura_vertices_in_ghost" );
        // access enriched ip mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // vertices in ghost that have interpolation
        moris::Cell< mtk::Vertex* > tGhostVerticesWithInterpolation;

        // vertices in ghost that do not have interpolation
        moris::Cell< mtk::Vertex* > tGhostVerticesWithoutInterpolation;

        // an interpolation cell that the ghost vertex without interpolation is connected to
        // this is needed for communication routine
        moris::Cell< mtk::Cell const * > tGhostIpCellConnectedToVertex;

        this->get_ip_vertices_in_ghost_sets(
                aGhostSetupData,
                tGhostVerticesWithInterpolation,
                tGhostVerticesWithoutInterpolation,
                tGhostIpCellConnectedToVertex );

        for ( uint iMeshIndex = 0; iMeshIndex < tEnrInterpMesh.mMeshIndices.numel(); iMeshIndex++ )
        {
            // current mesh index
            moris_index tMeshIndex = tEnrInterpMesh.mMeshIndices( iMeshIndex );

            // sort the ghost vertices without interpolation by proc
            Cell< Matrix< IndexMat > >               tNotOwnedIPVertIndsToProcs;
            Cell< Matrix< IndexMat > >               tNotOwnedBGIPVertsIdsToProcs;
            Cell< Matrix< IndexMat > >               tNotOwnedEnrichedCellIdToProcs;
            Cell< Matrix< IndexMat > >               tNotOwnedEnrichedCellBulkPhaseToProcs;    // for checking against
            Cell< uint >                             tProcRanks;
            std::unordered_map< moris_id, moris_id > tProcRankToDataIndex;
            this->prepare_interpolation_vertex_t_matrix_requests(
                    tGhostVerticesWithoutInterpolation,
                    tGhostIpCellConnectedToVertex,
                    tNotOwnedIPVertIndsToProcs,
                    tNotOwnedBGIPVertsIdsToProcs,
                    tNotOwnedEnrichedCellIdToProcs,
                    tNotOwnedEnrichedCellBulkPhaseToProcs,
                    tProcRanks,
                    tProcRankToDataIndex );

            // send requests
            uint tMPITag = 3001;

            // send the background vertex id
            mXTKModel->send_outward_requests( tMPITag, tProcRanks, tNotOwnedBGIPVertsIdsToProcs );

            // send the enriched interpolation cell id
            mXTKModel->send_outward_requests( tMPITag + 1, tProcRanks, tNotOwnedEnrichedCellIdToProcs );

            // send the enriched interpolation cell bulk phase ids
            mXTKModel->send_outward_requests( tMPITag + 2, tProcRanks, tNotOwnedEnrichedCellBulkPhaseToProcs );

            barrier();

            // receive requests for the t-matrices
            Cell< Matrix< IndexMat > > tReceivedVertexIds;
            Cell< Matrix< IndexMat > > tReceivedEnrichedCellId;
            Cell< Matrix< IndexMat > > tReceivedEnrichedCellBulkPhase;
            Cell< uint >               tProcsReceivedFrom1;
            Cell< uint >               tProcsReceivedFrom2;
            Cell< uint >               tProcsReceivedFrom3;
            mXTKModel->inward_receive_requests( tMPITag, 1, tReceivedVertexIds, tProcsReceivedFrom1 );                    // receive the requests ofr BG VertexIds
            mXTKModel->inward_receive_requests( tMPITag + 1, 1, tReceivedEnrichedCellId, tProcsReceivedFrom2 );           // recieve the requests for Enriched IP Cell Id
            mXTKModel->inward_receive_requests( tMPITag + 2, 1, tReceivedEnrichedCellBulkPhase, tProcsReceivedFrom3 );    // recieve the requests for Enriched IP Cell Bulk phases

            // prepare the t-matrices for sending
            Cell< Matrix< DDRMat > >   tTMatrixWeights;
            Cell< Matrix< IndexMat > > tTMatrixIndices;
            Cell< Matrix< IndexMat > > tTMatrixOwners;
            Cell< Matrix< IndexMat > > tTMatrixOffsets;
            this->prepare_t_matrix_request_answers(
                    tMeshIndex,
                    tReceivedVertexIds,
                    tReceivedEnrichedCellId,
                    tReceivedEnrichedCellBulkPhase,
                    tTMatrixWeights,
                    tTMatrixIndices,
                    tTMatrixOwners,
                    tTMatrixOffsets );

            // send information
            mXTKModel->return_request_answers_reals( tMPITag + 3, tTMatrixWeights, tProcsReceivedFrom1 );
            mXTKModel->return_request_answers( tMPITag + 4, tTMatrixIndices, tProcsReceivedFrom1 );
            mXTKModel->return_request_answers( tMPITag + 5, tTMatrixOwners, tProcsReceivedFrom1 );
            mXTKModel->return_request_answers( tMPITag + 6, tTMatrixOffsets, tProcsReceivedFrom1 );

            // wait
            barrier();

            // receive
            Cell< Matrix< DDRMat > >   tRequestedTMatrixWeights;
            Cell< Matrix< IndexMat > > tRequestedTMatrixIndices;
            Cell< Matrix< IndexMat > > tRequestedTMatrixOwners;
            Cell< Matrix< IndexMat > > tRequestedTMatrixOffsets;

            // receive the answers
            mXTKModel->inward_receive_request_answers_reals( tMPITag + 3, 1, tProcRanks, tRequestedTMatrixWeights );
            mXTKModel->inward_receive_request_answers( tMPITag + 4, 1, tProcRanks, tRequestedTMatrixIndices );
            mXTKModel->inward_receive_request_answers( tMPITag + 5, 1, tProcRanks, tRequestedTMatrixOwners );
            mXTKModel->inward_receive_request_answers( tMPITag + 6, 1, tProcRanks, tRequestedTMatrixOffsets );

            barrier();

            // commit it to my data
            this->handle_received_interpolation_data(
                    tMeshIndex,
                    tNotOwnedIPVertIndsToProcs,
                    tNotOwnedEnrichedCellBulkPhaseToProcs,
                    tRequestedTMatrixWeights,
                    tRequestedTMatrixIndices,
                    tRequestedTMatrixOwners,
                    tRequestedTMatrixOffsets );

            // wait
            barrier();
        }
    }

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
                // get the master ip cell
                moris::mtk::Cell const & tMasterIpCell = tDblSideSetClusters( iC )->get_interpolation_cell( mtk::Master_Slave::MASTER );

                // get the slave ip cell
                moris::mtk::Cell const & tSlaveIpCell = tDblSideSetClusters( iC )->get_interpolation_cell( mtk::Master_Slave::SLAVE );

                // get the vertices attached to master/slave cells
                moris::Cell< mtk::Vertex* > tMasterVertices = tMasterIpCell.get_vertex_pointers();
                moris::Cell< mtk::Vertex* > tSlaveVertices  = tSlaveIpCell.get_vertex_pointers();

                // iterate through master vertices and place them in the correct list
                for ( uint iV = 0; iV < tMasterVertices.size(); iV++ )
                {
                    // get the ID of the current unzipped vertex
                    moris_index tVertexId = tMasterVertices( iV )->get_id();

                    // find out whether a T-matrix exists wrt. each B-spline mesh
                    bool tHasInterpolation = true;
                    for ( uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
                    {
                        tHasInterpolation = ( tHasInterpolation && tMasterVertices( iV )->has_interpolation( mMeshIndices( iBspMesh ) ) );
                    }

                    // add to vertices without interpolation
                    if ( !tHasInterpolation )
                    {
                        if ( tGhostVerticesWithOutInterpolationMap.find( tVertexId ) == tGhostVerticesWithOutInterpolationMap.end() )
                        {
                            tGhostVerticesWithOutInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithoutInterpolation.push_back( tMasterVertices( iV ) );
                            aGhostIpCellConnectedToVertex.push_back( &tMasterIpCell );
                        }
                    }

                    // add to vertices with interpolation
                    else if ( tHasInterpolation )
                    {
                        if ( tGhostVerticesWithInterpolationMap.find( tVertexId ) == tGhostVerticesWithInterpolationMap.end() )
                        {
                            tGhostVerticesWithInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithInterpolation.push_back( tMasterVertices( iV ) );
                        }
                    }
                }

                // iterate through slave vertices and place them in the correct list
                for ( uint iV = 0; iV < tSlaveVertices.size(); iV++ )
                {
                    // get the ID of the current unzipped vertex
                    moris_index tVertexId = tSlaveVertices( iV )->get_id();

                    // find out whether a T-matrix exists wrt. each B-spline mesh
                    bool tHasInterpolation = true;
                    for ( uint iBspMesh = 0; iBspMesh < mMeshIndices.numel(); iBspMesh++ )
                    {
                        tHasInterpolation = ( tHasInterpolation && tSlaveVertices( iV )->has_interpolation( mMeshIndices( iBspMesh ) ) );
                    }

                    // add to vertices without interpolation
                    if ( !tHasInterpolation )
                    {
                        if ( tGhostVerticesWithOutInterpolationMap.find( tVertexId ) == tGhostVerticesWithOutInterpolationMap.end() )
                        {
                            tGhostVerticesWithOutInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithoutInterpolation.push_back( tSlaveVertices( iV ) );
                            aGhostIpCellConnectedToVertex.push_back( &tSlaveIpCell );
                        }
                    }

                    // add to vertices with interpolation
                    else if ( tHasInterpolation )
                    {
                        if ( tGhostVerticesWithInterpolationMap.find( tVertexId ) == tGhostVerticesWithInterpolationMap.end() )
                        {
                            tGhostVerticesWithInterpolationMap[ tVertexId ] = true;
                            aGhostVerticesWithInterpolation.push_back( tSlaveVertices( iV ) );
                        }
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_interpolation_vertex_t_matrix_requests(
            moris::Cell< mtk::Vertex* >&              aGhostVerticesWithoutInterpolation,
            moris::Cell< mtk::Cell const * >&         aGhostIpCellConnectedToVertex,
            Cell< Matrix< IndexMat > >&               aNotOwnedIPVertIndsToProcs,
            Cell< Matrix< IndexMat > >&               aNotOwnedBGIPVertsIdsToProcs,
            Cell< Matrix< IndexMat > >&               aNotOwnedIpCellIdToProcs,
            Cell< Matrix< IndexMat > >&               aNotOwnedEnrichedCellBulkPhaseToProcs,
            Cell< uint >&                             aProcRanks,
            std::unordered_map< moris_id, moris_id >& aProcRankToDataIndex )
    {
        // access the enriched interpolation mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        Cell< Interpolation_Cell_Unzipped* >& tEnrIpCells = tEnrInterpMesh.get_enriched_interpolation_cells();

        // Counter and current proc index
        Cell< moris_id > tCounts( 0 );

        // temporary cell of cells which will be converted to matrices later
        Cell< Cell< moris_id > > tNotOwnedIPVertIndsToProcs;
        Cell< Cell< moris_id > > tNotOwnedBGIPVertsIdsToProcs;
        Cell< Cell< moris_id > > tNotOwnedIpCellIdToProcs;
        Cell< Cell< moris_id > > tNotOwnedIpCellBulkPhase;

        // get the communication table
        Matrix< IndexMat > tCommTable = tEnrInterpMesh.get_communication_table();

        // size
        aProcRanks.resize( tCommTable.numel() );

        for ( uint i = 0; i < tCommTable.numel(); i++ )
        {
            aProcRankToDataIndex[ tCommTable( i ) ] = i;
            aProcRanks( i )                         = ( tCommTable( i ) );

            // add cell for verts
            tNotOwnedIPVertIndsToProcs.push_back( Cell< moris_id >( 0 ) );
            tNotOwnedBGIPVertsIdsToProcs.push_back( Cell< moris_id >( 0 ) );
            tNotOwnedIpCellIdToProcs.push_back( Cell< moris_id >( 0 ) );
            tNotOwnedIpCellBulkPhase.push_back( Cell< moris_id >( 0 ) );
        }

        for ( uint iV = 0; iV < aGhostVerticesWithoutInterpolation.size(); iV++ )
        {
            // get the background xtk vertex
            Interpolation_Vertex_Unzipped* tXTKIpVert = tEnrInterpMesh.get_unzipped_vertex_pointer( aGhostVerticesWithoutInterpolation( iV )->get_index() );

            // get the xtk cell
            Interpolation_Cell_Unzipped* tEnrIpCell = tEnrIpCells( aGhostIpCellConnectedToVertex( iV )->get_index() );

            // unzipped cell owner always knows the T-matrices for all the vertices
            moris_index tOwnerProc = tEnrIpCell->get_owner();

            // get the index of this proc
            auto tProcIndexInData = aProcRankToDataIndex.find( tOwnerProc );

// check that the vertices are found
#ifdef MORIS_HAVE_DEBUG
            if( tProcIndexInData == aProcRankToDataIndex.end() )
            {
                std::cout << "Proc #" << par_rank() << ": Trying to communicate T-matrices on unzipped vertex which should be constructed by proc #" << 
                        tOwnerProc << " which is not part of the communication list." << std::endl;
                std::cout << "UIPV #" << tXTKIpVert->get_index() << ", ID: " << tXTKIpVert->get_id() << std::endl;
                moris::print_as_row_vector( tXTKIpVert->get_coords(), "Coords" );
                std::cout << "Vertex is part of UIPC #" << tEnrIpCell->get_index() << ", ID: " << tEnrIpCell->get_id() << 
                        " owned by proc #" << tEnrIpCell->get_owner() << std::endl;

                MORIS_ASSERT( false, "Ghost_Stabilization::prepare_interpolation_vertex_t_matrix_requests() - "
                        "Error in T-matrix communication. See information above." );
            }
#endif

            tNotOwnedIPVertIndsToProcs( tProcIndexInData->second ).push_back( tXTKIpVert->get_index() );
            tNotOwnedBGIPVertsIdsToProcs( tProcIndexInData->second ).push_back( tXTKIpVert->get_base_vertex()->get_id() );

            tNotOwnedIpCellIdToProcs( tProcIndexInData->second ).push_back( tEnrIpCell->get_id() );
            tNotOwnedIpCellBulkPhase( tProcIndexInData->second ).push_back( tEnrIpCell->get_bulkphase_index() );
        }

        // populate matrix in input data
        aNotOwnedIPVertIndsToProcs.clear();
        aNotOwnedIPVertIndsToProcs.resize( tNotOwnedIPVertIndsToProcs.size() );
        aNotOwnedBGIPVertsIdsToProcs.resize( tNotOwnedBGIPVertsIdsToProcs.size() );
        aNotOwnedIpCellIdToProcs.resize( tNotOwnedIpCellIdToProcs.size() );
        aNotOwnedEnrichedCellBulkPhaseToProcs.resize( tNotOwnedIpCellBulkPhase.size() );

        for ( uint iD = 0; iD < tNotOwnedIPVertIndsToProcs.size(); iD++ )
        {
            aNotOwnedIPVertIndsToProcs( iD ).resize( 1, tNotOwnedIPVertIndsToProcs( iD ).size() );
            aNotOwnedBGIPVertsIdsToProcs( iD ).resize( 1, tNotOwnedBGIPVertsIdsToProcs( iD ).size() );
            aNotOwnedIpCellIdToProcs( iD ).resize( 1, tNotOwnedIpCellIdToProcs( iD ).size() );
            aNotOwnedEnrichedCellBulkPhaseToProcs( iD ).resize( 1, tNotOwnedIpCellBulkPhase( iD ).size() );

            for ( uint jD = 0; jD < tNotOwnedIPVertIndsToProcs( iD ).size(); jD++ )
            {
                aNotOwnedIPVertIndsToProcs( iD )( jD )            = tNotOwnedIPVertIndsToProcs( iD )( jD );
                aNotOwnedBGIPVertsIdsToProcs( iD )( jD )          = tNotOwnedBGIPVertsIdsToProcs( iD )( jD );
                aNotOwnedIpCellIdToProcs( iD )( jD )              = tNotOwnedIpCellIdToProcs( iD )( jD );
                aNotOwnedEnrichedCellBulkPhaseToProcs( iD )( jD ) = tNotOwnedIpCellBulkPhase( iD )( jD );
            }

            if ( tNotOwnedIPVertIndsToProcs( iD ).size() == 0 )
            {
                aNotOwnedIPVertIndsToProcs( iD ).resize( 1, 1 );
                aNotOwnedBGIPVertsIdsToProcs( iD ).resize( 1, 1 );
                aNotOwnedIpCellIdToProcs( iD ).resize( 1, 1 );
                aNotOwnedEnrichedCellBulkPhaseToProcs( iD ).resize( 1, 1 );
                aNotOwnedIPVertIndsToProcs( iD )( 0 )            = MORIS_INDEX_MAX;
                aNotOwnedBGIPVertsIdsToProcs( iD )( 0 )          = MORIS_INDEX_MAX;
                aNotOwnedIpCellIdToProcs( iD )( 0 )              = MORIS_INDEX_MAX;
                aNotOwnedEnrichedCellBulkPhaseToProcs( iD )( 0 ) = MORIS_INDEX_MAX;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_t_matrix_request_answers(
            moris_index const &                aMeshIndex,
            Cell< Matrix< IndexMat > > const & aRequestedBgVertexIds,
            Cell< Matrix< IndexMat > > const & aRequestedIpCellIds,
            Cell< Matrix< IndexMat > > const & aIpCellBulkPhases,
            Cell< Matrix< DDRMat > >&          aTMatrixWeights,
            Cell< Matrix< IndexMat > >&        aTMatrixIndices,
            Cell< Matrix< IndexMat > >&        aBasisOwners,
            Cell< Matrix< IndexMat > >&        aTMatrixOffsets )
    {
        // access enriched ip mesh
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        MORIS_ASSERT( aRequestedBgVertexIds.size() == aRequestedIpCellIds.size(), "Mismatch in received communication information." );

        // information about size of interpolation mats
        Cell< uint > tSizes( aRequestedBgVertexIds.size(), 0 );

        // resize the input data
        aTMatrixWeights.resize( aRequestedBgVertexIds.size() );
        aTMatrixIndices.resize( aRequestedBgVertexIds.size() );
        aBasisOwners.resize( aRequestedBgVertexIds.size() );
        aTMatrixOffsets.resize( aRequestedBgVertexIds.size() );

        // collect the vertex interpolations
        Cell< Cell< Vertex_Enrichment* > > tVertexInterpolations( aRequestedBgVertexIds.size() );

        // collect size information throughout loop
        Cell< moris_index > tDataSizes( aRequestedBgVertexIds.size(), 0 );

        // iterate through and figure out how big to make the weights and indices mats
        // also collect vertex interpolations
        for ( uint iP = 0; iP < aRequestedBgVertexIds.size(); iP++ )
        {
            // no information requested
            if ( aRequestedBgVertexIds( iP ).numel() == 1 and aRequestedBgVertexIds( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                aTMatrixWeights( iP ).resize( 1, 1 );
                aTMatrixIndices( iP ).resize( 1, 1 );
                aBasisOwners( iP ).resize( 1, 1 );

                aTMatrixWeights( iP )( 0 ) = MORIS_REAL_MAX;
                aTMatrixIndices( iP )( 0 ) = MORIS_INDEX_MAX;
                aBasisOwners( iP )( 0 )    = MORIS_INDEX_MAX;
                continue;
            }

            // size the tmatrix offset for each vertex requested (num verts +1)
            aTMatrixOffsets( iP ).resize( 1, aRequestedBgVertexIds( iP ).numel() + 1 );

            // set the first one to 0
            aTMatrixOffsets( iP )( 0 ) = 0;

            // iterate through the vertices and get their interpolations and figure out
            // how big it is
            for ( uint iV = 0; iV < aRequestedBgVertexIds( iP ).numel(); iV++ )
            {
                // check that the bulk phases are consistent
                moris_index tCellIndex = tEnrInterpMesh.get_loc_entity_ind_from_entity_glb_id( aRequestedIpCellIds( iP )( iV ), EntityRank::ELEMENT, 0 );

                // get the cell
                Interpolation_Cell_Unzipped* tEnrIpCell = tEnrInterpMesh.get_enriched_interpolation_cells()( tCellIndex );

                // verifythat the bulk phases are consistent across procs
                MORIS_ERROR( tEnrIpCell->get_bulkphase_index() == aIpCellBulkPhases( iP )( iV ), "Parallel bulkphase mismatch." );

                // get the vertex
                moris_index tVertexIndex = this->get_enriched_interpolation_vertex( aRequestedBgVertexIds( iP )( iV ), aRequestedIpCellIds( iP )( iV ) );

                // get the vertex interpolation
                Interpolation_Vertex_Unzipped* tVertex = tEnrInterpMesh.get_unzipped_vertex_pointer( tVertexIndex );

                // get the vertex interpolation
                Vertex_Enrichment* tVertexInterp = tVertex->get_xtk_interpolation( aMeshIndex );

                MORIS_ASSERT( tVertexInterp->get_base_vertex_interpolation() != nullptr, "Owning proc has a nullptr for the vertex interpolation." );

                tVertexInterpolations( iP ).push_back( tVertexInterp );

                // number of basis functions interpolating into this vertex
                moris_index tNumBasis = tVertexInterp->get_basis_indices().numel();

                // offsets
                aTMatrixOffsets( iP )( iV + 1 ) = aTMatrixOffsets( iP )( iV ) + tNumBasis;

                // add to size
                tDataSizes( iP ) = tDataSizes( iP ) + tNumBasis;
            }
        }

        //  iterate through and size data
        for ( uint iP = 0; iP < aRequestedBgVertexIds.size(); iP++ )
        {
            if ( aRequestedBgVertexIds( iP ).numel() == 1 and aRequestedBgVertexIds( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                continue;
            }

            aTMatrixWeights( iP ).resize( 1, tDataSizes( iP ) );
            aTMatrixIndices( iP ).resize( 1, tDataSizes( iP ) );
            aBasisOwners( iP ).resize( 1, tDataSizes( iP ) );
        }

        // populate the data
        for ( uint iP = 0; iP < aRequestedBgVertexIds.size(); iP++ )
        {

            if ( aRequestedBgVertexIds( iP ).numel() == 1 and aRequestedBgVertexIds( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                aTMatrixOffsets( iP ).resize( 1, 1 );
                aTMatrixOffsets( iP )( 0 ) = MORIS_INDEX_MAX;
                continue;
            }

            uint tCount = 0;

            for ( uint iV = 0; iV < tVertexInterpolations( iP ).size(); iV++ )
            {
                this->add_vertex_interpolation_to_communication_data(
                        tCount,
                        tVertexInterpolations( iP )( iV ),
                        aTMatrixWeights( iP ),
                        aTMatrixIndices( iP ),
                        aBasisOwners( iP ),
                        aTMatrixOffsets( iP ) );
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::handle_received_interpolation_data(
            moris_index const &                aMeshIndex,
            Cell< Matrix< IndexMat > > const & aNotOwnedIPVertIndsToProcs,
            Cell< Matrix< IndexMat > > const & aNotOwnedEnrichedCellBulkPhaseToProcs,
            Cell< Matrix< DDRMat > > const &   aRequestedTMatrixWeights,
            Cell< Matrix< IndexMat > > const & aRequestedTMatrixIndices,
            Cell< Matrix< IndexMat > > const & aRequestedBasisOwners,
            Cell< Matrix< IndexMat > > const & aRequestedTMatrixOffsets )
    {
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // access the communication
        Matrix< IdMat > tCommTable = tEnrInterpMesh.get_communication_table();

        std::unordered_map< moris_id, moris_id > tProcRankToIndexInData;

        uint tCount = tCommTable.numel();

        // resize proc ranks and setup map to comm table
        for ( uint i = 0; i < tCommTable.numel(); i++ )
        {
            tProcRankToIndexInData[ tCommTable( i ) ] = i;
        }

        // iterate through returned information
        for ( uint iP = 0; iP < aNotOwnedIPVertIndsToProcs.size(); iP++ )
        {

            if ( aNotOwnedIPVertIndsToProcs( iP ).numel() == 1 and aNotOwnedIPVertIndsToProcs( iP )( 0 ) == MORIS_INDEX_MAX )
            {
                // do nothing for this iP
            }

            // standard case
            else
            {

                // extract the t-matrices and basis ids/owners for the proc ip
                Cell< Matrix< DDRMat > >   tExtractedTMatrixWeights;
                Cell< Matrix< IndexMat > > tExtractedTMatrixIds;
                Cell< Matrix< IndexMat > > tExtractedTBasisOwners;

                this->extract_vertex_interpolation_from_communication_data(
                        aNotOwnedIPVertIndsToProcs( iP ).numel(),
                        aRequestedTMatrixWeights( iP ),
                        aRequestedTMatrixIndices( iP ),
                        aRequestedBasisOwners( iP ),
                        aRequestedTMatrixOffsets( iP ),
                        tExtractedTMatrixWeights,
                        tExtractedTMatrixIds,
                        tExtractedTBasisOwners );

                // verify consistent sizes
                MORIS_ASSERT( aNotOwnedIPVertIndsToProcs( iP ).numel() == tExtractedTMatrixWeights.size(), "Size mismatch in t-matrix weights." );
                MORIS_ASSERT( aNotOwnedIPVertIndsToProcs( iP ).numel() == tExtractedTMatrixIds.size(), "Size mismatch in t-matrix ids." );
                MORIS_ASSERT( aNotOwnedIPVertIndsToProcs( iP ).numel() == tExtractedTBasisOwners.size(), "Size mismatch in basis owners." );
                MORIS_ASSERT( aNotOwnedIPVertIndsToProcs( iP ).numel() == aNotOwnedEnrichedCellBulkPhaseToProcs( iP ).numel(), "Size mismatch in bulk phases." );

                // iterate through vertices and set their interpolation weights and basis ids
                for ( uint iV = 0; iV < aNotOwnedIPVertIndsToProcs( iP ).numel(); iV++ )
                {
                    // get the vertex
                    moris_index tVertexIndex = aNotOwnedIPVertIndsToProcs( iP )( iV );

                    Interpolation_Vertex_Unzipped& tVertex = tEnrInterpMesh.get_xtk_interp_vertex( tVertexIndex );

                    // get the enriched vertex interpolation
                    Vertex_Enrichment* tVertexInterp = tVertex.get_xtk_interpolation( aMeshIndex );

                    // iterate through basis functions and find local indices
                    moris::Matrix< IndexMat > tBasisIndices( tExtractedTMatrixIds( iV ).n_rows(), tExtractedTMatrixIds( iV ).n_cols() );

                    for ( uint iBs = 0; iBs < tExtractedTMatrixIds( iV ).numel(); iBs++ )
                    {
                        // basis id
                        moris_id tId = tExtractedTMatrixIds( iV )( iBs );

                        // add this basis to the mesh if it doesnt exists on the current partition
                        if ( !tEnrInterpMesh.basis_exists_on_partition( aMeshIndex, tId ) )
                        {
                            MORIS_ASSERT( tExtractedTBasisOwners( iV )( iBs ) != par_rank(), "Owned basis should already exist on partition." );

                            tEnrInterpMesh.add_basis_function( aMeshIndex, tId, tExtractedTBasisOwners( iV )( iBs ), aNotOwnedEnrichedCellBulkPhaseToProcs( iP )( iV ) );
                        }

                        tBasisIndices( iBs ) = tEnrInterpMesh.get_enr_basis_index_from_enr_basis_id( aMeshIndex, tId );

                        moris_id tBasisOwner = tExtractedTBasisOwners( iV )( iBs );

                        MORIS_ASSERT( tEnrInterpMesh.get_basis_owner( tBasisIndices( iBs ), aMeshIndex ) == tBasisOwner, "Ownership discrepency." );

                        // TODO: this check is non-sense in the case of SPG based enrichment
                        // MORIS_ASSERT( tEnrInterpMesh.get_basis_bulk_phase( tBasisIndices( iBs ), aMeshIndex ) == aNotOwnedEnrichedCellBulkPhaseToProcs( iP )( iV ), "Bulkphase discrepency." );

                        // if the basis has an owning proc that is not in the comm table, add it to the comm table
                        if ( tProcRankToIndexInData.find( tBasisOwner ) == tProcRankToIndexInData.end() && tBasisOwner != par_rank() )
                        {
                            tEnrInterpMesh.add_proc_to_comm_table( tBasisOwner );
                            tProcRankToIndexInData[ tBasisOwner ] = tCount;
                            tCount++;
                        }
                    }

                    // iterate through basis in the base vertex interpolation
                    uint tNumCoeffs = tExtractedTMatrixIds( iV ).numel();

                    // Setup the map in the basis function
                    std::unordered_map< moris::moris_index, moris::moris_index >& tVertEnrichMap = tVertexInterp->get_basis_map();

                    for ( uint iB = 0; iB < tNumCoeffs; iB++ )
                    {
                        moris::moris_index tBasisIndex = tBasisIndices( iB );

                        tVertEnrichMap[ tBasisIndex ] = iB;
                    }

                    // get the basis indices from the basis ids
                    tVertexInterp->add_basis_information( tBasisIndices, tExtractedTMatrixIds( iV ) );
                    tVertexInterp->add_basis_weights( tBasisIndices, tExtractedTMatrixWeights( iV ) );
                    tVertexInterp->add_basis_owners( tBasisIndices, tExtractedTBasisOwners( iV ) );
                    tVertexInterp->add_base_vertex_interpolation( nullptr );    // base vertex interpolation does not exists (other  proc)
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
        moris_index tCellIndex = tEnrInterpMesh.get_loc_entity_ind_from_entity_glb_id( aEnrichedIpCellId, EntityRank::ELEMENT, 0 );

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

        // inialize container with Ghost Set Names
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

        // inialize container with Ghost Set Names
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
    Ghost_Stabilization::add_vertex_interpolation_to_communication_data(
            uint&               aCount,
            Vertex_Enrichment*  aInterpolation,
            Matrix< DDRMat >&   aTMatrixWeights,
            Matrix< IndexMat >& aTMatrixIndices,
            Matrix< IndexMat >& aTMatrixOwners,
            Matrix< IndexMat >& aTMatrixOffsets )
    {
        // access the basis indices and weights
        moris::Matrix< moris::IndexMat > const & tBasisIndices = aInterpolation->get_basis_ids();
        moris::Matrix< moris::DDRMat > const &   tBasisWeights = aInterpolation->get_basis_weights();
        moris::Matrix< moris::IndexMat >         tBasisOwners  = aInterpolation->get_owners();

        for ( uint i = 0; i < tBasisIndices.numel(); i++ )
        {
            aTMatrixIndices( aCount ) = tBasisIndices( i );
            aTMatrixWeights( aCount ) = tBasisWeights( i );

            MORIS_ASSERT( tBasisOwners( i ) < par_size() || tBasisOwners( i ) > 0, "Bad ownership for basis function." );

            aTMatrixOwners( aCount ) = tBasisOwners( i );
            aCount++;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::extract_vertex_interpolation_from_communication_data(
            uint const &                aNumVerts,
            Matrix< DDRMat > const &    aTMatrixWeights,
            Matrix< IndexMat > const &  aTMatrixIndices,
            Matrix< IndexMat > const &  aTMatrixOwners,
            Matrix< IndexMat > const &  aTMatrixOffsets,
            Cell< Matrix< DDRMat > >&   aExtractedTMatrixWeights,
            Cell< Matrix< IndexMat > >& aExtractedTMatrixIndices,
            Cell< Matrix< IndexMat > >& aExtractedBasisOwners )
    {

        // size output data
        aExtractedTMatrixWeights.resize( aNumVerts );
        aExtractedTMatrixIndices.resize( aNumVerts );
        aExtractedBasisOwners.resize( aNumVerts );

        // current starting index
        moris_index tStart = 0;

        // extract the data into the cells
        for ( uint iV = 0; iV < aNumVerts; iV++ )
        {
            // number of basis interpolating into the vertex
            moris::moris_index tNumBasis = aTMatrixOffsets( iV + 1 ) - tStart;

            aExtractedTMatrixWeights( iV ).resize( tNumBasis, 1 );
            aExtractedTMatrixIndices( iV ).resize( tNumBasis, 1 );
            aExtractedBasisOwners( iV ).resize( 1, tNumBasis );

            // itere and grab  data
            for ( moris::moris_index iIp = 0; iIp < tNumBasis; iIp++ )
            {
                aExtractedTMatrixWeights( iV )( iIp ) = aTMatrixWeights( tStart + iIp );
                aExtractedTMatrixIndices( iV )( iIp ) = aTMatrixIndices( tStart + iIp );
                aExtractedBasisOwners( iV )( iIp )    = aTMatrixOwners( tStart + iIp );
            }

            tStart = aTMatrixOffsets( iV + 1 );
        }
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
        aGhostSetupData.mMasterSideIpCells.reserve( tReserveSize );
        aGhostSetupData.mSlaveSideIpCells.reserve( tReserveSize );
        aGhostSetupData.mMasterSideIgCellSideOrds.reserve( tReserveSize );
        aGhostSetupData.mSlaveSideIgCellSideOrds.reserve( tReserveSize );
        aGhostSetupData.mNonTrivialFlag.reserve( tReserveSize );
        aGhostSetupData.mTransitionLocation.reserve( tReserveSize );

        aGhostSetupData.mMasterSideIpCells.resize( tNumBulkPhases );
        aGhostSetupData.mSlaveSideIpCells.resize( tNumBulkPhases );
        aGhostSetupData.mMasterSideIgCellSideOrds.resize( tNumBulkPhases );
        aGhostSetupData.mSlaveSideIgCellSideOrds.resize( tNumBulkPhases );
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
                    aGhostSetupData.mMasterSideIpCells( tBulkPhase ).push_back( tFirstInterpCell->get_index() );
                    aGhostSetupData.mSlaveSideIpCells( tBulkPhase ).push_back( tSecondInterpCell->get_index() );

                    // setup ig cells in ghost set up data
                    aGhostSetupData.mMasterSideIgCellSideOrds( tBulkPhase ).push_back( ( *tSubphaseToSubphaseMySideOrds( iSP ) )( jSpNeighbor ) );
                    aGhostSetupData.mSlaveSideIgCellSideOrds( tBulkPhase ).push_back( ( *tSubphaseToSubphaseNeighborSideOrds( iSP ) )( jSpNeighbor ) );

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
        MORIS_ASSERT( aGhostSetupData.mMasterSideIpCells.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mMasterSideIpCells too small, increase by %f\n",
                aGhostSetupData.mMasterSideIpCells.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mSlaveSideIpCells.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mSlaveSideIpCells too small, increase by %f\n",
                aGhostSetupData.mSlaveSideIpCells.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mMasterSideIgCellSideOrds.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mMasterSideIgCellSideOrds too small, increase by %f\n",
                aGhostSetupData.mMasterSideIgCellSideOrds.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mSlaveSideIgCellSideOrds.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mSlaveSideIgCellSideOrds too small, increase by %f\n",
                aGhostSetupData.mSlaveSideIgCellSideOrds.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mNonTrivialFlag.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mTrivialFlag too small, increase by %f\n",
                aGhostSetupData.mNonTrivialFlag.size() / tReserveSize );

        MORIS_ASSERT( aGhostSetupData.mTransitionLocation.size() < tReserveSize,
                "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of mTransitionLocation too small, increase by %f\n",
                aGhostSetupData.mTransitionLocation.size() / tReserveSize );

        // free up unused memory
        aGhostSetupData.mMasterSideIpCells.shrink_to_fit();
        aGhostSetupData.mSlaveSideIpCells.shrink_to_fit();
        aGhostSetupData.mMasterSideIgCellSideOrds.shrink_to_fit();
        aGhostSetupData.mSlaveSideIgCellSideOrds.shrink_to_fit();
        aGhostSetupData.mNonTrivialFlag.shrink_to_fit();
        aGhostSetupData.mTransitionLocation.shrink_to_fit();

        // allocate ids for non-trivial integration cells
        moris_id tCurrentId = 0;

        // reserve global entity ids for additional elements that need to be created to get dbl. sided facets ...
        if ( aGhostSetupData.mLinearBackgroundMesh )
        {
            // ... for non-trivial element transitions only in the linear case
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tNonTrivialCount, EntityRank::ELEMENT );
        }
        else
        {
            // ... // TODO: I don't understand this bit
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tEnrIntegMesh.get_num_entities( EntityRank::ELEMENT ), EntityRank::ELEMENT );
        }

        moris_id tCurrentIndex = tEnrIntegMesh.get_num_entities( EntityRank::ELEMENT );

        // get total number of ghost facets for each bulk phase on current processor
        Matrix< DDUMat > tLocalNumberOfGhostFacets( aGhostSetupData.mMasterSideIpCells.size(), 1 );
        for ( uint iBulkPhase = 0; iBulkPhase < aGhostSetupData.mMasterSideIpCells.size(); iBulkPhase++ )
        {
            tLocalNumberOfGhostFacets( iBulkPhase ) = aGhostSetupData.mMasterSideIpCells( iBulkPhase ).size();
        }

        // get global (i.e. across all procs) number of ghost facets for each bulk phase
        Matrix< DDUMat > tTotalNumberOfGhostFacets = sum_all_matrix( tLocalNumberOfGhostFacets );

        // build list of double side set indices
        moris::Cell< moris_index > tDoubleSideSetIndexList;
        tDoubleSideSetIndexList.reserve( aGhostSetupData.mMasterSideIpCells.size() );

        // ----------------------------------------------------------------------------------
        // loop constructing the actual ghost facets
        // iterate through bulk phases
        for ( uint iBulkPhase = 0; iBulkPhase < aGhostSetupData.mMasterSideIpCells.size(); iBulkPhase++ )
        {
            // allocate space in the integration mesh double side sets
            tEnrIntegMesh.mDoubleSideSets( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) ).resize( aGhostSetupData.mMasterSideIpCells( iBulkPhase ).size() );

            MORIS_LOG_SPEC( "Total Ghost Facets for Bulk Phase " + std::to_string( iBulkPhase ), tTotalNumberOfGhostFacets( iBulkPhase ) );

            // iterate through double sides in this bulk phase
            for ( uint iGhostFacet = 0; iGhostFacet < aGhostSetupData.mMasterSideIpCells( iBulkPhase ).size(); iGhostFacet++ )
            {
                // create a new side cluster for each of the pairs
                std::shared_ptr< Side_Cluster > tSlaveSideCluster =
                        this->create_slave_side_cluster( aGhostSetupData, tEnrIpCells, iBulkPhase, iGhostFacet, tCurrentIndex, tCurrentId );

                std::shared_ptr< Side_Cluster > tMasterSideCluster =
                        this->create_master_side_cluster( aGhostSetupData, tEnrIpCells, iBulkPhase, iGhostFacet, tSlaveSideCluster.get(), tCurrentIndex, tCurrentId );

                // verify the subphase cluster
                MORIS_ASSERT( tSlaveSideCluster->mInterpolationCell->get_bulkphase_index() == (moris_index)iBulkPhase,
                        "Bulk phase mismatch on slave side of double side set cluster" );

                MORIS_ASSERT( tMasterSideCluster->mInterpolationCell->get_bulkphase_index() == (moris_index)iBulkPhase,
                        "Bulk phase mismatch on master side of double side set cluster" );

                // add to side clusters the integration mesh
                tEnrIntegMesh.mDoubleSideSetsMasterIndex( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() );

                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tMasterSideCluster );

                tEnrIntegMesh.mDoubleSideSetsSlaveIndex( aGhostSetupData.mDblSideSetIndexInMesh( iBulkPhase ) ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() );

                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tSlaveSideCluster );

                // create double side cluster
                std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster = std::make_shared< mtk::Double_Side_Cluster >(
                        tMasterSideCluster.get(),
                        tSlaveSideCluster.get(),
                        tMasterSideCluster->mVerticesInCluster );

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

        tEnrIntegMesh.commit_double_side_set( tDoubleSideSetIndexList );

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
        aGhostSetupData.mMasterSideIpCellsNew.resize( tNumBspMeshes );
        aGhostSetupData.mSlaveSideIpCellsNew.resize( tNumBspMeshes );
        aGhostSetupData.mMasterSideIgCellSideOrdsNew.resize( tNumBspMeshes );
        aGhostSetupData.mSlaveSideIgCellSideOrdsNew.resize( tNumBspMeshes );
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
            aGhostSetupData.mMasterSideIpCellsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mSlaveSideIpCellsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mMasterSideIgCellSideOrdsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mSlaveSideIgCellSideOrdsNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mNonTrivialFlagNew( iBspMesh ).reserve( tReserveSize );
            aGhostSetupData.mTransitionLocationNew( iBspMesh ).reserve( tReserveSize );

            aGhostSetupData.mMasterSideIpCellsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mSlaveSideIpCellsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mMasterSideIgCellSideOrdsNew( iBspMesh ).resize( tNumBulkPhases );
            aGhostSetupData.mSlaveSideIgCellSideOrdsNew( iBspMesh ).resize( tNumBulkPhases );
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

                    // sanity check that both SPGs have the same bulk-phase indes
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

                                // store in Ghost setup data the Slave UIPCs used and their side ordinals
                                aGhostSetupData.mSlaveSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tUipcIndex );
                                aGhostSetupData.mSlaveSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tMySideOrdinal );

                                // store in Ghost setup data the Master UIPCs used and their side ordinals
                                aGhostSetupData.mMasterSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborUipcIndex );
                                aGhostSetupData.mMasterSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborSideOrdinal );

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

                                    // store in Ghost setup data the Slave UIPCs used and their side ordinals
                                    aGhostSetupData.mMasterSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tUipcIndex );
                                    aGhostSetupData.mMasterSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tMySideOrdinal );

                                    // store in Ghost setup data the Master UIPCs used and their side ordinals
                                    aGhostSetupData.mSlaveSideIpCellsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborUipcIndex );
                                    aGhostSetupData.mSlaveSideIgCellSideOrdsNew( iBspMesh )( tBulkPhaseIndex ).push_back( tNeighborSideOrdinal );

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

                    }    // end if: mark elements for ghost side construction
                }        // end for: loop over neighboring SPGs
            }            // end for: loop over all SPGs

            // check that reserved size was appropriate
            MORIS_ASSERT( aGhostSetupData.mMasterSideIpCellsNew( iBspMesh ).size() < tReserveSize,
                    "Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh: initial reservation of aGhostSetupData too small, increased by %f\n",
                    aGhostSetupData.mMasterSideIpCellsNew( iBspMesh ).size() / tReserveSize );

            // free up unused memory
            aGhostSetupData.mMasterSideIpCellsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mSlaveSideIpCellsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mMasterSideIgCellSideOrdsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mSlaveSideIgCellSideOrdsNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mNonTrivialFlagNew( iBspMesh ).shrink_to_fit();
            aGhostSetupData.mTransitionLocationNew( iBspMesh ).shrink_to_fit();
        }

        // allocate ids for non-trivial integration cells
        moris_id tCurrentId = 0;

        // reserve global entity ids for additional elements that need to be created to get dbl. sided facets ...
        if ( aGhostSetupData.mLinearBackgroundMesh )
        {
            // ... for non-trivial element transitions only in the linear case
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tNonTrivialCount, EntityRank::ELEMENT );
        }
        else
        {
            // ... ?
            tCurrentId = tEnrIntegMesh.allocate_entity_ids( tEnrIntegMesh.get_num_entities( EntityRank::ELEMENT ), EntityRank::ELEMENT );
        }

        // Get next free index for new cells/elements
        moris_id tCurrentIndex = tEnrIntegMesh.get_num_entities( EntityRank::ELEMENT );

        // get total number of ghost facets for each B-spline mesh and bulk phase on current processor
        Matrix< DDUMat > tLocalNumberOfGhostFacets( tNumBspMeshes, tNumBulkPhases, 1 );
        for ( uint iBspMesh = 0; iBspMesh < tNumBspMeshes; iBspMesh++ )
        {
            for ( uint iBulkPhase = 0; iBulkPhase < tNumBulkPhases; iBulkPhase++ )
            {
                tLocalNumberOfGhostFacets( iBspMesh, iBulkPhase ) = aGhostSetupData.mMasterSideIpCellsNew( iBspMesh )( iBulkPhase ).size();
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
                uint tNumGhostFacetsInDblSS = aGhostSetupData.mMasterSideIpCellsNew( iBspMesh )( iBulkPhase ).size();

                // resize list of Double sided clusters for current Dbl.-SS. in enr. IG mesh
                tEnrIntegMesh.mDoubleSideSets( tCurrentDblSsIndexInEnrIgMesh ).resize( tNumGhostFacetsInDblSS );

                // print number of Ghost facets to console
                MORIS_LOG_SPEC(
                        "Total Ghost Facets for B-spline mesh " + std::to_string( mMeshIndices( iBspMesh ) ) + " Bulk Phase " + std::to_string( iBulkPhase ),
                        tLocalNumberOfGhostFacets( iBspMesh, iBulkPhase ) );

                // iterate through double sides in this bulk phase
                for ( uint iGhostFacet = 0; iGhostFacet < tNumGhostFacetsInDblSS; iGhostFacet++ )
                {
                    // create a new single side cluster - slave side
                    std::shared_ptr< Side_Cluster > tSlaveSideCluster =
                            this->create_slave_side_cluster_new(
                                    aGhostSetupData,
                                    tEnrIpCells,
                                    iBspMesh,
                                    iBulkPhase,
                                    iGhostFacet,
                                    tCurrentIndex,
                                    tCurrentId );

                    // create a new single side cluster - master side
                    std::shared_ptr< Side_Cluster > tMasterSideCluster =
                            this->create_master_side_cluster_new(
                                    aGhostSetupData,
                                    tEnrIpCells,
                                    iBspMesh,
                                    iBulkPhase,
                                    iGhostFacet,
                                    tSlaveSideCluster.get(),
                                    tCurrentIndex,
                                    tCurrentId );

                    // add single side cluster for master side to the integration mesh
                    tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tMasterSideCluster );

                    // store the index of the single side clusters associated with the double side clusters
                    tEnrIntegMesh.mDoubleSideSetsMasterIndex( tCurrentDblSsIndexInEnrIgMesh ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() - 1 );

                    // add single side cluster for slave side to the integration mesh
                    tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back( tSlaveSideCluster );

                    // store the index of the single side clusters associated with the double side clusters
                    tEnrIntegMesh.mDoubleSideSetsSlaveIndex( tCurrentDblSsIndexInEnrIgMesh ).push_back( tEnrIntegMesh.mDoubleSideSingleSideClusters.size() - 1 );

                    // create double side cluster
                    std::shared_ptr< mtk::Double_Side_Cluster > tDblSideCluster = std::make_shared< mtk::Double_Side_Cluster >(
                            tMasterSideCluster.get(),
                            tSlaveSideCluster.get(),
                            tMasterSideCluster->mVerticesInCluster );

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
        // 2. The owning processor of the master (first) subphase constructs the ghost facet.
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
    Ghost_Stabilization::create_slave_side_cluster(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBulkIndex,
            uint const &                          aCellIndex,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // cut integration mesh access
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create a new side cluster for the slave
        std::shared_ptr< Side_Cluster > tSlaveSideCluster = std::make_shared< Side_Cluster >();

        // give the cluster the enriched interpolation cell
        tSlaveSideCluster->mInterpolationCell = aEnrIpCells( aGhostSetupData.mSlaveSideIpCells( aBulkIndex )( aCellIndex ) );

        // slave cluster is always trivial because the small facet is always the slave in the case of HMR hanging nodes
        tSlaveSideCluster->mTrivial = true;

        // create a linear IG cell over the whole Lagrange interpolation element and store this as being part of the new side cluster
        tSlaveSideCluster->mIntegrationCells = {
            this->get_linear_ig_cell(
                    aGhostSetupData,
                    aEnrIpCells( aGhostSetupData.mSlaveSideIpCells( aBulkIndex )( aCellIndex ) ),
                    aCurrentIndex,
                    aCurrentId )
        };

        // allocate space in integration cell side ordinals
        tSlaveSideCluster->mIntegrationCellSideOrdinals = { { aGhostSetupData.mSlaveSideIgCellSideOrds( aBulkIndex )( aCellIndex ) } };

        // add geometric vertices to the cluster // note: why is this
        // tSlaveSideCluster->mVerticesInCluster =
        //         tSlaveSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( tSlaveSideCluster->mIntegrationCellSideOrdinals( 0 ) );

        // add the corresponding cell (i.e. bulk) cluster to the new side cluster
        tSlaveSideCluster->mAssociatedCellCluster =
                &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tSlaveSideCluster->mInterpolationCell );

        // find and add the group of vertices to the new side cluster
        std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tSlaveSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );
        tSlaveSideCluster->set_ig_vertex_group( tVertexGroup );

        // return the shared pointer to the newly generated ghost slave side cluster
        return tSlaveSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_slave_side_cluster_new(
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

        // create a new side cluster for the slave
        std::shared_ptr< Side_Cluster > tSlaveSideCluster = std::make_shared< Side_Cluster >();

        // get the index to the Enr. IP cell the side cluster is attached to
        moris_index tSlaveEnrIpCellIndex = aGhostSetupData.mSlaveSideIpCellsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );

        // give the cluster the enriched interpolation cell
        tSlaveSideCluster->mInterpolationCell = aEnrIpCells( tSlaveEnrIpCellIndex );

        // slave cluster is always trivial because the small facet is always the slave in the case of HMR hanging nodes
        tSlaveSideCluster->mTrivial = true;

        // add integration cell
        tSlaveSideCluster->mIntegrationCells = {
            this->get_linear_ig_cell(
                    aGhostSetupData,
                    aEnrIpCells( tSlaveEnrIpCellIndex ),
                    aCurrentIndex,
                    aCurrentId )
        };

        // allocate space in integration cell side ordinals
        moris_index tSideOrdinal                        = aGhostSetupData.mSlaveSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );
        tSlaveSideCluster->mIntegrationCellSideOrdinals = { { tSideOrdinal } };

        // assign bulk cluster corresponding to new side cluster
        tSlaveSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tSlaveSideCluster->mInterpolationCell );

        // leader vertex group
        std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tSlaveSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

        tSlaveSideCluster->set_ig_vertex_group( tVertexGroup );

        return tSlaveSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_master_side_cluster(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBulkIndex,
            uint const &                          aCellIndex,
            Side_Cluster*                         aSlaveSideCluster,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // get access to the cut integration mesh
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create a master side cluster
        std::shared_ptr< Side_Cluster > tMasterSideCluster = std::make_shared< Side_Cluster >();

        // assign the UIPC to the new side cluster
        tMasterSideCluster->mInterpolationCell = aEnrIpCells( aGhostSetupData.mMasterSideIpCells( aBulkIndex )( aCellIndex ) );

        // Case: non-trivial Lagrange element transition
        if ( aGhostSetupData.mNonTrivialFlag( aBulkIndex )( aCellIndex ) > 0 )
        {
            // flag the master side as non-trivial
            tMasterSideCluster->mTrivial = false;

            // create new integration cell using the vertices on the slave facet and the adjacent vertices of the base interpolation cell
            std::shared_ptr< xtk::Cell_XTK_No_CM > tNewIgCell = this->create_non_trivial_master_ig_cell(
                    aGhostSetupData,
                    aBulkIndex,
                    aCellIndex,
                    tMasterSideCluster.get(),
                    aSlaveSideCluster,
                    aCurrentIndex,
                    aCurrentId );

            // get the local coordinates of the IP cell vertices including the hanging vertex
            Cell< Matrix< DDRMat > > tLocCoords;
            this->get_local_coords_on_transition_side(
                    aGhostSetupData.mMasterSideIgCellSideOrds( aBulkIndex )( aCellIndex ),
                    aGhostSetupData.mTransitionLocation( aBulkIndex )( aCellIndex ),
                    tLocCoords );

            // add integration cell
            tMasterSideCluster->mIntegrationCells.push_back( tNewIgCell.get() );

            // add side ordinal relative to the integration cell
            tMasterSideCluster->mIntegrationCellSideOrdinals = { { this->get_side_ordinals_for_non_trivial_master() } };

            // get the vertex group
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tMasterSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

            // Add the slave vertex to the vertex group
            moris::Cell< moris::mtk::Vertex const * > tSlaveVertices =
                    aSlaveSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aSlaveSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // iterate through vertices
            for ( uint iV = 0; iV < tSlaveVertices.size(); iV++ )
            {
                if ( !tVertexGroup->vertex_is_in_group( tSlaveVertices( iV )->get_index() ) )
                {
                    tVertexGroup->add_vertex( tSlaveVertices( iV ), std::make_shared< Matrix< DDRMat > >( tLocCoords( iV ) ) );
                }

                MORIS_ASSERT(
                        moris::norm( ( *tVertexGroup->get_vertex_local_coords( tSlaveVertices( iV )->get_index() ) - tLocCoords( iV ) ) ) < 1e-8,
                        "Ghost_Stabilization::create_master_side_cluster() - Local coord issue" );
            }

            // set the vertex group
            tMasterSideCluster->set_ig_vertex_group( tVertexGroup );

            // // get the vertices on the side ordinal
            // tMasterSideCluster->mVerticesInCluster = tSlaveVertices;

            // add the local coordinates
            tMasterSideCluster->mVertexLocalCoords = tLocCoords;

            // associated cell cluster
            tMasterSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tMasterSideCluster->mInterpolationCell );

            // verify new cluster
            mtk::Mesh_Checker tCheck;
            MORIS_ASSERT(
                    tCheck.verify_side_cluster( tMasterSideCluster.get(), mtk::Master_Slave::MASTER ),
                    "Ghost_Stabilization::create_master_side_cluster() - "
                    "Invalid Side Cluster Created Check the local coordinates" );

            // place the new ig cell in the background mesh
            mXTKModel->get_cut_integration_mesh()->add_integration_cell( tNewIgCell->get_index(), tNewIgCell );
        }

        // Case: trivial lagrange element transition
        else
        {
            // flag the master side as trivial
            tMasterSideCluster->mTrivial = true;

            // create a linear IG cell over the whole Lagrange interpolation element and store this as being part of the new side cluster
            tMasterSideCluster->mIntegrationCells = {
                this->get_linear_ig_cell(
                        aGhostSetupData,
                        aEnrIpCells( aGhostSetupData.mMasterSideIpCells( aBulkIndex )( aCellIndex ) ),
                        aCurrentIndex,
                        aCurrentId )
            };

            // add side ordinal on which the ghost facet sits relative to the master cluster
            moris_index tSideOrdinal                         = aGhostSetupData.mMasterSideIgCellSideOrds( aBulkIndex )( aCellIndex );
            tMasterSideCluster->mIntegrationCellSideOrdinals = { { tSideOrdinal } };

            // get and add the vertices on that side ordinal
            tMasterSideCluster->mVerticesInCluster =
                    tMasterSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( tMasterSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // get and set the associated Cell (i.e. bulk) cluster
            tMasterSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tMasterSideCluster->mInterpolationCell );

            // get and set the vertex group from the IP cell
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tMasterSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );
            tMasterSideCluster->set_ig_vertex_group( tVertexGroup );
        }

        // return the shared pointer to the master side cluster contstructed
        return tMasterSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< Side_Cluster >
    Ghost_Stabilization::create_master_side_cluster_new(
            Ghost_Setup_Data&                     aGhostSetupData,
            Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
            uint const &                          aBsplineMeshListIndex,
            uint const &                          aBulkPhaseIndex,
            uint const &                          aGhostFacetIndexInSet,
            Side_Cluster*                         aSlaveSideCluster,
            moris_index&                          aCurrentIndex,
            moris_index&                          aCurrentId )
    {
        // get access to the cut integration mesh
        Cut_Integration_Mesh* tCutIGMesh = mXTKModel->get_cut_integration_mesh();

        // create the master side cluster
        std::shared_ptr< Side_Cluster > tMasterSideCluster = std::make_shared< Side_Cluster >();

        // get the index to the Enr. IP cell the side cluster is attached to
        moris_index tMasterEnrIpCellIndex = aGhostSetupData.mMasterSideIpCellsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );

        // set the enr. IP cell the side cluster is attached to
        tMasterSideCluster->mInterpolationCell = aEnrIpCells( tMasterEnrIpCellIndex );

        // Setup the master side cluster
        if ( aGhostSetupData.mNonTrivialFlagNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet ) > 0 )
        {
            // flag the master side as non-trivial
            tMasterSideCluster->mTrivial = false;

            // create new integration cell using the vertices on the slave facet and the adjacent vertices of the base interpolation cell
            std::shared_ptr< xtk::Cell_XTK_No_CM > tNewIgCell = this->create_non_trivial_master_ig_cell_new(
                    aGhostSetupData,
                    aBsplineMeshListIndex,
                    aBulkPhaseIndex,
                    aGhostFacetIndexInSet,
                    tMasterSideCluster.get(),
                    aSlaveSideCluster,
                    aCurrentIndex,
                    aCurrentId );

            // get the local coordinates of the IP cell vertices including the hanging vertex
            Cell< Matrix< DDRMat > > tLocCoords;
            this->get_local_coords_on_transition_side(
                    aGhostSetupData.mMasterSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet ),
                    aGhostSetupData.mTransitionLocationNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet ),
                    tLocCoords );

            // add integration cell
            tMasterSideCluster->mIntegrationCells.push_back( tNewIgCell.get() );

            // add side ordinal relative to the integration cell
            tMasterSideCluster->mIntegrationCellSideOrdinals = { { this->get_side_ordinals_for_non_trivial_master() } };

            // get the vertex group
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tMasterSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

            // Add the slave vertex to the vertex group
            moris::Cell< moris::mtk::Vertex const * > tSlaveVertices =
                    aSlaveSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aSlaveSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // iterate through vertices
            for ( uint iV = 0; iV < tSlaveVertices.size(); iV++ )
            {
                if ( !tVertexGroup->vertex_is_in_group( tSlaveVertices( iV )->get_index() ) )
                {
                    tVertexGroup->add_vertex( tSlaveVertices( iV ), std::make_shared< Matrix< DDRMat > >( tLocCoords( iV ) ) );
                }

                MORIS_ASSERT(
                        moris::norm( ( *tVertexGroup->get_vertex_local_coords( tSlaveVertices( iV )->get_index() ) - tLocCoords( iV ) ) ) < 1e-8,
                        "Ghost_Stabilization::create_master_side_cluster() - Local coord issue" );
            }

            // set the vertex group
            tMasterSideCluster->set_ig_vertex_group( tVertexGroup );

            // add the local coordinates
            tMasterSideCluster->mVertexLocalCoords = tLocCoords;

            // associated cell cluster
            tMasterSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tMasterSideCluster->mInterpolationCell );

            // verify new cluster
            mtk::Mesh_Checker tCheck;
            MORIS_ASSERT(
                    tCheck.verify_side_cluster( tMasterSideCluster.get(), mtk::Master_Slave::MASTER ),
                    "Ghost_Stabilization::create_master_side_cluster_new() - "
                    "Invalid Side Cluster Created Check the local coordinates" );

            // place the new ig cell in the background mesh
            mXTKModel->get_cut_integration_mesh()->add_integration_cell( tNewIgCell->get_index(), tNewIgCell );
        }

        // Case: trivial Lagrange element transition
        else
        {
            // flag the master side as trivial
            tMasterSideCluster->mTrivial = true;

            // add integration cell
            tMasterSideCluster->mIntegrationCells = {
                this->get_linear_ig_cell(
                        aGhostSetupData,
                        aEnrIpCells( tMasterEnrIpCellIndex ),
                        aCurrentIndex,
                        aCurrentId )
            };

            // add side ordinal relative to the integration cell
            moris_index tSideOrdinal                         = aGhostSetupData.mMasterSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkPhaseIndex )( aGhostFacetIndexInSet );
            tMasterSideCluster->mIntegrationCellSideOrdinals = { { tSideOrdinal } };

            // add the vertices on the side ordinal
            tMasterSideCluster->mVerticesInCluster =
                    tMasterSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( tMasterSideCluster->mIntegrationCellSideOrdinals( 0 ) );

            // get and set the cell cluster
            tMasterSideCluster->mAssociatedCellCluster = &mXTKModel->get_enriched_integ_mesh().get_xtk_cell_cluster( *tMasterSideCluster->mInterpolationCell );

            // get the vertex group from the IP cell
            std::shared_ptr< IG_Vertex_Group > tVertexGroup =
                    tCutIGMesh->get_vertex_group( tCutIGMesh->get_parent_cell_group_index( tMasterSideCluster->mInterpolationCell->get_base_cell()->get_index() ) );

            // set Side cluster to use vertex group from IP cell
            tMasterSideCluster->set_ig_vertex_group( tVertexGroup );
        }

        // return side cluster contstructed
        return tMasterSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< xtk::Cell_XTK_No_CM >
    Ghost_Stabilization::create_non_trivial_master_ig_cell(
            Ghost_Setup_Data& aGhostSetupData,
            uint const &      aBulkIndex,
            uint const &      aCellIndex,
            Side_Cluster*     aMasterSideCluster,
            Side_Cluster*     aSlaveSideCluster,
            moris_index&      aCurrentIndex,
            moris_index&      aCurrentId )
    {
        // make sure the corresponding slave side cluster is valid
        MORIS_ASSERT( aSlaveSideCluster->get_cells_in_side_cluster().size() == 1,
                "Ghost_Stabilization::create_non_trivial_master_ig_cell() - "
                "Slave side cluster should have exactly one integration cell." );

        // get the vertices on the side for the slave side cluster
        moris::Cell< moris::mtk::Vertex const * > tSlaveVertices =
                aSlaveSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aSlaveSideCluster->mIntegrationCellSideOrdinals( 0 ) );

        // get the master UIPC
        Interpolation_Cell_Unzipped const * tMasterIpCell = aMasterSideCluster->mInterpolationCell;

        // base master IP cell
        moris::mtk::Cell const * tBaseMasterCell = tMasterIpCell->get_base_cell();

        // get the connectivity information from the cell
        std::shared_ptr< moris::mtk::Cell_Info > tCellInfo = tMasterIpCell->get_cell_info_sp();

        // adjacent side ordinal on master
        uint tAdjFacetOrd = tCellInfo->get_adjacent_side_ordinal(
                aGhostSetupData.mMasterSideIgCellSideOrds( aBulkIndex )( aCellIndex ) );

        // setup the vertices and local coordinates of the vertices relative to the cell
        moris::Cell< moris::mtk::Vertex const * > tAdjVertices =
                tBaseMasterCell->get_geometric_vertices_on_side_ordinal( tAdjFacetOrd );

        // line up the ordering of the master and slave side vertices such that they match
        moris::Cell< moris::mtk::Vertex const * > tPermutedSlaveVertices;
        moris::Cell< moris::mtk::Vertex const * > tPermutedAdjVertices;
        this->permute_slave_vertices(
                tSlaveVertices,
                tAdjVertices,
                tPermutedSlaveVertices,
                tPermutedAdjVertices );

        // collect all vertices on the ghost facet
        moris::Cell< moris::mtk::Vertex* > tCellVertices( tPermutedAdjVertices.size() + tPermutedSlaveVertices.size() );

        uint tCount = 0;
        for ( auto iV : tPermutedAdjVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }
        for ( auto iV : tPermutedSlaveVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }

        // linear cell info
        mtk::Cell_Info_Factory                   tCellInfoFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo =
                tCellInfoFactory.create_cell_info_sp( tMasterIpCell->get_geometry_type(), mtk::Interpolation_Order::LINEAR );

        // create a new integration cell that does not have a child mesh association
        std::shared_ptr< xtk::Cell_XTK_No_CM > tIgCell = std::make_shared< xtk::Cell_XTK_No_CM >(
                aCurrentId,
                aCurrentIndex,
                tMasterIpCell->get_owner(),
                tLinearCellInfo,
                tCellVertices );

        // increment current id and index
        aCurrentId++;
        aCurrentIndex++;

        return tIgCell;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr< xtk::Cell_XTK_No_CM >
    Ghost_Stabilization::create_non_trivial_master_ig_cell_new(
            Ghost_Setup_Data& aGhostSetupData,
            uint const &      aBsplineMeshListIndex,
            uint const &      aBulkIndex,
            uint const &      aCellIndex,
            Side_Cluster*     aMasterSideCluster,
            Side_Cluster*     aSlaveSideCluster,
            moris_index&      aCurrentIndex,
            moris_index&      aCurrentId )
    {
        // make sure the corresponding slave side cluster is valid
        MORIS_ASSERT( aSlaveSideCluster->get_cells_in_side_cluster().size() == 1,
                "Ghost_Stabilization::create_non_trivial_master_ig_cell() - "
                "Slave side cluster should have exactly one integration cell." );

        // get the vertices on the side for the slave side cluster
        moris::Cell< moris::mtk::Vertex const * > tSlaveVertices =
                aSlaveSideCluster->mIntegrationCells( 0 )->get_geometric_vertices_on_side_ordinal( aSlaveSideCluster->mIntegrationCellSideOrdinals( 0 ) );

        // get the master UIPC
        Interpolation_Cell_Unzipped const * tMasterIpCell = aMasterSideCluster->mInterpolationCell;

        // base master IP cell
        moris::mtk::Cell const * tBaseMasterCell = tMasterIpCell->get_base_cell();

        // get the connectivity information from the cell
        std::shared_ptr< moris::mtk::Cell_Info > tCellInfo = tMasterIpCell->get_cell_info_sp();

        // adjacent side ordinal on master
        uint tAdjFacetOrd = tCellInfo->get_adjacent_side_ordinal(
                aGhostSetupData.mMasterSideIgCellSideOrdsNew( aBsplineMeshListIndex )( aBulkIndex )( aCellIndex ) );

        // setup the vertices and local coordinates of the vertices relative to the cell
        moris::Cell< moris::mtk::Vertex const * > tAdjVertices =
                tBaseMasterCell->get_geometric_vertices_on_side_ordinal( tAdjFacetOrd );

        // line up the ordering of the master and slave side vertices such that they match
        moris::Cell< moris::mtk::Vertex const * > tPermutedSlaveVertices;
        moris::Cell< moris::mtk::Vertex const * > tPermutedAdjVertices;
        this->permute_slave_vertices(
                tSlaveVertices,
                tAdjVertices,
                tPermutedSlaveVertices,
                tPermutedAdjVertices );

        // collect all vertices on the ghost facet
        moris::Cell< moris::mtk::Vertex* > tCellVertices( tPermutedAdjVertices.size() + tPermutedSlaveVertices.size() );

        uint tCount = 0;
        for ( auto iV : tPermutedAdjVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }
        for ( auto iV : tPermutedSlaveVertices )
        {
            tCellVertices( tCount++ ) = &mXTKModel->get_background_mesh().get_mtk_vertex( iV->get_index() );
        }

        // linear cell info
        mtk::Cell_Info_Factory                   tCellInfoFactory;
        std::shared_ptr< moris::mtk::Cell_Info > tLinearCellInfo =
                tCellInfoFactory.create_cell_info_sp( tMasterIpCell->get_geometry_type(), mtk::Interpolation_Order::LINEAR );

        // create a new integration cell that does not have a child mesh association
        std::shared_ptr< xtk::Cell_XTK_No_CM > tIgCell = std::make_shared< xtk::Cell_XTK_No_CM >(
                aCurrentId,
                aCurrentIndex,
                tMasterIpCell->get_owner(),
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
    Ghost_Stabilization::permute_slave_vertices(
            moris::Cell< moris::mtk::Vertex const * > const & aSlaveVertices,
            moris::Cell< moris::mtk::Vertex const * > const & aMasterVertices,
            moris::Cell< moris::mtk::Vertex const * >&        aPermutedSlaveVertices,
            moris::Cell< moris::mtk::Vertex const * >&        aPermutedAdjMastVertices )
    {
        uint tSpatialDim = mXTKModel->get_spatial_dim();

        switch ( tSpatialDim )
        {
            case 2:
            {
                aPermutedSlaveVertices   = { aSlaveVertices( 1 ), aSlaveVertices( 0 ) };
                aPermutedAdjMastVertices = { aMasterVertices( 0 ), aMasterVertices( 1 ) };
                break;
            }
            default:
            {
                aPermutedSlaveVertices   = { aSlaveVertices( 0 ), aSlaveVertices( 3 ), aSlaveVertices( 2 ), aSlaveVertices( 1 ) };
                aPermutedAdjMastVertices = { aMasterVertices( 0 ), aMasterVertices( 3 ), aMasterVertices( 2 ), aMasterVertices( 1 ) };
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
    Ghost_Stabilization::get_side_ordinals_for_non_trivial_master()
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
                    "Ghost_Stabilization::get_side_ordinals_for_non_trivial_master() - "
                    "Invalid spatial dimension for ghost" );

            return 0;
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Ghost_Stabilization::is_linear_ip_mesh()
    {
        Enriched_Interpolation_Mesh& tEnrIpMesh = mXTKModel->get_enriched_interp_mesh( 0 );

        MORIS_ERROR( tEnrIpMesh.get_num_entities( EntityRank::ELEMENT ) > 0, "Cannot deduce type on empty mesh." );

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
