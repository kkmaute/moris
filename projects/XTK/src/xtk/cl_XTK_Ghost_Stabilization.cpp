/*
 * cl_XTK_Ghost_Stabilization.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "fn_generate_shared_face_element_graph.hpp"

namespace xtk
{
    Ghost_Stabilization::Ghost_Stabilization()
    :mXTKModel(nullptr)
    {}

    // ----------------------------------------------------------------------------------

    Ghost_Stabilization::Ghost_Stabilization( Model* aXTKModel )
    :mXTKModel(aXTKModel)
    {}

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::setup_ghost_stabilization()
    {
        Ghost_Setup_Data tGhostSetupData;

        //
        //        mtk::Mesh_Checker tMeshChecker(0,&mXTKModel->get_enriched_interp_mesh(0),&mXTKModel->get_enriched_integ_mesh(0));
        //        tMeshChecker.perform();
        //        tMeshChecker.print_diagnostics();
        //



        // construct trivial subphase interpolation cells
        this->construct_ip_ig_cells_for_ghost_side_clusters(tGhostSetupData);

        // identify and communicate vertices in aura that interpolate into ghost facets
        this->identify_and_setup_aura_vertices_in_ghost();

        // setup the side sets
        this->declare_ghost_double_side_sets_in_mesh(tGhostSetupData);

        // Construct Ghost Double Side Clusters and Sets
        this->construct_ghost_double_side_sets_in_mesh(tGhostSetupData);

        MORIS_ERROR(mXTKModel->get_enriched_interp_mesh(0).verify_basis_support(),"Issue detected in basis support after ghost stabilization."); // TODO: change to assert once done debugging

    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::visualize_ghost_on_mesh(moris_index aBulkPhase)
    {
        // get the enriched integration mesh
        Enriched_Integration_Mesh & tEnrIgMesh = mXTKModel->get_enriched_integ_mesh(0);

        // get the set name
        std::string tGhostName = this->get_ghost_dbl_side_set_name(aBulkPhase);

        // get the double side set index in the mesh
        moris_index tDSSIndexInMesh = tEnrIgMesh.get_double_sided_set_index(tGhostName);

        // topology
        enum CellTopology tFacetTopo = CellTopology::QUAD4;
        if(mXTKModel->get_spatial_dim() == 2)
        {
            tFacetTopo = CellTopology::QUAD4;
        }

        moris_index tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set(tDSSIndexInMesh,"ghost_ss_" + std::to_string(aBulkPhase));
        tEnrIgMesh.create_block_set_from_cells_of_side_set(tSSIndex,"ghost_bs_" + std::to_string(aBulkPhase), tFacetTopo);

        tEnrIgMesh.setup_color_to_set();
        tEnrIgMesh.collect_all_sets();
    }

    // ----------------------------------------------------------------------------------

    std::string
    Ghost_Stabilization::get_ghost_dbl_side_set_name(moris_index aBulkPhase)
    {
        MORIS_ASSERT(aBulkPhase < (moris_index)mXTKModel->get_geom_engine()->get_num_bulk_phase(),"Bulk Phase index out of bounds.");

        return "ghost_p" + std::to_string(aBulkPhase);
    }
    // ----------------------------------------------------------------------------------

    //FIXME: Keenan - code below needs to be split up in smaller units
    void
    Ghost_Stabilization::construct_ip_ig_cells_for_ghost_side_clusters(Ghost_Setup_Data & aGhostSetupData)
    {
        // access enriched integration mesh and enriched interp mesh
        Enriched_Integration_Mesh &   tEnrIgMesh = mXTKModel->get_enriched_integ_mesh(0);
        Enriched_Interpolation_Mesh & tEnrIpMesh = mXTKModel->get_enriched_interp_mesh(0);

        // get the groupings of interpolation cells
        Cell<Interpolation_Cell_Unzipped *>       tOwnedInterpCells;
        Cell<Cell<Interpolation_Cell_Unzipped *>> tNotOwnedInterpCells;
        Cell<uint>                                tProcRanks;
        tEnrIpMesh.get_owned_and_not_owned_enriched_interpolation_cells(tOwnedInterpCells,tNotOwnedInterpCells, tProcRanks);

        // number of interpolation cells
        moris::uint tCurrentNewInterpCellIndex = tEnrIpMesh.get_num_entities(EntityRank::ELEMENT);

        // allocate data in ghost setup data
        aGhostSetupData.mSubphaseIndexToInterpolationCellIndex.resize(mXTKModel->get_subphase_to_subphase().size(),MORIS_INDEX_MAX);

        // figure out how many of the interpolation cells are not trivial (this value is the number of new interpolion cells)
        // iterate through owned interpolation cells and give them an id if they are not trivial clusters
        // Also, we're only interested in the nontrivial cell clusters so keep track of those.
        moris::uint tNumNewInterpCellsNotOwned = 0;
        moris::uint tNumNewInterpCellsOwned    = 0;

        Cell<Interpolation_Cell_Unzipped *> tNonTrivialOwnedInterpCells;

        for(moris::size_t i = 0; i<tOwnedInterpCells.size(); i++)
        {
            xtk::Cell_Cluster const & tCluster = tEnrIgMesh.get_xtk_cell_cluster(*tOwnedInterpCells(i));
            moris_index tSubphase = tOwnedInterpCells(i)->get_subphase_index();

            if(!tCluster.is_trivial())
            {
                tNumNewInterpCellsOwned++;
                tNonTrivialOwnedInterpCells.push_back(tOwnedInterpCells(i));
                aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tOwnedInterpCells(i)->get_index();
                tCurrentNewInterpCellIndex++;
            }
            else
            {
                aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tOwnedInterpCells(i)->get_index();
            }
        }

        // allocate new interpolation cell ids
        moris_id tCurrentId = tEnrIpMesh.allocate_entity_ids(tNumNewInterpCellsOwned,EntityRank::ELEMENT);

        // assign new non trivial owned Interp cell Ids
        Cell<moris_id> tNewNonTrivialOwnedInterpCellsIds(tNumNewInterpCellsOwned);
        std::unordered_map<moris_id, moris_id> tBaseEnrIdToIndexInNonTrivialOwned;

        for(moris::size_t i = 0; i<tNonTrivialOwnedInterpCells.size(); i++)
        {
            tNewNonTrivialOwnedInterpCellsIds(i) = tCurrentId++;
            tBaseEnrIdToIndexInNonTrivialOwned[tNonTrivialOwnedInterpCells(i)->get_id()] = i;
        }

        Cell<Cell<Interpolation_Cell_Unzipped *>> tNonTrivialNotOwnedInterpCells(tNotOwnedInterpCells.size());

        // setup requests for not owned not trivial
        for(moris::size_t iP = 0; iP<tNotOwnedInterpCells.size(); iP++)
        {
            for(moris::size_t iC = 0; iC<tNotOwnedInterpCells(iP).size(); iC++)
            {

                xtk::Cell_Cluster const & tCluster = tEnrIgMesh.get_xtk_cell_cluster(*tNotOwnedInterpCells(iP)(iC));
                moris_index tSubphase = tNotOwnedInterpCells(iP)(iC)->get_subphase_index();

                if(!tCluster.is_trivial())
                {
                    tNumNewInterpCellsNotOwned++;
                    tNonTrivialNotOwnedInterpCells(iP).push_back(tNotOwnedInterpCells(iP)(iC));
                    aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tNotOwnedInterpCells(iP)(iC)->get_index();
                    tCurrentNewInterpCellIndex++;
                }
                else
                {
                    aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tNotOwnedInterpCells(iP)(iC)->get_index();
                }
            }
        }

        //        // prepare data to request
        //        Cell<Matrix<IndexMat>> tBaseEnrInterpCellId(tNonTrivialNotOwnedInterpCells.size());
        //
        //        for(moris::size_t iP = 0; iP<tNonTrivialNotOwnedInterpCells.size(); iP++)
        //        {
        //            tBaseEnrInterpCellId(iP).resize(1,tNonTrivialNotOwnedInterpCells(iP).size());
        //
        //            for(moris::size_t iC = 0; iC<tNonTrivialNotOwnedInterpCells(iP).size(); iC++)
        //            {
        //                tBaseEnrInterpCellId(iP)(iC) = tNonTrivialNotOwnedInterpCells(iP)(iC)->get_id();
        //            }
        //
        //            if(tNonTrivialNotOwnedInterpCells(iP).size() == 0)
        //            {
        //                tBaseEnrInterpCellId(iP).resize(1,1);
        //                tBaseEnrInterpCellId(iP)(0) = MORIS_INDEX_MAX;
        //            }
        //        }
        //
        //        // send requests
        //        moris::uint tMPITag = 301;
        //        for(moris::size_t iP = 0; iP<tProcRanks.size(); iP++)
        //        {
        //            mXTKModel->send_outward_requests(tMPITag, tProcRanks,tBaseEnrInterpCellId);
        //        }
        //
        //        barrier();
        //
        //        // receive requests
        //        Cell<Matrix<IndexMat>> tReceivedBaseEnrCellIds;
        //        Cell<uint> tProcsReceivedFrom;
        //        mXTKModel->inward_receive_requests(tMPITag, 1, tReceivedBaseEnrCellIds, tProcsReceivedFrom);
        //
        //        // prepare answers
        //        Cell<Matrix<IndexMat>> tEnrCellIds;
        //        this->prepare_ip_cell_id_answers(tReceivedBaseEnrCellIds,tNewNonTrivialOwnedInterpCellsIds,tEnrCellIds, tBaseEnrIdToIndexInNonTrivialOwned);
        //
        //        // return information
        //        mXTKModel->return_request_answers(tMPITag+1, tEnrCellIds, tProcsReceivedFrom);
        //
        //        barrier();
        //
        //        // receive the answers
        //        Cell<Matrix<IndexMat>> tReceivedEnrCellIds;
        //        mXTKModel->inward_receive_request_answers(tMPITag+1,1,tProcRanks,tReceivedEnrCellIds);
        //
        //        // allocate space in interpolation cells of enriched interpolation mesh
        //        tEnrIpMesh.mEnrichedInterpCells.resize(tNumNewInterpCellsNotOwned + tNumNewInterpCellsOwned + tEnrIpMesh.mEnrichedInterpCells.size());
        //
        //        // add interpolation cells for ghost to enriched interp mesh
        //        this->create_not_owned_ghost_ip_cells(aGhostSetupData,tEnrIpMesh,tNonTrivialNotOwnedInterpCells,tReceivedEnrCellIds);
        //
        //        // create owned ghost ip cells
        //        this->create_owned_ghost_ip_cells(aGhostSetupData,tEnrIpMesh,tNonTrivialOwnedInterpCells,tNewNonTrivialOwnedInterpCellsIds);
        //
        //        barrier();

    }
    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_ip_cell_id_answers( Cell<Matrix<IndexMat>> & aReceivedEnrCellIds,
            Cell<moris_id>                         & aNewInterpCellIds,
            Cell<Matrix<IndexMat>>                 & aEnrCellIds,
            std::unordered_map<moris_id, moris_id> & aBaseEnrIdToIndexInNonTrivialOwned)
    {
        aEnrCellIds.resize(aReceivedEnrCellIds.size());

        for(moris::uint iP = 0; iP< aReceivedEnrCellIds.size(); iP++)
        {
            aEnrCellIds(iP).resize(1,aReceivedEnrCellIds(iP).numel());

            if(aReceivedEnrCellIds(iP)(0,0) != MORIS_INDEX_MAX)
            {

                for(moris::uint iC = 0; iC< aReceivedEnrCellIds(iP).numel(); iC++)
                {
                    auto tIter = aBaseEnrIdToIndexInNonTrivialOwned.find(aReceivedEnrCellIds(iP)(iC));

                    if(tIter == aBaseEnrIdToIndexInNonTrivialOwned.end())
                    {
                        MORIS_ASSERT(tIter != aBaseEnrIdToIndexInNonTrivialOwned.end(),"Enriched cell id not in map");
                    }
                    aEnrCellIds(iP)(iC) = aNewInterpCellIds(tIter->second);

                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::create_not_owned_ghost_ip_cells(
            Ghost_Setup_Data                                & aGhostSetupData,
            Enriched_Interpolation_Mesh                     & aEnrInterpMesh,
            Cell<Cell<Interpolation_Cell_Unzipped *>> const & aNonTrivialNotOwnedInterpCells,
            Cell<Matrix<IndexMat>>                    const & aReceivedEnrCellIds)
    {
        //        // iterate through received data
        //        for(moris::uint i = 0; i < aNonTrivialNotOwnedInterpCells.size(); i++)
        //        {
        //            uint tNumReceivedReqs = aNonTrivialNotOwnedInterpCells(i).size();
        //
        //            // iterate through received requests
        //            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
        //            {
        //                moris_index tReceivedEnrCellId = aReceivedEnrCellIds(i)(j);
        //                moris_index tSubphase          = aNonTrivialNotOwnedInterpCells(i)(j)->get_subphase_index();
        //                moris_index tGhostCellIpIndex  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase);
        //
        //                MORIS_ASSERT(tGhostCellIpIndex != MORIS_INDEX_MAX,"Max index on not owned ghost ip cell. This could be a communication error.");
        //
        //                aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex) = new Interpolation_Cell_Unzipped(aNonTrivialNotOwnedInterpCells(i)(j)->get_base_cell(),
        //                        aNonTrivialNotOwnedInterpCells(i)(j)->get_subphase_index(),
        //                        aNonTrivialNotOwnedInterpCells(i)(j)->get_bulkphase_index(),
        //                        aNonTrivialNotOwnedInterpCells(i)(j)->get_id(),
        //                        tGhostCellIpIndex,
        //                        tReceivedEnrCellId,
        //                        aNonTrivialNotOwnedInterpCells(i)(j)->get_connectivity());
        //
        //                aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex)->set_vertices(aNonTrivialNotOwnedInterpCells(i)(j)->get_xtk_interpolation_vertices());
        //
        //            }
        //        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::create_owned_ghost_ip_cells(
            Ghost_Setup_Data                    & aGhostSetupData,
            Enriched_Interpolation_Mesh         & aEnrInterpMesh,
            Cell<Interpolation_Cell_Unzipped *> & aNonTrivialOwnedInterpCells,
            Cell<moris_id>                      & aEnrCellIds)
    {
        //        // iterate through received data
        //        for(moris::uint i = 0; i < aNonTrivialOwnedInterpCells.size(); i++)
        //        {
        //            moris_index tEnrCellId = aEnrCellIds(i);
        //            moris_index tSubphase          = aNonTrivialOwnedInterpCells(i)->get_subphase_index();
        //            moris_index tGhostCellIpIndex  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase);
        //
        //            aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex) = new Interpolation_Cell_Unzipped(aNonTrivialOwnedInterpCells(i)->get_base_cell(),
        //                    aNonTrivialOwnedInterpCells(i)->get_subphase_index(),
        //                    aNonTrivialOwnedInterpCells(i)->get_bulkphase_index(),
        //                    tEnrCellId,
        //                    tGhostCellIpIndex,
        //                    aNonTrivialOwnedInterpCells(i)->get_owner(),
        //                    aNonTrivialOwnedInterpCells(i)->get_connectivity());
        //
        //            aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex)->set_vertices(aNonTrivialOwnedInterpCells(i)->get_xtk_interpolation_vertices());
        //        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::identify_and_setup_aura_vertices_in_ghost()
    {
        // access the enriched interpolation mesh
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // access the enriched interpolation cells
        Cell<Interpolation_Cell_Unzipped*> & tEnrIpCells = tEnrInterpMesh.get_enriched_interpolation_cells();

        // access subphase neighborhood information
        moris::Cell<moris::Cell<moris_index>>  const & tSubphaseToSubphase = mXTKModel->get_subphase_to_subphase();

        // number of subphases in mesh
        moris::uint tNumSubphases = tSubphaseToSubphase.size();

        // par rank
        moris_index tParRank = moris::par_rank();

        // list of enriched interp
        std::unordered_map<moris_index,bool> tIpCellsWithVerticesToComm;

        // iterate through subphases and collect all not owned subphases and neighbors that may have ghost
        for(moris::uint i = 0; i < tNumSubphases; i++)
        {
            // get my enriched Ip cell (ghost one)
//            Interpolation_Cell_Unzipped* tMyEnrichedIpCell = tEnrIpCells((moris_index)i);

            // only continue into the neighbors if I do not own this cell/subphase
//            if(tParRank != tMyEnrichedIpCell->get_owner())
//            {
                // iterate through subphase neighbors
                for(moris::uint j = 0; j < tSubphaseToSubphase(i).size(); j++)
                {
                    // neighbor index
                    moris_index tNeighborSubphaseIndex = tSubphaseToSubphase(i)(j);

                    // add it to ip cells to consider if i am supposed to construct ghost at this neighbor relationship
                    if(tParRank != tEnrIpCells(tNeighborSubphaseIndex)->get_owner())
                    {
                        if(tIpCellsWithVerticesToComm.find(tNeighborSubphaseIndex) == tIpCellsWithVerticesToComm.end() )
                        {
                            tIpCellsWithVerticesToComm[tNeighborSubphaseIndex] = true;
                        }
                    }
                }
//            }
        }

        // iterate through interpolation cells that are in the aura and collect vertex interpolation that are not owned
        // start setting up data for communication
        // prepare outward requests
        Cell<Matrix<IndexMat>> tNotOwnedIPVertIndsToProcs;
        Cell<Matrix<IndexMat>> tNotOwnedBGIPVertsIdsToProcs;
        Cell<Matrix<IndexMat>> tNotOwnedEnrichedCellIdToProcs;
        Cell<Matrix<IndexMat>> tNotOwnedEnrichedCellBulkPhaseToProcs; // for checking against
        Cell<uint>           tProcRanks;
        std::unordered_map<moris_id,moris_id>  tProcRankToDataIndex;
        this->prepare_interpolation_vertex_t_matrix_requests(tIpCellsWithVerticesToComm, tNotOwnedIPVertIndsToProcs,
                tNotOwnedBGIPVertsIdsToProcs, tNotOwnedEnrichedCellIdToProcs, tNotOwnedEnrichedCellBulkPhaseToProcs,
                tProcRanks,tProcRankToDataIndex);




        // send requests
        moris::uint tMPITag = 3001;
        // send the background vertex id
        mXTKModel->send_outward_requests(tMPITag, tProcRanks,tNotOwnedBGIPVertsIdsToProcs);

        // send the enriched interpolation cell id
        mXTKModel->send_outward_requests(tMPITag+1, tProcRanks,tNotOwnedEnrichedCellIdToProcs);

        // send the enriched interpolation cell bulk phase ids
        mXTKModel->send_outward_requests(tMPITag+2, tProcRanks,tNotOwnedEnrichedCellBulkPhaseToProcs);


        barrier();

        // receive requests for the t-matrices
        Cell<Matrix<IndexMat>> tReceivedVertexIds;
        Cell<Matrix<IndexMat>> tReceivedEnrichedCellId;
        Cell<Matrix<IndexMat>> tReceivedEnrichedCellBulkPhase;
        Cell<uint> tProcsReceivedFrom1;
        Cell<uint> tProcsReceivedFrom2;
        Cell<uint> tProcsReceivedFrom3;
        mXTKModel->inward_receive_requests(tMPITag, 1, tReceivedVertexIds, tProcsReceivedFrom1); // receive the requests ofr BG VertexIds
        mXTKModel->inward_receive_requests(tMPITag+1, 1, tReceivedEnrichedCellId, tProcsReceivedFrom2);// recieve the requests for Enriched IP Cell Id
        mXTKModel->inward_receive_requests(tMPITag+2, 1, tReceivedEnrichedCellBulkPhase, tProcsReceivedFrom3);// recieve the requests for Enriched IP Cell Bulk phases

        // prepare the t-matrices for sending
        Cell<Matrix<DDRMat>>   tTMatrixWeights;
        Cell<Matrix<IndexMat>> tTMatrixIndices;
        Cell<Matrix<IndexMat>> tTMatrixOwners;
        Cell<Matrix<IndexMat>> tTMatrixOffsets;
        this->prepare_t_matrix_request_answers(tReceivedVertexIds,tReceivedEnrichedCellId,tReceivedEnrichedCellBulkPhase,
                tTMatrixWeights,tTMatrixIndices,tTMatrixOwners,tTMatrixOffsets);


        // send information
        mXTKModel->return_request_answers_reals(tMPITag+3, tTMatrixWeights, tProcsReceivedFrom1);
        mXTKModel->return_request_answers(tMPITag+4, tTMatrixIndices, tProcsReceivedFrom1);
        mXTKModel->return_request_answers(tMPITag+5, tTMatrixOwners, tProcsReceivedFrom1);
        mXTKModel->return_request_answers(tMPITag+6, tTMatrixOffsets, tProcsReceivedFrom1);

        // wait
        barrier();

        // receive
        Cell<Matrix<DDRMat>>   tRequestedTMatrixWeights;
        Cell<Matrix<IndexMat>> tRequestedTMatrixIndices;
        Cell<Matrix<IndexMat>> tRequestedTMatrixOwners;
        Cell<Matrix<IndexMat>> tRequestedTMatrixOffsets;

        // receive the answers
        mXTKModel->inward_receive_request_answers_reals(tMPITag+3,1,tProcRanks,tRequestedTMatrixWeights);
        mXTKModel->inward_receive_request_answers(tMPITag+4,1,tProcRanks,tRequestedTMatrixIndices);
        mXTKModel->inward_receive_request_answers(tMPITag+5,1,tProcRanks,tRequestedTMatrixOwners);
        mXTKModel->inward_receive_request_answers(tMPITag+6,1,tProcRanks,tRequestedTMatrixOffsets);

        barrier();

        // commit it to my data
        this->handle_received_interpolation_data(tNotOwnedIPVertIndsToProcs,tNotOwnedEnrichedCellBulkPhaseToProcs,
                tRequestedTMatrixWeights,tRequestedTMatrixIndices,tRequestedTMatrixOwners,tRequestedTMatrixOffsets);

        //wait
        barrier();
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_interpolation_vertex_t_matrix_requests(
            std::unordered_map<moris_index,bool> const & aNotOwnedIpCellsInGhost,
            Cell<Matrix<IndexMat>>                     & aNotOwnedIPVertIndsToProcs,
            Cell<Matrix<IndexMat>>                     & aNotOwnedBGIPVertsIdsToProcs,
            Cell<Matrix<IndexMat>>                     & aNotOwnedIpCellIdToProcs,
            Cell<Matrix<IndexMat>>                     & aNotOwnedEnrichedCellBulkPhaseToProcs,
            Cell<uint>                                 & aProcRanks,
            std::unordered_map<moris_id,moris_id>      & aProcRankToDataIndex)
    {
        // access the enriched interpolation mesh
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        Cell<Interpolation_Cell_Unzipped*> & tEnrIpCells = tEnrInterpMesh.get_enriched_interpolation_cells();

        // Counter and current proc index
        Cell<moris_id> tCounts(0);

        // temporary cell of cells which will be converted to matrices later
        Cell<Cell<moris_id>> tNotOwnedIPVertIndsToProcs;
        Cell<Cell<moris_id>> tNotOwnedBGIPVertsIdsToProcs;
        Cell<Cell<moris_id>> tNotOwnedIpCellIdToProcs;
        Cell<Cell<moris_id>> tNotOwnedIpCellBulkPhase;

        // get the communication table
        Matrix<IndexMat> tCommTable  = tEnrInterpMesh.get_communication_table();

        // size
        aProcRanks.resize(tCommTable.numel());

        for(moris::uint i = 0; i <tCommTable.numel(); i++)
        {
            aProcRankToDataIndex[tCommTable(i)] = i;
            aProcRanks(i) = (tCommTable(i));

            // add cell for verts
            tNotOwnedIPVertIndsToProcs.push_back(Cell<moris_id>(0));
            tNotOwnedBGIPVertsIdsToProcs.push_back(Cell<moris_id>(0));
            tNotOwnedIpCellIdToProcs.push_back(Cell<moris_id>(0));
            tNotOwnedIpCellBulkPhase.push_back(Cell<moris_id>(0));
        }

        for (auto tIpCell : aNotOwnedIpCellsInGhost)
        {
            // access the ip cell index
            moris_index tIpCellIndex = tIpCell.first;

            // get the ip cell
            Interpolation_Cell_Unzipped* tEnrIpCell = tEnrIpCells(tIpCellIndex);

            // get the vertices
            moris::Cell< xtk::Interpolation_Vertex_Unzipped* > & tVertexPointers = tEnrIpCell->get_xtk_interpolation_vertices();

            // iterate through vertices check if the owner is there
            for(moris::uint iV = 0; iV<tVertexPointers.size(); iV++)
            {
                // vertex owner
                moris_index tOwnerProc = tVertexPointers(iV)->get_owner();

                // If I don't own this vertex, and its on a cell that has a vertex
                if(!tVertexPointers(iV)->get_base_vertex()->has_interpolation(1))
                {
                    // get the index of this proc
                    auto tProcIndexInData = aProcRankToDataIndex.find(tOwnerProc);

                    tNotOwnedIPVertIndsToProcs(tProcIndexInData->second).push_back(tVertexPointers(iV)->get_index());
                    tNotOwnedBGIPVertsIdsToProcs(tProcIndexInData->second).push_back(tVertexPointers(iV)->get_base_vertex()->get_id());
                    tNotOwnedIpCellIdToProcs(tProcIndexInData->second).push_back(tEnrIpCell->get_id());
                    tNotOwnedIpCellBulkPhase(tProcIndexInData->second).push_back(tEnrIpCell->get_bulkphase_index());
                }
            }
        }

        // populate matrix in input data
        aNotOwnedIPVertIndsToProcs.clear();
        aNotOwnedIPVertIndsToProcs.resize(tNotOwnedIPVertIndsToProcs.size());
        aNotOwnedBGIPVertsIdsToProcs.resize(tNotOwnedBGIPVertsIdsToProcs.size());
        aNotOwnedIpCellIdToProcs.resize(tNotOwnedIpCellBulkPhase.size());
        aNotOwnedEnrichedCellBulkPhaseToProcs.resize(tNotOwnedIpCellBulkPhase.size());
        for(moris::uint iD = 0; iD< tNotOwnedIPVertIndsToProcs.size(); iD++)
        {
            aNotOwnedIPVertIndsToProcs(iD).resize(1,tNotOwnedIPVertIndsToProcs(iD).size());
            aNotOwnedBGIPVertsIdsToProcs(iD).resize(1,tNotOwnedBGIPVertsIdsToProcs(iD).size());
            aNotOwnedIpCellIdToProcs(iD).resize(1,tNotOwnedIpCellIdToProcs(iD).size());
            aNotOwnedEnrichedCellBulkPhaseToProcs(iD).resize(1,tNotOwnedIpCellBulkPhase(iD).size());
            for(moris::uint jD = 0; jD< tNotOwnedIPVertIndsToProcs(iD).size(); jD++)
            {
                aNotOwnedIPVertIndsToProcs(iD)(jD) = tNotOwnedIPVertIndsToProcs(iD)(jD);
                aNotOwnedBGIPVertsIdsToProcs(iD)(jD) = tNotOwnedBGIPVertsIdsToProcs(iD)(jD);
                aNotOwnedIpCellIdToProcs(iD)(jD) = tNotOwnedIpCellIdToProcs(iD)(jD);
                aNotOwnedEnrichedCellBulkPhaseToProcs(iD)(jD) = tNotOwnedIpCellBulkPhase(iD)(jD);
            }

            if(tNotOwnedIPVertIndsToProcs(iD).size() == 0)
            {
                aNotOwnedIPVertIndsToProcs(iD).resize(1,1);
                aNotOwnedBGIPVertsIdsToProcs(iD).resize(1,1);
                aNotOwnedIpCellIdToProcs(iD).resize(1,1);
                aNotOwnedEnrichedCellBulkPhaseToProcs(iD).resize(1,1);
                aNotOwnedIPVertIndsToProcs(iD)(0) =MORIS_INDEX_MAX;
                aNotOwnedBGIPVertsIdsToProcs(iD)(0) = MORIS_INDEX_MAX;
                aNotOwnedIpCellIdToProcs(iD)(0) = MORIS_INDEX_MAX;
                aNotOwnedEnrichedCellBulkPhaseToProcs(iD)(0) = MORIS_INDEX_MAX;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::prepare_t_matrix_request_answers(
            Cell<Matrix<IndexMat>> const & aRequestedBgVertexIds,
            Cell<Matrix<IndexMat>> const & aRequestedIpCellIds,
            Cell<Matrix<IndexMat>> const & aIpCellBulkPhases,
            Cell<Matrix<DDRMat>>         & aTMatrixWeights,
            Cell<Matrix<IndexMat>>       &  aTMatrixIndices,
            Cell<Matrix<IndexMat>>       &  aBasisOwners,
            Cell<Matrix<IndexMat>>       &  aTMatrixOffsets)
    {
        // access enriched ip mesh
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        MORIS_ASSERT(aRequestedBgVertexIds.size() == aRequestedIpCellIds.size(),"Mismatch in received communication information.");

        // information about size of interpolation mats
        Cell<uint> tSizes(aRequestedBgVertexIds.size(),0);

        //resize the input data
        aTMatrixWeights.resize(aRequestedBgVertexIds.size());
        aTMatrixIndices.resize(aRequestedBgVertexIds.size());
        aBasisOwners.resize(aRequestedBgVertexIds.size());
        aTMatrixOffsets.resize(aRequestedBgVertexIds.size());

        // collect the vertex interpolations
        Cell<Cell<Vertex_Enrichment* >> tVertexInterpolations(aRequestedBgVertexIds.size());

        // collect size information throughout loop
        Cell<moris_index> tDataSizes(aRequestedBgVertexIds.size(),0);

        // iterate through and figure out how big to make the weights and indices mats
        // also collect vertex interpolations
        for(moris::uint iP = 0; iP < aRequestedBgVertexIds.size(); iP++)
        {
            // size the tmatrix offset for each vertex requested (num verts +1)
            aTMatrixOffsets(iP).resize(1,aRequestedBgVertexIds(iP).numel()+1 );

            // set the first one to 0
            aTMatrixOffsets(iP)(0) = 0;

            // iterate through the vertices and get their interpolations and figure out
            // how big it is
            for(moris::uint iV = 0; iV < aRequestedBgVertexIds(iP).numel(); iV++)
            {

                // check that the bulk phases are consistent
                moris_index tCellIndex = tEnrInterpMesh.get_loc_entity_ind_from_entity_glb_id( aRequestedIpCellIds(iP)(iV), EntityRank::ELEMENT, 0);

                // get the cell
                Interpolation_Cell_Unzipped* tEnrIpCell = tEnrInterpMesh.get_enriched_interpolation_cells()(tCellIndex);

                // verifythat the bulk phases are consistent across procs
                MORIS_ERROR(tEnrIpCell->get_bulkphase_index() == aIpCellBulkPhases(iP)(iV),"Parallel bulkphase mismatch.");


                // get the vertex
                moris_index tVertexIndex = this->get_enriched_interpolation_vertex(aRequestedBgVertexIds(iP)(iV), aRequestedIpCellIds(iP)(iV));

                // get the vertex interpolation
                Vertex_Enrichment* tVertexInterp = tEnrInterpMesh.get_vertex_enrichment(0,tVertexIndex);

                MORIS_ASSERT(tVertexInterp->get_base_vertex_interpolation() != nullptr,"Owning proc has a nullptr for the vertex interpolation.");

                tVertexInterpolations(iP).push_back(tVertexInterp);

                // number of basis functions interpolating into this vertex
                moris_index tNumBasis = tVertexInterp->get_basis_indices().numel();

                // offsets
                aTMatrixOffsets(iP)(iV+1) = aTMatrixOffsets(iP)(iV) + tNumBasis;

                // add to size
                tDataSizes(iP) = tDataSizes(iP) + tNumBasis;
            }
        }

        //  iterate through and size data
        for(moris::uint iP = 0; iP < aRequestedBgVertexIds.size(); iP++)
        {
            aTMatrixWeights(iP).resize(1,tDataSizes(iP));
            aTMatrixIndices(iP).resize(1,tDataSizes(iP));
            aBasisOwners(iP).resize(1,tDataSizes(iP));
        }

        // populate the data
        for(moris::uint iP = 0; iP < aRequestedBgVertexIds.size(); iP++)
        {
            moris::uint tCount = 0;

            for(moris::uint iV = 0; iV < tVertexInterpolations(iP).size(); iV++)
            {
                this->add_vertex_interpolation_to_communication_data(
                        tCount,
                        tVertexInterpolations(iP)(iV),
                        aTMatrixWeights(iP),
                        aTMatrixIndices(iP),
                        aBasisOwners(iP),
                        aTMatrixOffsets(iP));
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::handle_received_interpolation_data(
            Cell<Matrix<IndexMat>> const & aNotOwnedIPVertIndsToProcs,
            Cell<Matrix<IndexMat>> const & aNotOwnedEnrichedCellBulkPhaseToProcs,
            Cell<Matrix<DDRMat>>   const & aRequestedTMatrixWeights,
            Cell<Matrix<IndexMat>> const & aRequestedTMatrixIndices,
            Cell<Matrix<IndexMat>> const & aRequestedBasisOwners,
            Cell<Matrix<IndexMat>> const & aRequestedTMatrixOffsets)
    {
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // access the communication
        Matrix<IdMat> tCommTable = tEnrInterpMesh.get_communication_table();

        std::unordered_map<moris_id,moris_id> tProcRankToIndexInData;

        moris::uint tCount = tCommTable.numel();

        // resize proc ranks and setup map to comm table
        for(moris::uint i = 0; i <tCommTable.numel(); i++)
        {
            tProcRankToIndexInData[tCommTable(i)] = i;
        }

        // iterate through returned information
        for(moris::uint iP = 0; iP < aNotOwnedIPVertIndsToProcs.size(); iP++)
        {

            // extract the t-matrices and basis ids/owners for the proc ip
            Cell<Matrix<DDRMat>>   tExtractedTMatrixWeights;
            Cell<Matrix<IndexMat>> tExtractedTMatrixIds;
            Cell<Matrix<IndexMat>> tExtractedTBasisOwners;

            this->extract_vertex_interpolation_from_communication_data(
                    aNotOwnedIPVertIndsToProcs(iP).numel(),
                    aRequestedTMatrixWeights(iP),
                    aRequestedTMatrixIndices(iP),
                    aRequestedBasisOwners(iP),
                    aRequestedTMatrixOffsets(iP),
                    tExtractedTMatrixWeights,
                    tExtractedTMatrixIds,
                    tExtractedTBasisOwners);

            // verify consistent sizes
            MORIS_ASSERT(aNotOwnedIPVertIndsToProcs(iP).numel() == tExtractedTMatrixWeights.size(),"Size mismatch in t-matrix weights.");
            MORIS_ASSERT(aNotOwnedIPVertIndsToProcs(iP).numel() == tExtractedTMatrixIds.size(),"Size mismatch in t-matrix ids.");
            MORIS_ASSERT(aNotOwnedIPVertIndsToProcs(iP).numel() == tExtractedTBasisOwners.size(),"Size mismatch in basis owners.");
            MORIS_ASSERT(aNotOwnedIPVertIndsToProcs(iP).numel() == aNotOwnedEnrichedCellBulkPhaseToProcs(iP).numel(),"Size mismatch in bulk phases.");

            // iterate through vertices and set their interpolation weights and basis ids
            for(moris::uint iV = 0; iV < aNotOwnedIPVertIndsToProcs(iP).numel(); iV++)
            {
                // get the vertex
                moris_index tVertexIndex = aNotOwnedIPVertIndsToProcs(iP)(iV);

                Interpolation_Vertex_Unzipped & tVertex = tEnrInterpMesh.get_xtk_interp_vertex(tVertexIndex);

                // get the enriched vertex interpolation
                Vertex_Enrichment * tVertexInterp = tVertex.get_xtk_interpolation( tEnrInterpMesh.mMeshIndices(0) );

                // iterate through basis functions and find local indices
                moris::Matrix<IndexMat> tBasisIndices(tExtractedTMatrixIds(iV).n_rows(),tExtractedTMatrixIds(iV).n_cols());

                for(moris::uint iBs = 0; iBs < tExtractedTMatrixIds(iV).numel(); iBs++)
                {
                    // basis id
                    moris_id tId = tExtractedTMatrixIds(iV)(iBs);

                    // add this basis to the mesh if it doesnt exists on the current partition
                    if(!tEnrInterpMesh.basis_exists_on_partition(tEnrInterpMesh.mMeshIndices(0),tId))
                    {
                        MORIS_ASSERT(tExtractedTBasisOwners(iV)(iBs) != par_rank(),"Owned basis should already exist on partition.");

                        tEnrInterpMesh.add_basis_function(tEnrInterpMesh.mMeshIndices(0),tId,
                                tExtractedTBasisOwners(iV)(iBs),
                                aNotOwnedEnrichedCellBulkPhaseToProcs(iP)(iV));
                    }

                    tBasisIndices(iBs) = tEnrInterpMesh.get_enr_basis_index_from_enr_basis_id(tEnrInterpMesh.mMeshIndices(0),tId);

                    moris_id tBasisOwner = tExtractedTBasisOwners(iV)(iBs);

                    MORIS_ASSERT(tEnrInterpMesh.get_basis_owner(tBasisIndices(iBs),tEnrInterpMesh.mMeshIndices(0)) == tBasisOwner,"Ownership discrepency.");
                    MORIS_ASSERT(tEnrInterpMesh.get_basis_bulk_phase(tBasisIndices(iBs),tEnrInterpMesh.mMeshIndices(0)) == aNotOwnedEnrichedCellBulkPhaseToProcs(iP)(iV),"Bulkphase discrepency.");

                    // if the basis has an owning proc that is not in the comm table, add it to the comm table
                    if(tProcRankToIndexInData.find(tBasisOwner) == tProcRankToIndexInData.end() && tBasisOwner != par_rank())
                    {
                        tEnrInterpMesh.add_proc_to_comm_table(tBasisOwner);
                        tProcRankToIndexInData[tBasisOwner] = tCount;
                        tCount++;
                    }
                }

                // iterate through basis in the base vertex interpolation
                moris::uint tNumCoeffs = tExtractedTMatrixIds(iV).numel();

                // Setup the map in the basis function
                std::unordered_map<moris::moris_index,moris::moris_index> & tVertEnrichMap = tVertexInterp->get_basis_map();

                for(moris::uint iB = 0; iB<tNumCoeffs; iB++)
                {
                    moris::moris_index tBasisIndex = tBasisIndices(iB);

                    tVertEnrichMap[tBasisIndex] = iB;
                }

                // get the basis indices from the basis ids
                tVertexInterp->add_basis_information(tBasisIndices,tExtractedTMatrixIds(iV));
                tVertexInterp->add_basis_weights(tBasisIndices,tExtractedTMatrixWeights(iV));
                tVertexInterp->add_basis_owners(tBasisIndices,tExtractedTBasisOwners(iV));
                tVertexInterp->add_base_vertex_interpolation(nullptr); // base vertex interpolation does not exists (other  proc)
            }
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Ghost_Stabilization::get_enriched_interpolation_vertex(
            moris_index const & aBGVertId,
            moris_index const & aEnrichedIpCellId)
    {
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // get the interpolation cell index using the id
        moris_index tCellIndex = tEnrInterpMesh.get_loc_entity_ind_from_entity_glb_id(aEnrichedIpCellId,EntityRank::ELEMENT,0);

        // get the cell
        Interpolation_Cell_Unzipped* tEnrIpCell = tEnrInterpMesh.get_enriched_interpolation_cells()(tCellIndex);

        // get the vertices
        moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tVertexPointers = tEnrIpCell->get_xtk_interpolation_vertices();

        moris_index tVertexPointerInd = 0;
        uint tCount=0;

        for(moris::uint i  = 0; i < tVertexPointers.size(); i ++ )
        {
            if(tVertexPointers(i)->get_base_vertex()->get_id() == aBGVertId)
            {
                tVertexPointerInd = tVertexPointers(i)->get_index();
                tCount++;
            }
        }

        MORIS_ERROR(tCount==1,"Enriched interpolation vertex not found or found more than once");

        return tVertexPointerInd;
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::declare_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData)
    {
        uint tNumBulkPhases = mXTKModel->get_geom_engine()->get_num_bulk_phase();

        Cell<std::string> tGhostDoubleSideNames(tNumBulkPhases);

        for(moris::moris_index iP0 = 0; iP0 <(moris_index) tNumBulkPhases; iP0++)
        {

            std::string tGhostSideSetName =this->get_ghost_dbl_side_set_name(iP0);

            tGhostDoubleSideNames(iP0) = tGhostSideSetName;

        }
        aGhostSetupData.mDblSideSetIndexInMesh = mXTKModel->get_enriched_integ_mesh(0).register_double_side_set_names(tGhostDoubleSideNames);
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::add_vertex_interpolation_to_communication_data(
            moris::uint      & aCount,
            Vertex_Enrichment* aInterpolation,
            Matrix<DDRMat>   & aTMatrixWeights,
            Matrix<IndexMat> & aTMatrixIndices,
            Matrix<IndexMat> & aTMatrixOwners,
            Matrix<IndexMat> & aTMatrixOffsets)
    {
        // access the basis indices and weights
        moris::Matrix< moris::IndexMat > const & tBasisIndices = aInterpolation->get_basis_ids();
        moris::Matrix< moris::DDRMat >   const & tBasisWeights = aInterpolation->get_basis_weights();
        moris::Matrix< moris::IndexMat >         tBasisOwners  = aInterpolation->get_owners();

        for(moris::uint i = 0 ; i < tBasisIndices.numel(); i ++ )
        {
            aTMatrixIndices(aCount) = tBasisIndices(i);
            aTMatrixWeights(aCount) = tBasisWeights(i);

            MORIS_ASSERT(tBasisOwners(i) < par_size() || tBasisOwners(i) > 0, "Bad ownership for basis function.");

            aTMatrixOwners(aCount)  = tBasisOwners(i);
            aCount++;
        }

    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::extract_vertex_interpolation_from_communication_data(
            moris::uint      const & aNumVerts,
            Matrix<DDRMat>   const & aTMatrixWeights,
            Matrix<IndexMat> const & aTMatrixIndices,
            Matrix<IndexMat> const & aTMatrixOwners,
            Matrix<IndexMat> const & aTMatrixOffsets,
            Cell<Matrix<DDRMat>>   & aExtractedTMatrixWeights,
            Cell<Matrix<IndexMat>> & aExtractedTMatrixIndices,
            Cell<Matrix<IndexMat>> & aExtractedBasisOwners)
    {

        // size output data
        aExtractedTMatrixWeights.resize(aNumVerts);
        aExtractedTMatrixIndices.resize(aNumVerts);
        aExtractedBasisOwners.resize(aNumVerts);

        // current starting index
        moris_index tStart = 0;

        // extract the data into the cells
        for(moris::uint iV =0; iV < aNumVerts; iV++)
        {
            // number of basis interpolating into the vertex
            moris::moris_index tNumBasis = aTMatrixOffsets(iV+1) - tStart;

            aExtractedTMatrixWeights(iV).resize(tNumBasis,1);
            aExtractedTMatrixIndices(iV).resize(tNumBasis,1);
            aExtractedBasisOwners(iV).resize(1,tNumBasis);

            // itere and grab  data
            for(moris::moris_index iIp = 0; iIp < tNumBasis; iIp++ )
            {
                aExtractedTMatrixWeights(iV)(iIp) = aTMatrixWeights(tStart+iIp);
                aExtractedTMatrixIndices(iV)(iIp) = aTMatrixIndices(tStart+iIp);
                aExtractedBasisOwners(iV)(iIp)    = aTMatrixOwners(tStart+iIp);
            }

            tStart = aTMatrixOffsets(iV+1);
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::construct_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData)
    {
        // enriched interpolation mesh
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();
        Enriched_Integration_Mesh   & tEnrIntegMesh = mXTKModel->get_enriched_integ_mesh();

        // all interpolation cells
        Cell<Interpolation_Cell_Unzipped*> & tEnrIpCells = tEnrInterpMesh.get_enriched_interpolation_cells();

        // access subphase neighborhood information
        moris::Cell<moris::Cell<moris_index>>  const & tSubphaseToSubphase                 = mXTKModel->get_subphase_to_subphase();
        moris::Cell<moris::Cell<moris_index>>  const & tSubphaseToSubphaseMySideOrds       = mXTKModel->get_subphase_to_subphase_my_side_ords();
        moris::Cell<moris::Cell<moris_index>>  const & tSubphaseToSubphaseNeighborSideOrds = mXTKModel->get_subphase_to_subphase_neighbor_side_ords();
        moris::Cell<moris::Cell<moris_index>>  const & tTransitionLocation                 = mXTKModel->get_subphase_to_subphase_transition_loc();

        // number of bulk phases in the mesh
        moris::uint tNumBulkPhases = mXTKModel->get_geom_engine()->get_num_bulk_phase();

        // number of subphases in mesh
        moris::uint tNumSubphases = tSubphaseToSubphase.size();

        // reserve space in data
        aGhostSetupData.mMasterSideIpCells.reserve(tNumBulkPhases*tNumSubphases);
        aGhostSetupData.mSlaveSideIpCells.reserve(tNumBulkPhases*tNumSubphases);
        aGhostSetupData.mMasterSideIgCellSideOrds.reserve(tNumBulkPhases*tNumSubphases);
        aGhostSetupData.mSlaveSideIgCellSideOrds.reserve(tNumBulkPhases*tNumSubphases);
        aGhostSetupData.mTrivialFlag.reserve(tNumBulkPhases*tNumSubphases);
        aGhostSetupData.mTransitionLocation.reserve(tNumBulkPhases*tNumSubphases);

        aGhostSetupData.mMasterSideIpCells.resize(tNumBulkPhases);
        aGhostSetupData.mSlaveSideIpCells.resize(tNumBulkPhases);
        aGhostSetupData.mMasterSideIgCellSideOrds.resize(tNumBulkPhases);
        aGhostSetupData.mSlaveSideIgCellSideOrds.resize(tNumBulkPhases);
        aGhostSetupData.mTrivialFlag.resize(tNumBulkPhases);
        aGhostSetupData.mTransitionLocation.resize(tNumBulkPhases);

        //        mXTKModel->print_subphase_neighborhood();

        // flag indicating whether its trivial or not
        moris_index tTrivial = 0;
        moris_index tNonTrivialCount = 0;

        // iterate through subphases
        for(moris::uint i = 0; i < tNumSubphases; i++)
        {
            // iterate through subphase neighbors
            for(moris::uint j = 0; j < tSubphaseToSubphase(i).size(); j++)
            {
                // if I am the one constructing this subphase then add it to ghost setup data
                if(this->create_ghost(aGhostSetupData,(moris_index)i,tSubphaseToSubphase(i)(j),tTrivial))
                {
                    //                    std::cout<<"Create ghost between" << " SP "<<std::setw(6)<<i<<" and SP " <<std::setw(6)<< tSubphaseToSubphase(i)(j)<<" Trivial? "<< tTrivial <<std::endl;
                    moris_index tFirstInterpCellIndex  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(i);
                    moris_index tSecondInterpCellIndex = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphaseToSubphase(i)(j));

                    Interpolation_Cell_Unzipped* tFirstInterpCell = tEnrIpCells(tFirstInterpCellIndex);
                    Interpolation_Cell_Unzipped* tSecondInterpCell = tEnrIpCells(tSecondInterpCellIndex);

                    //get the bulk phase
                    moris_index tBulkPhase = tFirstInterpCell->get_bulkphase_index();
                    MORIS_ASSERT(tBulkPhase == tSecondInterpCell->get_bulkphase_index(),"Bulk phase between neighboring subphases does not match");

                    // setup ip cell indices in ghost setup data
                    aGhostSetupData.mMasterSideIpCells(tBulkPhase).push_back(tFirstInterpCell->get_index());
                    aGhostSetupData.mSlaveSideIpCells(tBulkPhase).push_back(tSecondInterpCell->get_index());

                    // setup ig cells in ghost set up data
                    aGhostSetupData.mMasterSideIgCellSideOrds(tBulkPhase).push_back(tSubphaseToSubphaseMySideOrds(i)(j));
                    aGhostSetupData.mSlaveSideIgCellSideOrds(tBulkPhase).push_back(tSubphaseToSubphaseNeighborSideOrds(i)(j));

                    // mark as trivial or not in ghost setup data
                    aGhostSetupData.mTrivialFlag(tBulkPhase).push_back(tTrivial);

                    // mark the transition location
                    aGhostSetupData.mTransitionLocation(tBulkPhase).push_back(tTransitionLocation(i)(j));

                    // increment count
                    if(tTrivial > 0)
                    {
                        tNonTrivialCount++;
                    }
                }
            }
        }

        aGhostSetupData.mMasterSideIpCells.shrink_to_fit();
        aGhostSetupData.mSlaveSideIpCells.shrink_to_fit();
        aGhostSetupData.mMasterSideIgCellSideOrds.shrink_to_fit();
        aGhostSetupData.mSlaveSideIgCellSideOrds.shrink_to_fit();
        aGhostSetupData.mTrivialFlag.shrink_to_fit();
        aGhostSetupData.mTransitionLocation.shrink_to_fit();

        // allocate ids for non-trivial integration cells
        moris_id tCurrentId    = tEnrIntegMesh.allocate_entity_ids(tNonTrivialCount,EntityRank::ELEMENT);
        moris_id tCurrentIndex = tEnrIntegMesh.get_num_entities(EntityRank::ELEMENT);

        // iterate through bulk phases
        for(moris::uint i = 0; i < aGhostSetupData.mMasterSideIpCells.size(); i++)
        {
            // allocate space in the integration mesh double side sets
            tEnrIntegMesh.mDoubleSideSets(aGhostSetupData.mDblSideSetIndexInMesh(i)).resize(aGhostSetupData.mMasterSideIpCells(i).size());

            // iterate through double sides in this bulk phase
            for(moris::uint j = 0; j < aGhostSetupData.mMasterSideIpCells(i).size(); j++)
            {
                // create a new side cluster for each of the pairs
                std::shared_ptr<Side_Cluster> tSlaveSideCluster  = this->create_slave_side_cluster(aGhostSetupData,tEnrIpCells,i,j);
                std::shared_ptr<Side_Cluster>tMasterSideCluster = this->create_master_side_cluster(aGhostSetupData,tEnrIpCells,i,j,tSlaveSideCluster.get(),tCurrentIndex,tCurrentId);

                // verify the subphase cluster
                MORIS_ASSERT(tSlaveSideCluster->mInterpolationCell->get_bulkphase_index() == (moris_index)i,"Bulk phase mismatch on slave side of double side set cluster");
                MORIS_ASSERT(tMasterSideCluster->mInterpolationCell->get_bulkphase_index() == (moris_index)i,"Bulk phase mismatch on master side of double side set cluster");

                MORIS_ASSERT(tEnrInterpMesh.verify_basis_interpolating_into_cluster(*tSlaveSideCluster.get(),0),"Basis issue in cluster");
                MORIS_ASSERT(tEnrInterpMesh.verify_basis_interpolating_into_cluster(*tMasterSideCluster.get(),0),"Basis issue in cluster");

                // add to side clusters the integration mesh
                tEnrIntegMesh.mDoubleSideSetsMasterIndex(aGhostSetupData.mDblSideSetIndexInMesh(i)).push_back(tEnrIntegMesh.mDoubleSideSingleSideClusters.size());
                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back(tMasterSideCluster);
                tEnrIntegMesh.mDoubleSideSetsSlaveIndex(aGhostSetupData.mDblSideSetIndexInMesh(i)).push_back(tEnrIntegMesh.mDoubleSideSingleSideClusters.size());
                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back(tSlaveSideCluster);

                // create double side cluster
                mtk::Double_Side_Cluster* tDblSideCluster  = new mtk::Double_Side_Cluster(tMasterSideCluster.get(),tSlaveSideCluster.get(),tMasterSideCluster->mVerticesInCluster);

                // add to integration mesh
                tEnrIntegMesh.mDoubleSideClusters.push_back(tDblSideCluster);

                // add to the integration mesh
                tEnrIntegMesh.mDoubleSideSets(aGhostSetupData.mDblSideSetIndexInMesh(i))(j) = tDblSideCluster;

            }

            tEnrIntegMesh.commit_double_side_set(aGhostSetupData.mDblSideSetIndexInMesh(i));
            tEnrIntegMesh.set_double_side_set_colors(aGhostSetupData.mDblSideSetIndexInMesh(i),{{(moris_index)i}},{{(moris_index)i}});
        }

        tEnrIntegMesh.collect_all_sets();
    }

    // ----------------------------------------------------------------------------------

    bool
    Ghost_Stabilization::create_ghost(
            Ghost_Setup_Data  & aGhostSetupData,
            moris_index const & aFirstSubphase,
            moris_index const & aSecondSubphase,
            moris_index       & aTrivialFlag)
    {
        // Rules:
        // 1. Only create ghost facets between a subphase created inside an intersected
        //    cell and its neighbors.
        // 2. The owning processor of the master (first) subphase constructs the ghost facet.
        // 3. Construct from coarse to fine in HMR

        // make sure flag is set to true
        aTrivialFlag = 0;

        moris_index tFirstSubphaseId  = mXTKModel->get_subphase_id(aFirstSubphase);
        moris_index tSecondSubphaseId = mXTKModel->get_subphase_id(aSecondSubphase);

        MORIS_ASSERT(tFirstSubphaseId != tSecondSubphaseId,
                "Subphase neighbor relation inconsistent\n");

        // interpolation cell for this subphase
        moris_index tFirstInterpCell  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(aFirstSubphase);
        moris_index tSecondInterpCell = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(aSecondSubphase);

        // interpolation mesh
        Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKModel->get_enriched_interp_mesh();

        // MTK Cells
        mtk::Cell const & tFirstCell = tEnrInterpMesh.get_mtk_cell(tFirstInterpCell);
        mtk::Cell const & tSecondCell = tEnrInterpMesh.get_mtk_cell(tSecondInterpCell);

        // Levels
        uint tFirstLevel = tFirstCell.get_level();
        uint tSecondLevel = tSecondCell.get_level();

        // owners of interpolation cells
        moris_index tFirstOwnerIndex  = tFirstCell.get_owner();
        //        moris_index tSecondOwnerIndex = tSecondCell.get_owner();

        // proc rank
        moris_index tProcRank = par_rank();

        // if both subphases are not in child meshes return false
        if(!mXTKModel->subphase_is_in_child_mesh(aFirstSubphase) && !mXTKModel->subphase_is_in_child_mesh(aSecondSubphase))
        {
            return false;
        }

        // Check based on refinement level of subphases

        // if the first subphase is more coarse than the second subphase,
        // do not construct ghost
        if(tFirstLevel > tSecondLevel)
        {
            return false;
        }

        // if the first subphase is finer than the second subphase
        // do construct ghost if proc owns first subphase
        if(tFirstLevel < tSecondLevel)
        {
            if(tFirstOwnerIndex == tProcRank)
            {
                aTrivialFlag = 1;
                return true;
            }

            return false;
        }

        // first set of checks passed - subphases are on same refinement level

        // do not construct if first subphase ID is smaller than second subphase ID
        if(tFirstSubphaseId < tSecondSubphaseId)
        {
            return false;
        }

        // second set of checks passed - first subphase ID is larger than second one

        // do not construct if I do not own first subphase
        if(tFirstOwnerIndex != tProcRank )
        {
            return false;
        }

        // third set of checks passed - first subphase ID is larger than second one and it is owned by current proc

        return true;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr<Side_Cluster>
    Ghost_Stabilization::create_slave_side_cluster(
            Ghost_Setup_Data                         &  aGhostSetupData,
            Cell<Interpolation_Cell_Unzipped*>       & aEnrIpCells,
            uint                               const & aBulkIndex,
            uint                               const & aCellIndex)
    {
        // create a new side cluster for the slave
        std::shared_ptr<Side_Cluster> tSlaveSideCluster  = std::make_shared< Side_Cluster >();

        // give the cluster the enriched interpolation cell
        tSlaveSideCluster->mInterpolationCell  = aEnrIpCells(aGhostSetupData.mSlaveSideIpCells(aBulkIndex)(aCellIndex));

        // slave cluster is always trivial because the small facet is always the slave in the case of HMR hanging nodes
        tSlaveSideCluster->mTrivial = true;

        // add integration cell
        tSlaveSideCluster->mIntegrationCells  = {aEnrIpCells(aGhostSetupData.mSlaveSideIpCells(aBulkIndex)(aCellIndex))->get_base_cell()};

        // allocate space in integration cell side ordinals
        tSlaveSideCluster->mIntegrationCellSideOrdinals = {{aGhostSetupData.mSlaveSideIgCellSideOrds(aBulkIndex)(aCellIndex)}};

        // add vertices
        tSlaveSideCluster->mVerticesInCluster = tSlaveSideCluster->mIntegrationCells(0)->
                get_geometric_vertices_on_side_ordinal(tSlaveSideCluster->mIntegrationCellSideOrdinals(0));

        tSlaveSideCluster->mVerticesInCluster.size();

        return tSlaveSideCluster;
    }

    // ----------------------------------------------------------------------------------

    std::shared_ptr<Side_Cluster>
    Ghost_Stabilization::create_master_side_cluster(
            Ghost_Setup_Data                         & aGhostSetupData,
            Cell<Interpolation_Cell_Unzipped*>       & aEnrIpCells,
            uint                               const & aBulkIndex,
            uint                               const & aCellIndex,
            Side_Cluster* aSlaveSideCluster,
            moris_index & aCurrentIndex,
            moris_index & aCurrentId)
    {
        // create the master side cluster
        std::shared_ptr<Side_Cluster> tMasterSideCluster  = std::make_shared< Side_Cluster >();

        tMasterSideCluster->mInterpolationCell = aEnrIpCells(aGhostSetupData.mMasterSideIpCells(aBulkIndex)(aCellIndex));

        // Setup the master side cluster
        if(aGhostSetupData.mTrivialFlag(aBulkIndex)(aCellIndex) > 0)
        {
            // flag the master side as non-trivial
            tMasterSideCluster->mTrivial = false;

            // create new integration cell using the vertices on the slave facet and the adjacent vertices of the base interpolation cell
            moris::mtk::Cell* tNewIgCell = this->create_non_trivial_master_ig_cell(aGhostSetupData,aBulkIndex,aCellIndex,tMasterSideCluster.get(),aSlaveSideCluster,aCurrentIndex,aCurrentId);

            // get the local coordinates from a table
            Cell<Matrix<DDRMat>> tLocCoords;
            this->get_local_coords_on_transition_side(aGhostSetupData.mMasterSideIgCellSideOrds(aBulkIndex)(aCellIndex),
                    aGhostSetupData.mTransitionLocation(aBulkIndex)(aCellIndex),
                    tLocCoords);

            // add integration cell
            tMasterSideCluster->mIntegrationCells.push_back( tNewIgCell );

            // add side ordinal relative to the integration cell
            tMasterSideCluster->mIntegrationCellSideOrdinals = {{this->get_side_ordinals_for_non_trivial_master()}};

            // get the vertices on the side ordinal
            tMasterSideCluster->mVerticesInCluster = aSlaveSideCluster->mVerticesInCluster;

            // add the local coordinates
            tMasterSideCluster->mVertexLocalCoords = tLocCoords;

            // finalize
            tMasterSideCluster->finalize_setup();

            // verify  new cluster
            mtk::Mesh_Checker tCheck;
            MORIS_ASSERT(tCheck.verify_side_cluster(tMasterSideCluster.get(),mtk::Master_Slave::MASTER),"Invalid Side Cluster Created Check the local coordinates");

            // place the new ig cell in the background mesh
            mXTKModel->get_background_mesh().add_new_cell_to_mesh(tNewIgCell);
        }
        else
        {
            // flag the master side as trivial
            tMasterSideCluster->mTrivial = true;

            // add integration cell
            tMasterSideCluster->mIntegrationCells = {aEnrIpCells(aGhostSetupData.mMasterSideIpCells(aBulkIndex)(aCellIndex))->get_base_cell()};

            // add side ordinal relative to the integration cell
            tMasterSideCluster->mIntegrationCellSideOrdinals = {{aGhostSetupData.mMasterSideIgCellSideOrds(aBulkIndex)(aCellIndex)}};

            // add the vertices on the side ordinal
            tMasterSideCluster->mVerticesInCluster = tMasterSideCluster->mIntegrationCells(0)->get_geometric_vertices_on_side_ordinal(tMasterSideCluster->mIntegrationCellSideOrdinals(0));

            // finalize
            tMasterSideCluster->finalize_setup();
        }

        return tMasterSideCluster;
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell*
    Ghost_Stabilization::create_non_trivial_master_ig_cell(
            Ghost_Setup_Data       &  aGhostSetupData,
            uint             const & aBulkIndex,
            uint             const & aCellIndex,
            Side_Cluster           * aMasterSideCluster,
            Side_Cluster           * aSlaveSideCluster,
            moris_index            & aCurrentIndex,
            moris_index            & aCurrentId)
    {
        // get the vertices on the side for the slave side cluster
        moris::Cell<moris::mtk::Vertex const *> const & tSlaveVertices = aSlaveSideCluster->get_vertices_in_cluster();

        // get the side master interpolation cell
        Interpolation_Cell_Unzipped const *  tMasterIpCell =  aMasterSideCluster->mInterpolationCell;

        // base master cell
        moris::mtk::Cell const* tBaseMasterCell = tMasterIpCell->get_base_cell();

        // get the connectivity information from the cell
        moris::mtk::Cell_Info const * tCellInfo = tMasterIpCell->get_connectivity();

        // adjacent side ordinal on master
        uint tAdjFacetOrd = tCellInfo->get_adjacent_side_ordinal(aGhostSetupData.mMasterSideIgCellSideOrds(aBulkIndex)(aCellIndex));

        // setup the vertices and local coordinates of the vertices relative to the cell
        moris::Cell<moris::mtk::Vertex const *> tAdjVertices = tBaseMasterCell->get_geometric_vertices_on_side_ordinal(tAdjFacetOrd);

        //properly order the vertices
        moris::Cell<moris::mtk::Vertex const *> tPermutedSlaveVertices;
        moris::Cell<moris::mtk::Vertex const *> tPermutedAdjVertices;
        this->permute_slave_vertices(tSlaveVertices,tAdjVertices,tPermutedSlaveVertices,tPermutedAdjVertices);

        // New cell vertices setup (non-const)
        moris::Cell<moris::mtk::Vertex *> tCellVertices (tPermutedAdjVertices.size() + tPermutedSlaveVertices.size());

        moris::uint tCount = 0;
        for( auto iV:tPermutedAdjVertices)
        {
            tCellVertices(tCount++) = &mXTKModel->get_background_mesh().get_mtk_vertex(iV->get_index());
        }
        for( auto iV:tPermutedSlaveVertices)
        {
            tCellVertices(tCount++) = &mXTKModel->get_background_mesh().get_mtk_vertex(iV->get_index());
        }

        // create a new integration cell that does not have a child mesh association
        moris::mtk::Cell* tIgCell = new Cell_XTK_No_CM(aCurrentId, aCurrentIndex, tMasterIpCell->get_owner(), tCellInfo, tCellVertices);

        // increment current id and index
        aCurrentId++;
        aCurrentIndex++;

        return tIgCell;
    }

    // ----------------------------------------------------------------------------------

    void
    Ghost_Stabilization::permute_slave_vertices(
            moris::Cell<moris::mtk::Vertex const *> const & aSlaveVertices,
            moris::Cell<moris::mtk::Vertex const *> const & aMasterVertices,
            moris::Cell<moris::mtk::Vertex const *>       & aPermutedSlaveVertices,
            moris::Cell<moris::mtk::Vertex const *>       & aPermutedAdjMastVertices)
    {
        moris::uint tSpatialDim = mXTKModel->get_spatial_dim();

        switch(tSpatialDim)
        {
            case 2:
            {
                aPermutedSlaveVertices   = {aSlaveVertices(1),aSlaveVertices(0)};
                aPermutedAdjMastVertices = {aMasterVertices(0),aMasterVertices(1)};
                break;
            }
            default:
            {
                aPermutedSlaveVertices   = {aSlaveVertices(0),aSlaveVertices(3), aSlaveVertices(2), aSlaveVertices(1)};
                aPermutedAdjMastVertices = {aMasterVertices(0),aMasterVertices(3), aMasterVertices(2), aMasterVertices(1)};
            }
        }
    }
    // ----------------------------------------------------------------------------------
    void
    Ghost_Stabilization::get_local_coords_on_transition_side(
            moris_index const     & aMySideOrdinal,
            moris_index const     &  aTransitionLoc,
            Cell<Matrix<DDRMat>>  & aLocCoord)
    {
        moris::uint tTag = 100*mXTKModel->get_spatial_dim() + 10*aMySideOrdinal + aTransitionLoc;

        switch(tTag)
        {
            case 200:
                aLocCoord = {{{0,-1}},{{-1,-1}}};
                break;
            case 201:
                aLocCoord = {{{1,-1}},{{0,-1}}};
                break;
            case 211:
                aLocCoord = {{{1,0}},{{1,-1}}};
                break;
            case 213:
                aLocCoord = {{{1,1}},{{1,0}}};
                break;
            case 223:
                aLocCoord = {{{-1,1}},{{0,1}}};
                break;
            case 222:
                aLocCoord = {{{0,1}},{{1,1}}};
                break;
            case 230:
                aLocCoord = {{{-1,0}},{{-1,1}}};
                break;
            case 232:
                aLocCoord = {{{-1,-1}},{{-1,0}}};
                break;
            case 300:
                aLocCoord = {{{-1,-1,-1}},{{-1,-1,0}},{{0,-1,0}},{{0,-1,-1}}};
                break;
            case 301:
                aLocCoord = {{{0,-1,-1}},{{0,-1,0}},{{1,-1,0}},{{1,-1,-1}}};
                break;
            case 304:
                aLocCoord = {{{0,-1,0}},{{0,-1,1}},{{1,-1,1}},{{1,-1,0}}};
                break;
            case 305:
                aLocCoord = {{{-1,-1,0}},{{-1,-1,1}},{{0,-1,1}},{{0,-1,0}}};
                break;
            case 311:
                aLocCoord = {{{1,-1,-1}},{{1,-1,0}},{{1,0,0}},{{1,0,-1}}};
                break;
            case 313:
                aLocCoord = {{{1,0,-1}},{{1,0,0}},{{1,1,0}},{{1,1,-1}}};
                break;
            case 315:
                aLocCoord = {{{1,0,0}},{{1,0,1}},{{1,1,1}},{{1,1,0}}};
                break;
            case 317:
                aLocCoord = {{{1,-1,0}},{{1,-1,1}},{{1,0,1}},{{1,0,0}}};
                break;
            case 322:
                aLocCoord = {{{1,1,-1}},{{1,1,0}},{{0,1,0}},{{0,1,-1}}};
                break;
            case 323:
                aLocCoord = {{{0,1,-1}},{{0,1,0}},{{-1,1,0}},{{-1,1,-1}}};
                break;
            case 326:
                aLocCoord = {{{0,1,0}},{{0,1,1}},{{-1,1,1}},{{-1,1,0}}};
                break;
            case 327:
                aLocCoord = {{{1,1,0}},{{1,1,1}},{{0,1,1}},{{0,1,0}}};
                break;
            case 330:
                aLocCoord = {{{-1,0,-1}},{{-1,1,-1}},{{-1,1,0}},{{-1,0,0}}};
                break;
            case 332:
                aLocCoord = {{{-1,-1,-1}},{{-1,0,-1}},{{-1,0,0}},{{-1,-1,0}}};
                break;
            case 334:
                aLocCoord = {{{-1,-1,0}},{{-1,0,0}},{{-1,0,1}},{{-1,-1,1}}};
                break;
            case 336:
                aLocCoord = {{{-1,0,0}},{{-1,1,0}},{{-1,1,1}},{{-1,0,1}}};
                break;
            case 340:
                aLocCoord = {{{0,0,-1}},{{1,0,-1}},{{1,1,-1}},{{0,1,-1}}};
                break;
            case 341:
                aLocCoord = {{{-1,0,-1}},{{0,0,-1}},{{0,1,-1}},{{-1,1,-1}}};
                break;
            case 342:
                aLocCoord = {{{-1,-1,-1}},{{0,-1,-1}},{{0,0,-1}},{{-1,0,-1}}};
                break;
            case 343:
                aLocCoord = {{{0,-1,-1}},{{1,-1,-1}},{{1,0,-1}},{{0,0,-1}}};
                break;
            case 354:
                aLocCoord = {{{-1,-1, 1}},{{-1,0,1}},{{0,0,1}},{{0,-1,1}}};
                break;
            case 355:
                aLocCoord = {{{0,-1,1}},{{0,0,1}},{{1,0,1}},{{1,-1, 1}}};
                break;
            case 356:
                aLocCoord = {{{0,0,1}},{{0,1,1}},{{1,1,1}},{{1,0,1}}};
                break;
            case 357:
                aLocCoord = {{{-1,0,1}},{{-1,1,1}},{{0,1,1}},{{0,0,1}}};
                break;
            default:
                MORIS_ERROR(0,"Invalid tag (100*spatial dim + 10 * side ord + transition location)");
        }
    }
    // ----------------------------------------------------------------------------------

    moris_index
    Ghost_Stabilization::get_side_ordinals_for_non_trivial_master()
    {
        if(mXTKModel->get_spatial_dim() == 2)
        {
            return 2;
        }
        else if (mXTKModel->get_spatial_dim() == 3)
        {
            return 5;
        }
        else
        {
            MORIS_ERROR(0,"Invalid spatial dimension for ghost");
            return 0;
        }
    }
}
