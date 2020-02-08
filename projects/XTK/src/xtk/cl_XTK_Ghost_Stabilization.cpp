/*
 * cl_XTK_Ghost_Stabilization.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "fn_generate_shared_face_element_graph.hpp"

namespace xtk
{
    void
    Ghost_Stabilization::setup_ghost_stabilization()
    {
        Ghost_Setup_Data tGhostSetupData;

        // construct trivial subphase interpolation cells
        this->construct_ip_ig_cells_for_ghost_side_clusters(tGhostSetupData);

        // setup the side sets
        this->declare_ghost_double_side_sets_in_mesh(tGhostSetupData);

        // Construct Ghost Double Side Clusters and Sets
        this->construct_ghost_double_side_sets_in_mesh(tGhostSetupData);
    }

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
            tFacetTopo = CellTopology::LINE2;
        }

        moris_index tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set(tDSSIndexInMesh,"ghost_ss_" + std::to_string(aBulkPhase));
        tEnrIgMesh.create_block_set_from_cells_of_side_set(tSSIndex,"ghost_bs_" + std::to_string(aBulkPhase), tFacetTopo);

        tEnrIgMesh.setup_color_to_set();
        tEnrIgMesh.collect_all_sets();
    }

    std::string
    Ghost_Stabilization::get_ghost_dbl_side_set_name(moris_index aBulkPhase)
    {
        MORIS_ASSERT(aBulkPhase < (moris_index)mXTKModel->get_geom_engine().get_num_bulk_phase(),"Bulk Phase index out of bounds.");
        return "ghost_p" + std::to_string(aBulkPhase);
    }

    void
    Ghost_Stabilization::construct_ip_ig_cells_for_ghost_side_clusters(Ghost_Setup_Data & aGhostSetupData)
    {
        // access enriched integration mesh and enriched interp mesh
        Enriched_Integration_Mesh & tEnrIgMesh = mXTKModel->get_enriched_integ_mesh(0);
        Enriched_Interpolation_Mesh & tEnrIpMesh = mXTKModel->get_enriched_interp_mesh(0);

        // get the groupings of interpolation cells
        Cell<Interpolation_Cell_Unzipped *>       tOwnedInterpCells;
        Cell<Cell<Interpolation_Cell_Unzipped *>> tNotOwnedInterpCells;
        Cell<uint>                                tProcRanks;
        tEnrIpMesh.get_owned_and_not_owned_enriched_interpolation_cells(tOwnedInterpCells,tNotOwnedInterpCells, tProcRanks);

        // number of interpolation cells
        moris::uint tCurrentNewInterpCellIndex = tEnrIpMesh.get_num_entities(EntityRank::ELEMENT);

        // allocate data in ghost setup data
        aGhostSetupData.mSubphaseIndexToInterpolationCellIndex.resize(tCurrentNewInterpCellIndex,MORIS_INDEX_MAX);

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
                aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tCurrentNewInterpCellIndex;
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
                    aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tCurrentNewInterpCellIndex;
                    tCurrentNewInterpCellIndex++;
                }
                else
                {
                    aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase) = tNotOwnedInterpCells(iP)(iC)->get_index();
                }
            }
        }

        // prepare data to request
        Cell<Matrix<IndexMat>> tBaseEnrInterpCellId(tNonTrivialNotOwnedInterpCells.size());

        for(moris::size_t iP = 0; iP<tNonTrivialNotOwnedInterpCells.size(); iP++)
        {
            tBaseEnrInterpCellId(iP).resize(1,tNonTrivialNotOwnedInterpCells(iP).size());

            for(moris::size_t iC = 0; iC<tNonTrivialNotOwnedInterpCells(iP).size(); iC++)
            {
                tBaseEnrInterpCellId(iP)(iC) = tNonTrivialNotOwnedInterpCells(iP)(iC)->get_id();
            }
        }

        // send requests
        moris::uint tMPITag = 301;
        mXTKModel->send_outward_requests(tMPITag, tProcRanks,tBaseEnrInterpCellId);

        barrier();

        // receive requests
        Cell<Matrix<IndexMat>> tReceivedBaseEnrCellIds;
        Cell<uint> tProcsReceivedFrom;
        mXTKModel->inward_receive_requests(tMPITag, 1, tReceivedBaseEnrCellIds, tProcsReceivedFrom);

        // prepare answers
        Cell<Matrix<IndexMat>> tEnrCellIds;
        this->prepare_ip_cell_id_answers(tReceivedBaseEnrCellIds,tNewNonTrivialOwnedInterpCellsIds,tEnrCellIds, tBaseEnrIdToIndexInNonTrivialOwned);

        // return information
        mXTKModel->return_request_answers(tMPITag+1, tEnrCellIds, tProcsReceivedFrom);

        barrier();

        // receive the answers
        Cell<Matrix<IndexMat>> tReceivedEnrCellIds;
        mXTKModel->inward_receive_request_answers(tMPITag+1,1,tProcRanks,tReceivedEnrCellIds);

        // allocate space in interpolation cells of enriched interpolation mesh
        tEnrIpMesh.mEnrichedInterpCells.resize(tNumNewInterpCellsNotOwned + tNumNewInterpCellsOwned + tEnrIpMesh.mEnrichedInterpCells.size());

        // add interpolation cells for ghost to enriched interp mesh
        this->create_not_owned_ghost_ip_cells(aGhostSetupData,tEnrIpMesh,tNonTrivialNotOwnedInterpCells,tReceivedEnrCellIds);

        // create owned ghost ip cells
        this->create_owned_ghost_ip_cells(aGhostSetupData,tEnrIpMesh,tNonTrivialOwnedInterpCells,tNewNonTrivialOwnedInterpCellsIds);

    }

    void
    Ghost_Stabilization::prepare_ip_cell_id_answers( Cell<Matrix<IndexMat>> & aReceivedEnrCellIds,
                                                     Cell<moris_id>         & aNewInterpCellIds,
                                                     Cell<Matrix<IndexMat>> & aEnrCellIds,
                                                     std::unordered_map<moris_id, moris_id> & aBaseEnrIdToIndexInNonTrivialOwned)
    {
        aEnrCellIds.resize(aReceivedEnrCellIds.size());

        for(moris::uint iP = 0; iP< aReceivedEnrCellIds.size(); iP++)
        {
            aEnrCellIds(iP).resize(1,aReceivedEnrCellIds(iP).numel());

            for(moris::uint iC = 0; iC< aReceivedEnrCellIds(iP).numel(); iC++)
            {
                auto tIter = aBaseEnrIdToIndexInNonTrivialOwned.find(aReceivedEnrCellIds(iP)(iC));
                MORIS_ASSERT(tIter != aBaseEnrIdToIndexInNonTrivialOwned.end(),"Enriched cell id not in map");

                aEnrCellIds(iP)(iC) = aNewInterpCellIds(tIter->second);

            }
        }
    }

    void
    Ghost_Stabilization::create_not_owned_ghost_ip_cells( Ghost_Setup_Data &                                      aGhostSetupData,
                                                          Enriched_Interpolation_Mesh &                           aEnrInterpMesh,
                                                          Cell<Cell<Interpolation_Cell_Unzipped *>> const & aNonTrivialNotOwnedInterpCells,
                                                          Cell<Matrix<IndexMat>>                          const & aReceivedEnrCellIds)
    {
        // iterate through received data
        for(moris::uint i = 0; i < aNonTrivialNotOwnedInterpCells.size(); i++)
        {
            uint tNumReceivedReqs = aNonTrivialNotOwnedInterpCells(i).size();

            // iterate through received requests
            for(moris::uint j = 0; j < tNumReceivedReqs; j++)
            {
                moris_index tReceivedEnrCellId = aReceivedEnrCellIds(i)(j);
                moris_index tSubphase          = aNonTrivialNotOwnedInterpCells(i)(j)->get_subphase_index();
                moris_index tGhostCellIpIndex  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase);

                MORIS_ASSERT(tGhostCellIpIndex != MORIS_INDEX_MAX,"Max index on not owned ghost ip cell. This could be a communication error.");

                aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex) = new Interpolation_Cell_Unzipped(aNonTrivialNotOwnedInterpCells(i)(j)->get_base_cell(),
                                                                                                         aNonTrivialNotOwnedInterpCells(i)(j)->get_subphase_index(),
                                                                                                         aNonTrivialNotOwnedInterpCells(i)(j)->get_bulkphase_index(),
                                                                                                         aNonTrivialNotOwnedInterpCells(i)(j)->get_id(),
                                                                                                         tGhostCellIpIndex,
                                                                                                         tReceivedEnrCellId,
                                                                                                         aNonTrivialNotOwnedInterpCells(i)(j)->get_connectivity());

                aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex)->set_vertices(aNonTrivialNotOwnedInterpCells(i)(j)->get_xtk_interpolation_vertices());

            }
        }
    }

    void
    Ghost_Stabilization::create_owned_ghost_ip_cells( Ghost_Setup_Data &                          aGhostSetupData,
                                                      Enriched_Interpolation_Mesh &               aEnrInterpMesh,
                                                      Cell<Interpolation_Cell_Unzipped *> & aNonTrivialOwnedInterpCells,
                                                      Cell<moris_id>                      & aEnrCellIds)
    {
        // iterate through received data
        for(moris::uint i = 0; i < aNonTrivialOwnedInterpCells.size(); i++)
        {
                moris_index tEnrCellId = aEnrCellIds(i);
                moris_index tSubphase          = aNonTrivialOwnedInterpCells(i)->get_subphase_index();
                moris_index tGhostCellIpIndex  = aGhostSetupData.mSubphaseIndexToInterpolationCellIndex(tSubphase);

                aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex) = new Interpolation_Cell_Unzipped(aNonTrivialOwnedInterpCells(i)->get_base_cell(),
                                                                                                         aNonTrivialOwnedInterpCells(i)->get_subphase_index(),
                                                                                                         aNonTrivialOwnedInterpCells(i)->get_bulkphase_index(),
                                                                                                         tEnrCellId,
                                                                                                         tGhostCellIpIndex,
                                                                                                         aNonTrivialOwnedInterpCells(i)->get_owner(),
                                                                                                         aNonTrivialOwnedInterpCells(i)->get_connectivity());

                aEnrInterpMesh.mEnrichedInterpCells(tGhostCellIpIndex)->set_vertices(aNonTrivialOwnedInterpCells(i)->get_xtk_interpolation_vertices());
        }
    }

    void
    Ghost_Stabilization::declare_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData)
    {
        uint tNumBulkPhases = mXTKModel->get_geom_engine().get_num_bulk_phase();

        Cell<std::string> tGhostDoubleSideNames(tNumBulkPhases);

        for(moris::moris_index iP0 = 0; iP0 <(moris_index) tNumBulkPhases; iP0++)
        {

            std::string tGhostSideSetName =this->get_ghost_dbl_side_set_name(iP0);

            tGhostDoubleSideNames(iP0) = tGhostSideSetName;

        }
        aGhostSetupData.mDblSideSetIndexInMesh = mXTKModel->get_enriched_integ_mesh(0).register_double_side_set_names(tGhostDoubleSideNames);
    }

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
        moris::uint tNumBulkPhases = mXTKModel->get_geom_engine().get_num_bulk_phase();

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

                // add to side clusters the integration mesh
                tEnrIntegMesh.mDoubleSideSetsMasterIndex(aGhostSetupData.mDblSideSetIndexInMesh(i)).push_back(tEnrIntegMesh.mDoubleSideSingleSideClusters.size());
                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back(tMasterSideCluster);
                tEnrIntegMesh.mDoubleSideSetsSlaveIndex(aGhostSetupData.mDblSideSetIndexInMesh(i)).push_back(tEnrIntegMesh.mDoubleSideSingleSideClusters.size());
                tEnrIntegMesh.mDoubleSideSingleSideClusters.push_back(tSlaveSideCluster);

                 // create double side cluster
                mtk::Double_Side_Cluster* tDblSideCluster  = new mtk::Double_Side_Cluster(tMasterSideCluster.get(),tSlaveSideCluster.get(),tMasterSideCluster->mVerticesInCluster);

                // add to the integration mesh
                tEnrIntegMesh.mDoubleSideSets(aGhostSetupData.mDblSideSetIndexInMesh(i))(j) = tDblSideCluster;

            }

            tEnrIntegMesh.commit_double_side_set(aGhostSetupData.mDblSideSetIndexInMesh(i));
            tEnrIntegMesh.set_double_side_set_colors(aGhostSetupData.mDblSideSetIndexInMesh(i),{{(moris_index)i}},{{(moris_index)i}});
        }

    }

    bool
    Ghost_Stabilization::create_ghost(Ghost_Setup_Data &  aGhostSetupData,
                                      moris_index const & aFirstSubphase,
                                      moris_index const & aSecondSubphase,
                                      moris_index &       aTrivialFlag)
    {
        // make sure flag is set to true
        aTrivialFlag = 0;

        moris_index tFirstSubphaseId  = mXTKModel->get_subphase_id(aFirstSubphase);
        moris_index tSecondSubphaseId = mXTKModel->get_subphase_id(aSecondSubphase);

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
        moris_index tSecondOwnerIndex = tSecondCell.get_owner();

        // proc rank
        moris_index tProcRank = par_rank();

        if(!mXTKModel->subphase_is_in_child_mesh(aFirstSubphase) && !mXTKModel->subphase_is_in_child_mesh(aSecondSubphase))
        {
            return false;
        }

        // HMR large to small facet transition and I own it
        // always create from large to small facet
        if(tFirstLevel < tSecondLevel && tFirstOwnerIndex == tProcRank)
        {
            aTrivialFlag = 1;
            return true;
        }

        // if the second subphase is more coarse than the first subphase
        if(tFirstLevel > tSecondLevel)
        {
            return false;
        }

        // normal checks
        if(tFirstSubphaseId < tSecondSubphaseId)
        {
            return false;
        }

        if(tFirstOwnerIndex > tSecondOwnerIndex)
        {
            return false;
        }

        return true;
    }

    std::shared_ptr<Side_Cluster>
    Ghost_Stabilization::create_slave_side_cluster(Ghost_Setup_Data &  aGhostSetupData,
                                                   Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                                                   uint const & aBulkIndex,
                                                   uint const & aCellIndex)
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
        tSlaveSideCluster->mVerticesInCluster = tSlaveSideCluster->mIntegrationCells(0)->get_vertices_on_side_ordinal(tSlaveSideCluster->mIntegrationCellSideOrdinals(0));

        return tSlaveSideCluster;
    }

    std::shared_ptr<Side_Cluster>
    Ghost_Stabilization::create_master_side_cluster(Ghost_Setup_Data &  aGhostSetupData,
                                                    Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                                                    uint const & aBulkIndex,
                                                    uint const & aCellIndex,
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
            tMasterSideCluster->mVerticesInCluster = tMasterSideCluster->mIntegrationCells(0)->get_vertices_on_side_ordinal(tMasterSideCluster->mIntegrationCellSideOrdinals(0));

            // finalize
            tMasterSideCluster->finalize_setup();
        }


        return tMasterSideCluster;
    }

    moris::mtk::Cell*
    Ghost_Stabilization::create_non_trivial_master_ig_cell(Ghost_Setup_Data &  aGhostSetupData,
                                                           uint const & aBulkIndex,
                                                           uint const & aCellIndex,
                                                           Side_Cluster* aMasterSideCluster,
                                                           Side_Cluster* aSlaveSideCluster,
                                                           moris_index & aCurrentIndex,
                                                           moris_index & aCurrentId)
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
        moris::Cell<moris::mtk::Vertex const *> tAdjVertices = tBaseMasterCell->get_vertices_on_side_ordinal(tAdjFacetOrd);

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



    void
    Ghost_Stabilization::permute_slave_vertices(moris::Cell<moris::mtk::Vertex const *> const & aSlaveVertices,
                                                moris::Cell<moris::mtk::Vertex const *> const & aMasterVertices,
                                                moris::Cell<moris::mtk::Vertex const *>  & aPermutedSlaveVertices,
                                                moris::Cell<moris::mtk::Vertex const *>  & aPermutedAdjMastVertices)
    {
        moris::uint tSpatialDim = mXTKModel->get_spatial_dim();

        switch(tSpatialDim)
        {
            case(2):
                aPermutedSlaveVertices   = {aSlaveVertices(1),aSlaveVertices(0)};
                aPermutedAdjMastVertices = {aMasterVertices(0),aMasterVertices(1)};
                break;
            default:
                aPermutedSlaveVertices   = {aSlaveVertices(0),aSlaveVertices(3), aSlaveVertices(2), aSlaveVertices(1)};
                aPermutedAdjMastVertices = {aMasterVertices(0),aMasterVertices(3), aMasterVertices(2), aMasterVertices(1)};
                break;
        }
    }


    void
    Ghost_Stabilization::get_local_coords_on_transition_side(moris_index const & aMySideOrdinal,
                                                             moris_index const &  aTransitionLoc,
                                                             Cell<Matrix<DDRMat>>    & aLocCoord)
    {
        moris::uint tTag = 100*mXTKModel->get_spatial_dim() + 10*aMySideOrdinal + aTransitionLoc;
        switch(tTag)
        {
            case(200):
                    aLocCoord = {{{0,-1}},{{-1,-1}}};
                    break;
            case(201):
                    aLocCoord = {{{1,-1}},{{0,-1}}};
                    break;
            case(211):
                    aLocCoord = {{{1,0}},{{1,-1}}};
                    break;
            case(213):
                    aLocCoord = {{{1,1}},{{1,0}}};
                    break;
            case(223):
                    aLocCoord = {{{-1,1}},{{0,1}}};
                    break;
            case(222):
                    aLocCoord = {{{0,1}},{{1,1}}};
                    break;
            case(230):
                    aLocCoord = {{{-1,0}},{{-1,1}}};
                    break;
            case(232):
                    aLocCoord = {{{-1,-1}},{{-1,0}}};
                    break;
            case(300):
                    aLocCoord = {{{-1,-1,-1}},{{-1,-1,0}},{{0,-1,0}},{{0,-1,-1}}};
                    break;
            case(301):
                    aLocCoord = {{{0,-1,-1}},{{0,-1,0}},{{1,-1,0}},{{1,-1,-1}}};
                    break;
            case(304):
                    aLocCoord = {{{0,-1,0}},{{0,-1,1}},{{1,-1,1}},{{1,-1,0}}};
                    break;
            case(305):
                    aLocCoord = {{{-1,-1,0}},{{-1,-1,1}},{{0,-1,1}},{{0,-1,0}}};
                    break;
            case(311):
                    aLocCoord = {{{1,-1,-1}},{{1,-1,0}},{{1,0,0}},{{1,0,-1}}};
                    break;
            case(313):
                    aLocCoord = {{{1,0,-1}},{{1,0,0}},{{1,1,0}},{{1,1,-1}}};
                    break;
            case(315):
                    aLocCoord = {{{1,0,0}},{{1,0,1}},{{1,1,1}},{{1,1,0}}};
                    break;
            case(317):
                    aLocCoord = {{{1,-1,0}},{{1,-1,1}},{{1,0,1}},{{1,0,0}}};
                    break;
            case(322):
                    aLocCoord = {{{1,1,-1}},{{1,1,0}},{{0,1,0}},{{0,1,-1}}};
                    break;
            case(323):
                    aLocCoord = {{{0,1,-1}},{{0,1,0}},{{-1,1,0}},{{-1,1,-1}}};
                    break;
            case(326):
                    aLocCoord = {{{0,1,0}},{{0,1,1}},{{-1,1,1}},{{-1,1,0}}};
                    break;
            case(327):
                    aLocCoord = {{{1,1,0}},{{1,1,1}},{{0,1,1}},{{0,1,0}}};
                    break;
            case(330):
                    aLocCoord = {{{-1,0,-1}},{{-1,1,-1}},{{-1,1,0}},{{-1,0,0}}};
                    break;
            case(332):
                    aLocCoord = {{{-1,-1,-1}},{{-1,0,-1}},{{-1,0,0}},{{-1,-1,0}}};
                    break;
            case(334):
                    aLocCoord = {{{-1,-1,0}},{{-1,0,0}},{{-1,0,1}},{{-1,-1,1}}};
                    break;
            case(336):
                    aLocCoord = {{{-1,0,0}},{{-1,1,0}},{{-1,1,1}},{{-1,0,1}}};
                    break;
            case(340):
                    aLocCoord = {{{0,0,-1}},{{1,0,-1}},{{1,1,-1}},{{0,1,-1}}};
                    break;
            case(341):
                    aLocCoord = {{{-1,0,-1}},{{0,0,-1}},{{0,1,-1}},{{-1,1,-1}}};
                    break;
            case(342):
                    aLocCoord = {{{-1,-1,-1}},{{0,-1,-1}},{{0,0,-1}},{{-1,0,-1}}};
                    break;
            case(343):
                    aLocCoord = {{{0,-1,-1}},{{1,-1,-1}},{{1,0,-1}},{{0,0,-1}}};
                    break;
            case(354):
                    aLocCoord = {{{-1,-1, 1}},{{-1,0,1}},{{0,0,1}},{{0,-1,1}}};
                    break;
            case(355):
                    aLocCoord = {{{0,-1,1}},{{0,0,1}},{{1,0,1}},{{1,-1, 1}}};
                    break;
            case(356):
                    aLocCoord = {{{0,0,1}},{{0,1,1}},{{1,1,1}},{{1,0,1}}};
                    break;
            case(357):
                    aLocCoord = {{{-1,0,1}},{{-1,1,1}},{{0,1,1}},{{0,0,1}}};
                    break;
            default:
                MORIS_ERROR(0,"Invalid tag (100*spatial dim + 10 * side ord + transition location)");
                break;
        }
    }


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
