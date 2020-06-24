/*
 * cl_XTK_Ghost_Penalization.hpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_


#include "cl_XTK_Model.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
using namespace moris;

namespace moris
{
    namespace mtk
    {
        class Cell;
    }
}
namespace xtk
{
    /*
     * Data used while setting up the ghost stabilization (to not have to pass everything around in function inputs)
     */
    struct Ghost_Setup_Data
    {
            Cell<moris_index> mSubphaseIndexToInterpolationCellIndex;
            Cell<moris_index> mDblSideSetIndexInMesh;


            // interpolation cells
            Cell<Cell<moris_index>> mMasterSideIpCells;
            Cell<Cell<moris_index>> mSlaveSideIpCells;

            // Side ordinals
            Cell<Cell<moris_index>> mMasterSideIgCellSideOrds;
            Cell<Cell<moris_index>> mSlaveSideIgCellSideOrds;

            // trivial clusters are ones that are not on a hanging node interface
            Cell<Cell<moris_index>> mTrivialFlag;
            Cell<Cell<moris_index>> mTransitionLocation;
    };

    class Model;

    class Ghost_Stabilization
    {
        public:
            // ----------------------------------------------------------------------------------
            Ghost_Stabilization();

            // ----------------------------------------------------------------------------------
            Ghost_Stabilization( Model* aXTKModel );

            // ----------------------------------------------------------------------------------
            void
            setup_ghost_stabilization();

            // ----------------------------------------------------------------------------------
            void
            visualize_ghost_on_mesh(moris_index aBulkPhase);

            // ----------------------------------------------------------------------------------
            std::string
            get_ghost_dbl_side_set_name(moris_index aBulkPhase);

        private:
            // ----------------------------------------------------------------------------------
            /*!
             * Construct interpolation cell for trivial cell cluster built on child mesh subphases.
             *
             * This is done so that we can still have 1 cluster per interpolation cell.
             */
            void
            construct_ip_ig_cells_for_ghost_side_clusters(Ghost_Setup_Data & aGhostSetupData);
            // ----------------------------------------------------------------------------------
            /*
             * identifies the vertices in the aura that will interpolate into a ghost facet.
             * In general, these vertices do not have a vertex interpolation so the vertex
             * interpolation information needs to be communicated. This function handles that.
             */
            void
            identify_and_setup_aura_vertices_in_ghost();
            // ----------------------------------------------------------------------------------
            /*!
             * Prepare requests for t-matrices in the ghost.
             * @param[in]  aNotOwnedIpCellsInGhost List of ip cells that have nodes on the aura and are in ghost
             * @param[out] aNotOwnedIPVertsToProcs - Vertices that need to have t-matrix communicated
             * @param[out] aProcRanks - process ranks
             * @param[out] aProcRankToDataIndex - proc rank in comm to index in aProcRanks
             */
            void
            prepare_interpolation_vertex_t_matrix_requests(
                    std::unordered_map<moris_index,bool> const & aNotOwnedIpCellsInGhost,
                    Cell<Matrix<IndexMat>>                     & aNotOwnedIPVertIndsToProcs,
                    Cell<Matrix<IndexMat>>                     & aNotOwnedBGIPVertsIdsToProcs,
                    Cell<Matrix<IndexMat>>                     & aNotOwnedIpCellIdToProcs,
                    Cell<Matrix<IndexMat>>                     & aNotOwnedEnrichedCellBulkPhaseToProcs,
                    Cell<uint>                                 & aProcRanks,
                    std::unordered_map<moris_id,moris_id>      & aProcRankToDataIndex);
            // ----------------------------------------------------------------------------------
            /*
             * Using a background vertex id, and a enriched interpolation cell id
             * find the corresponding enriched interpolation vertex
             */
            moris_index
            get_enriched_interpolation_vertex(moris_index const & aBGVertId,
                                              moris_index const & aEnrichedIpCellIndex);


            // ----------------------------------------------------------------------------------
            /*
             * Packages t-matrix  request into mpi ready data structure (i.e. sparse matrix for returning back
             * to requesting processor.
             */
            void
            prepare_t_matrix_request_answers(
                    Cell<Matrix<IndexMat>> const & aRequestedBgVertexIds,
                    Cell<Matrix<IndexMat>> const & aRequestedIpCellIds,
                    Cell<Matrix<IndexMat>> const & aIpCellBulkPhases,
                    Cell<Matrix<DDRMat>>   & aTMatrixWeights,
                    Cell<Matrix<IndexMat>> & aTMatrixIndices,
                    Cell<Matrix<IndexMat>> & aBasisOwners,
                    Cell<Matrix<IndexMat>> & aTMatrixOffsets);

            void
            add_vertex_interpolation_to_communication_data(moris::uint      & aCount,
                                                           Vertex_Enrichment* aInterpolation,
                                                           Matrix<DDRMat>   & aTMatrixWeights,
                                                           Matrix<IndexMat> & aTMatrixIndices,
                                                           Matrix<IndexMat> & aTMatrixOwners,
                                                           Matrix<IndexMat> & aTMatrixOffsets);

            void
            extract_vertex_interpolation_from_communication_data(
                    moris::uint      const & aNumVerts,
                    Matrix<DDRMat>   const & aTMatrixWeights,
                    Matrix<IndexMat> const & aTMatrixIndices,
                    Matrix<IndexMat> const & aTMatrixOwners,
                    Matrix<IndexMat> const & aTMatrixOffsets,
                    Cell<Matrix<DDRMat>>   & aExtractedTMatrixWeights,
                    Cell<Matrix<IndexMat>> & aExtractedTMatrixIndices,
                    Cell<Matrix<IndexMat>> & aExtractedBasisOwners);

            void
            handle_received_interpolation_data(
                    Cell<Matrix<IndexMat>> const & aNotOwnedIPVertIndsToProcs,
                    Cell<Matrix<IndexMat>> const & aNotOwnedEnrichedCellBulkPhaseToProcs,
                    Cell<Matrix<DDRMat>>   const & aRequestedTMatrixWeights,
                    Cell<Matrix<IndexMat>> const & aRequestedTMatrixIndices,
                    Cell<Matrix<IndexMat>> const & aRequestedBasisOwners,
                    Cell<Matrix<IndexMat>> const & aRequestedTMatrixOffsets);

            // ----------------------------------------------------------------------------------
            void
            prepare_ip_cell_id_answers(
                    Cell<Matrix<IndexMat>>                 & aReceivedEnrCellIds,
                    Cell<moris_id>                         & aNewInterpCellIds,
                    Cell<Matrix<IndexMat>>                 & aEnrCellIds,
                    std::unordered_map<moris_id, moris_id> & aBaseEnrIdToIndexInNonTrivialOwned);
            // ----------------------------------------------------------------------------------
            void
            create_not_owned_ghost_ip_cells(
                    Ghost_Setup_Data &                                 aGhostSetupData,
                    Enriched_Interpolation_Mesh &                      aEnrInterpMesh,
                    Cell<Cell<Interpolation_Cell_Unzipped  *>> const & aNonTrivialNotOwnedInterpCells,
                    Cell<Matrix<IndexMat>>                     const & aReceivedEnrCellIds);
            // ----------------------------------------------------------------------------------
            void
            create_owned_ghost_ip_cells(
                    Ghost_Setup_Data &                    aGhostSetupData,
                    Enriched_Interpolation_Mesh &         aEnrInterpMesh,
                    Cell<Interpolation_Cell_Unzipped *> & aNonTrivialOwnedInterpCells,
                    Cell<moris_id>                      & aEnrCellIds);
            // ----------------------------------------------------------------------------------
            bool
            verify_ghost_communication(Ghost_Setup_Data & aGhostSetupData);

            // ----------------------------------------------------------------------------------
            void
            declare_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData);

            // ----------------------------------------------------------------------------------
            void
            construct_ghost_double_side_sets_in_mesh(Ghost_Setup_Data & aGhostSetupData);

            // ----------------------------------------------------------------------------------
            bool
            create_ghost(
                    Ghost_Setup_Data &  aGhostSetupData,
                    moris_index const & aFirstSubphase,
                    moris_index const & aSecondSubphase,
                    moris_index &       aTrivialFlag);

            // ----------------------------------------------------------------------------------
            std::shared_ptr<Side_Cluster>
            create_slave_side_cluster(
                    Ghost_Setup_Data &  aGhostSetupData,
                    Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                    uint const & aBulkIndex,
                    uint const & aCellIndex);

            // ----------------------------------------------------------------------------------
            std::shared_ptr<Side_Cluster>
            create_master_side_cluster(
                    Ghost_Setup_Data &  aGhostSetupData,
                    Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                    uint const & aBulkIndex,
                    uint const & aCellIndex,
                    Side_Cluster* aSlaveSideCluster,
                    moris_index & aCurrentIndex,
                    moris_index & aCurrentId);

            // ----------------------------------------------------------------------------------
            moris::mtk::Cell*
            create_non_trivial_master_ig_cell(
                    Ghost_Setup_Data &  aGhostSetupData,
                    uint const & aBulkIndex,
                    uint const & aCellIndex,
                    Side_Cluster* aMasterSideCluster,
                    Side_Cluster* aSlaveSideCluster,
                    moris_index & aCurrentIndex,
                    moris_index & aCurrentId);

            // ----------------------------------------------------------------------------------
            void
            permute_slave_vertices(
                    moris::Cell<moris::mtk::Vertex const *> const & aSlaveVertices,
                    moris::Cell<moris::mtk::Vertex const *> const & aMasterVertices,
                    moris::Cell<moris::mtk::Vertex const *>  & aPermutedSlaveVertices,
                    moris::Cell<moris::mtk::Vertex const *>  & aPermutedAdjMastVertices);
            // ----------------------------------------------------------------------------------
            void
            get_local_coords_on_transition_side(moris_index const    & aMySideOrdinal,
                    moris_index const    & aTransitionLoc,
                    Cell<Matrix<DDRMat>> & aLocCoord);
            // ----------------------------------------------------------------------------------
            moris_index get_side_ordinals_for_non_trivial_master();

        private:
            Model* mXTKModel; /*Pointer to the model*/

    };

}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_ */
