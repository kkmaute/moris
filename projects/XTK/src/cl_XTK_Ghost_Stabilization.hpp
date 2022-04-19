/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Ghost_Penalization.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_


#include "cl_XTK_Model.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "fn_TOL_Capacities.hpp"
#include "cl_TOL_Memory_Map.hpp"
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
    class Cell_XTK_No_CM;
    class Model;

    // ----------------------------------------------------------------------------------
    
    /**
     * @brief Data used while setting up the ghost stabilization 
     * (to not have to pass everything around in function inputs)
     * 
     */
    struct Ghost_Setup_Data
    {
        Cell< moris_index > mSubphaseIndexToInterpolationCellIndex; // input: subphase index || output: enriched interpolation cell index
        Cell< moris_index > mDblSideSetIndexInMesh;                 // input: index of ghost side set || output: corresponding index of dbl. side set as saved on enriched IG mesh

        // ----------------------------------------------------------------------------------
        // Maps for old Ghost

        // interpolation cells
        Cell< Cell< moris_index > > mMasterSideIpCells; // input: bulk phase index || output: list of enriched master IP cells used for ghost facet construction 
        Cell< Cell< moris_index > > mSlaveSideIpCells;  // input: bulk phase index || output: list of corresponding enriched slave IP cells used for ghost facet construction  

        // Side ordinals
        Cell< Cell< moris_index > > mMasterSideIgCellSideOrds; // input: bulk phase index || output: list of side ordinals of the master IP cells used for ghost facet construction
        Cell< Cell< moris_index > > mSlaveSideIgCellSideOrds;  // input: bulk phase index || output: list of side ordinals of the corresponding slave IP cells used for ghost facet construction

        // trivial clusters are ones that are not on a hanging node interface
        Cell< Cell< moris_index > > mTrivialFlag;        // input: bulk phase index || output: list of trivial flags indicating whether corresponding ghost facet construction is trivial (i.e. no hanging nodes)
        Cell< Cell< moris_index > > mTransitionLocation; // input: bulk phase index || output: list of transition locations of corresponding ghost facet construction 

        // ----------------------------------------------------------------------------------
        // Revised Maps for new Ghost

        // interpolation cells
        Cell< Cell< Cell< moris_index > > > mMasterSideIpCellsNew; // input: B-spline mesh list index, bulk phase index || output: list of enriched master IP cells used for ghost facet construction 
        Cell< Cell< Cell< moris_index > > > mSlaveSideIpCellsNew;  // input: B-spline mesh list index, bulk phase index || output: list of corresponding enriched slave IP cells used for ghost facet construction  

        // Side ordinals
        Cell< Cell< Cell< moris_index > > > mMasterSideIgCellSideOrdsNew; // input: B-spline mesh list index, bulk phase index || output: list of side ordinals of the master IP cells used for ghost facet construction
        Cell< Cell< Cell< moris_index > > > mSlaveSideIgCellSideOrdsNew;  // input: B-spline mesh list index, bulk phase index || output: list of side ordinals of the corresponding slave IP cells used for ghost facet construction

        // trivial clusters are ones that are not on a hanging node interface
        Cell< Cell< Cell< moris_index > > > mTrivialFlagNew;        // input: B-spline mesh list index, bulk phase index || output: list of trivial flags indicating whether corresponding ghost facet construction is trivial (i.e. no hanging nodes)
        Cell< Cell< Cell< moris_index > > > mTransitionLocationNew; // input: B-spline mesh list index, bulk phase index || output: list of transition locations of corresponding ghost facet construction 

        // ----------------------------------------------------------------------------------

        // Linear integration cells
        bool                                           mLinearBackgroundMesh = false;
        std::unordered_map< moris_index, moris_index > mLinearIgCellIndex;
        Cell< mtk::Cell* >                             mLinearIgCells;

        size_t
        capacity()
        {
            size_t tTotal = 0;
            tTotal += mSubphaseIndexToInterpolationCellIndex.capacity();
            tTotal += mDblSideSetIndexInMesh.capacity();
            tTotal += moris::internal_capacity( mMasterSideIgCellSideOrds );
            tTotal += moris::internal_capacity( mMasterSideIgCellSideOrds );
            tTotal += moris::internal_capacity( mSlaveSideIgCellSideOrds );
            tTotal += moris::internal_capacity( mTrivialFlag );
            tTotal += moris::internal_capacity( mTransitionLocation );
            return tTotal;
        }
    };

    // ----------------------------------------------------------------------------------

    class Ghost_Stabilization
    {
        private:

            // Pointer to the XTK-model
            Model* mXTKModel;

            // B-spline mesh information
            Matrix< IndexMat > mMeshIndices;
            moris_index mMinMeshIndex; 
            moris::Cell< Bspline_Mesh_Info* > mBsplineMeshInfos;

            bool
            is_linear_ip_mesh();

        public:
            // ----------------------------------------------------------------------------------
            /*!
            *  @brief Default constructor
            */
            Ghost_Stabilization();

            // ----------------------------------------------------------------------------------
            /*!
            *  @brief Constructor
            */
            Ghost_Stabilization( Model* aXTKModel );

            // ----------------------------------------------------------------------------------
            /*!
            *  @brief Perform/setup ghost stabilization
            */
            void
            setup_ghost_stabilization();

            void
            setup_ghost_stabilization_new();

            // ----------------------------------------------------------------------------------
            /*!
            *  @brief Work with meshes, to add visualization set for provided bulk phase
            *  @param[in] aBulkPhase Bulk phase
            */
            void
            visualize_ghost_on_mesh( moris_index const& aBulkPhase );

            void
            visualize_ghost_on_mesh_new(
                    moris_index const& aBsplineMeshListIndex,
                    moris_index const& aBulkPhaseIndex );

            // ----------------------------------------------------------------------------------
            /*!
            *  @brief Get the dbl side set name of ghost stabilization
            *  @param[in] aBulkPhase Bulk phase
            *  @return Name of double side set
            */
            std::string
            get_ghost_dbl_side_set_name(moris_index const & aBulkPhase);

            std::string
            get_ghost_dbl_side_set_name( 
                    moris_index const& aBulkPhase,
                    moris_index const& aBspMeshIndex );

            // ----------------------------------------------------------------------------------
            /*!
            * @return Memory map of ghost stabilization
            */
           Memory_Map
           get_memory_usage();

        private:
            // ----------------------------------------------------------------------------------
            /**
             * @brief Construct interpolation cell for trivial cell cluster built on child mesh subphases.
             * Find the non-trivial clusters and mark the underlying enriched IP elements for 
             * construction of new IP cells for ghost; and give the new enr. IP elems. IDs
             * This is done so that we can still have 1 cluster per interpolation cell.
             * 
             * Also: fill Ghost Setup Data with Subphase-Index-to-Enriched-IP-cell-index map 
             * 
             */
            void
            construct_ip_ig_cells_for_ghost_side_clusters(Ghost_Setup_Data & aGhostSetupData);

            // ----------------------------------------------------------------------------------
            /**
             * @brief identifies the vertices in the aura that will interpolate into a ghost facet.
             * In general, these vertices do not have a vertex interpolation so the vertex
             * interpolation information needs to be communicated. This function handles that.
             */
            void
            identify_and_setup_aura_vertices_in_ghost(Ghost_Setup_Data &  aGhostSetupData);

            // ----------------------------------------------------------------------------------

            void
            get_ip_vertices_in_ghost_sets(
                    Ghost_Setup_Data                & aGhostSetupData,
                    moris::Cell<mtk::Vertex*>       & aGhostVerticesWithInterpolation,
                    moris::Cell<mtk::Vertex*>       & aGhostVerticesWithoutInterpolation,
                    moris::Cell<mtk::Cell  const *> & aGhostIpCellConnectedToVertex);

            // ----------------------------------------------------------------------------------

            void
            prepare_interpolation_vertex_t_matrix_requests(
                    moris::Cell<mtk::Vertex*>             & aGhostVerticesWithoutInterpolation,
                    moris::Cell<mtk::Cell  const *>       & aGhostIpCellConnectedToVertex,
                    Cell<Matrix<IndexMat>>                & aNotOwnedIPVertIndsToProcs,
                    Cell<Matrix<IndexMat>>                & aNotOwnedBGIPVertsIdsToProcs,
                    Cell<Matrix<IndexMat>>                & aNotOwnedIpCellIdToProcs,
                    Cell<Matrix<IndexMat>>                & aNotOwnedEnrichedCellBulkPhaseToProcs,
                    Cell<uint>                            & aProcRanks,
                    std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex);
            // ----------------------------------------------------------------------------------
            /**
             * @brief Using a background vertex id, and a enriched interpolation cell id
             * find the corresponding enriched interpolation vertex
             */
            moris_index
            get_enriched_interpolation_vertex(moris_index const & aBGVertId,
                                              moris_index const & aEnrichedIpCellIndex);

            // ----------------------------------------------------------------------------------
            /**
             * @brief Packages t-matrix  request into mpi ready data structure (i.e. sparse 
             * matrix for returning back to requesting processor.
             */
            void
            prepare_t_matrix_request_answers(
                    moris_index            const & aMeshIndex,
                    Cell<Matrix<IndexMat>> const & aRequestedBgVertexIds,
                    Cell<Matrix<IndexMat>> const & aRequestedIpCellIds,
                    Cell<Matrix<IndexMat>> const & aIpCellBulkPhases,
                    Cell<Matrix<DDRMat>>         & aTMatrixWeights,
                    Cell<Matrix<IndexMat>>       & aTMatrixIndices,
                    Cell<Matrix<IndexMat>>       & aBasisOwners,
                    Cell<Matrix<IndexMat>>       & aTMatrixOffsets);

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
                    moris_index            const & aMeshIndex,
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

            bool
            verify_ghost_communication(Ghost_Setup_Data & aGhostSetupData);

            // ----------------------------------------------------------------------------------

            void
            declare_ghost_double_side_sets_in_mesh( Ghost_Setup_Data& aGhostSetupData );

            void
            declare_ghost_double_side_sets_in_mesh_new( Ghost_Setup_Data& aGhostSetupData );

            // ----------------------------------------------------------------------------------

            void
            construct_ghost_double_side_sets_in_mesh( Ghost_Setup_Data& aGhostSetupData );

            void
            construct_ghost_double_side_sets_in_mesh_new( Ghost_Setup_Data& aGhostSetupData );

            // ----------------------------------------------------------------------------------
            /**
             * @brief figure out whether to trigger Ghost side set creation between two Subphases
             * Rules:
             * 1. Only create ghost facets between a subphase created inside an intersected
             *    cell and its neighbors.
             * 2. The owning processor of the master (first) subphase constructs the ghost facet.
             * 3. Construct from coarse to fine in HMR
             * 
             * @param aGhostSetupData 
             * @param aFirstSubphase 
             * @param aSecondSubphase 
             * @param aTrivialFlag output: flag whether the transition between the master and slave element is trivial 
             * @return bool whether Ghost side set should be created or not
             */
            bool
            create_ghost(
                    Ghost_Setup_Data &  aGhostSetupData,
                    moris_index const & aFirstSubphase,
                    moris_index const & aSecondSubphase,
                    moris_index &       aTrivialFlag);

            bool
            create_ghost_new(
                    moris_index const& aBspMeshListIndex,
                    moris_index const& aFirstSpgIndex, 
                    moris_index const& aSecondSpgIndex,
                    moris_index&       aTrivialFlag);

            // ----------------------------------------------------------------------------------

            std::shared_ptr<Side_Cluster>
            create_slave_side_cluster(
                    Ghost_Setup_Data &  aGhostSetupData,
                    Cell<Interpolation_Cell_Unzipped*> & aEnrIpCells,
                    uint const   & aBulkIndex,
                    uint const   & aCellIndex,
                    moris_index  & aCurrentIndex,
                    moris_index  & aCurrentId);

            std::shared_ptr< Side_Cluster >
            create_slave_side_cluster_new(
                    Ghost_Setup_Data&                     aGhostSetupData,
                    Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
                    uint const&                           aBsplineMeshListIndex,
                    uint const&                           aBulkPhaseIndex,
                    uint const&                           aGhostFacetIndexInSet,
                    moris_index&                          aCurrentIndex,
                    moris_index&                          aCurrentId );

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

            std::shared_ptr< Side_Cluster >
            create_master_side_cluster_new(
                    Ghost_Setup_Data&                     aGhostSetupData,
                    Cell< Interpolation_Cell_Unzipped* >& aEnrIpCells,
                    uint const&                           aBsplineMeshListIndex,
                    uint const&                           aBulkPhaseIndex,
                    uint const&                           aGhostFacetIndexInSet,
                    Side_Cluster*                         aSlaveSideCluster,
                    moris_index&                          aCurrentIndex,
                    moris_index&                          aCurrentId );

            // ----------------------------------------------------------------------------------

            std::shared_ptr<xtk::Cell_XTK_No_CM>
            create_non_trivial_master_ig_cell(
                    Ghost_Setup_Data &  aGhostSetupData,
                    uint const & aBulkIndex,
                    uint const & aCellIndex,
                    Side_Cluster* aMasterSideCluster,
                    Side_Cluster* aSlaveSideCluster,
                    moris_index & aCurrentIndex,
                    moris_index & aCurrentId);


            // ----------------------------------------------------------------------------------

            mtk::Cell*
            get_linear_ig_cell( Ghost_Setup_Data                  & aGhostSetupData,
                                Interpolation_Cell_Unzipped       * aInterpCell,
                                moris_index                       & aCurrentIndex,
                                moris_index                       & aCurrentId );

            // ----------------------------------------------------------------------------------

            mtk::Cell*
            create_linear_ig_cell( Ghost_Setup_Data                  & aGhostSetupData,
                                   Interpolation_Cell_Unzipped const * aInterpCell,
                                   moris_index                       & aCurrentIndex,
                                   moris_index                       & aCurrentId );

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
    };

}



#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_ */
