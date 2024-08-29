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

namespace moris::mtk
{
    class Cell;
    }

namespace moris::xtk
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
        Vector< moris_index > mSubphaseIndexToInterpolationCellIndex;    // input: subphase index || output: enriched interpolation cell index
        Vector< moris_index > mDblSideSetIndexInMesh;                    // input: index of ghost side set || output: corresponding index of dbl. side set as saved on enriched IG mesh

        // ----------------------------------------------------------------------------------
        // Maps for old Ghost

        // interpolation cells
        Vector< Vector< moris_index > > mLeaderSideIpCells;      // input: bulk phase index || output: list of enriched leader IP cells used for ghost facet construction
        Vector< Vector< moris_index > > mFollowerSideIpCells;    // input: bulk phase index || output: list of corresponding enriched follower IP cells used for ghost facet construction

        // Side ordinals
        Vector< Vector< moris_index > > mLeaderSideIgCellSideOrds;      // input: bulk phase index || output: list of side ordinals of the leader IP cells used for ghost facet construction
        Vector< Vector< moris_index > > mFollowerSideIgCellSideOrds;    // input: bulk phase index || output: list of side ordinals of the corresponding follower IP cells used for ghost facet construction

        // trivial clusters are ones that are not on a hanging node interface
        Vector< Vector< moris_index > > mNonTrivialFlag;        // input: bulk phase index || output: list of flags indicating whether corresponding ghost facet construction is Non-trivial (i.e. it has hanging nodes)
        Vector< Vector< moris_index > > mTransitionLocation;    // input: bulk phase index || output: list of transition locations of corresponding ghost facet construction

        // ----------------------------------------------------------------------------------
        // Revised Maps for new Ghost

        // interpolation cells
        Vector< Vector< Vector< moris_index > > > mLeaderSideIpCellsNew;      // input: B-spline mesh list index, bulk phase index || output: list of enriched leader IP cells used for ghost facet construction
        Vector< Vector< Vector< moris_index > > > mFollowerSideIpCellsNew;    // input: B-spline mesh list index, bulk phase index || output: list of corresponding enriched follower IP cells used for ghost facet construction

        // Side ordinals
        Vector< Vector< Vector< moris_index > > > mLeaderSideIgCellSideOrdsNew;      // input: B-spline mesh list index, bulk phase index || output: list of side ordinals of the leader IP cells used for ghost facet construction
        Vector< Vector< Vector< moris_index > > > mFollowerSideIgCellSideOrdsNew;    // input: B-spline mesh list index, bulk phase index || output: list of side ordinals of the corresponding follower IP cells used for ghost facet construction

        // trivial clusters are ones that are not on a hanging node interface
        Vector< Vector< Vector< moris_index > > > mNonTrivialFlagNew;        // input: B-spline mesh list index, bulk phase index || output: list of trivial flags indicating whether corresponding ghost facet construction is trivial (i.e. no hanging nodes)
        Vector< Vector< Vector< moris_index > > > mTransitionLocationNew;    // input: B-spline mesh list index, bulk phase index || output: list of transition locations of corresponding ghost facet construction

        // ----------------------------------------------------------------------------------

        // Linear integration cells
        bool                                           mLinearBackgroundMesh = false;
        std::unordered_map< moris_index, moris_index > mLinearIgCellIndex;
        Vector< mtk::Cell* >                           mLinearIgCells;

        size_t
        capacity()
        {
            size_t tTotal = 0;
            tTotal += mSubphaseIndexToInterpolationCellIndex.capacity();
            tTotal += mDblSideSetIndexInMesh.capacity();
            tTotal += moris::internal_capacity( mLeaderSideIgCellSideOrds );
            tTotal += moris::internal_capacity( mLeaderSideIgCellSideOrds );
            tTotal += moris::internal_capacity( mFollowerSideIgCellSideOrds );
            tTotal += moris::internal_capacity( mNonTrivialFlag );
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
        Matrix< IndexMat >           mMeshIndices;
        moris_index                  mMinMeshIndex;
        Vector< Bspline_Mesh_Info* > mBsplineMeshInfos;

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
        visualize_ghost_on_mesh( moris_index const & aBulkPhase );

        void
        visualize_ghost_on_mesh_new(
                moris_index const & aBsplineMeshListIndex,
                moris_index const & aBulkPhaseIndex );

        // ----------------------------------------------------------------------------------
        /*!
         *  @brief Get the dbl side set name of ghost stabilization
         *  @param[in] aBulkPhase Bulk phase
         *  @return Name of double side set
         */
        std::string
        get_ghost_dbl_side_set_name( moris_index const & aBulkPhase );

        std::string
        get_ghost_dbl_side_set_name(
                moris_index const & aBulkPhase,
                moris_index const & aBspMeshIndex );

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
        construct_ip_ig_cells_for_ghost_side_clusters( Ghost_Setup_Data& aGhostSetupData );

        // ----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        /**
         * @brief identifies the vertices in the aura that will interpolate into a ghost facet.
         * In general, these vertices do not have a vertex interpolation so the vertex
         * interpolation information needs to be communicated. This function handles that.
         *
         * @param aGhostSetupData
         */
        void
        identify_and_setup_aura_vertices_in_ghost( Ghost_Setup_Data& aGhostSetupData );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Find for which unzipped IP vertices the executing processor doesn't have a T-matrix for a given B-spline mesh
         * This happens as basis information is only known up to the processor boundaries, but not into the aura
         *
         * @param aGhostSetupData data which is used to construct ghost facets
         * @param aGhostVerticesWithoutInterpolation output: list of unzipped IP vertices for which a T-matrix is missing
         * @param aGhostIpCellConnectedToVertex output: unzipped IP cells that these IP vertices belong to
         * @param aMeshIndex input: B-spline mesh index for which the UIPVs with missing T-matrices are to be collected
         */
        void
        collect_unzipped_IP_vertices_without_interpolation(
                Ghost_Setup_Data&            aGhostSetupData,
                Vector< mtk::Vertex* >&      aGhostVerticesWithoutInterpolation,
                Vector< mtk::Cell const * >& aGhostIpCellConnectedToVertex,
                const moris_index            aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Prepare information by which a UIPV can be identified by another processor
         *
         * @param aGhostVerticesWithoutInterpolation input: list of unzipped IP vertices for which the T-matrix is known to be missing
         * @param aGhostIpCellConnectedToVertex input: unzipped IP cells that these IP vertices are attached to
         * @param aNotOwnedIPVertIndsInNotOwnedList output: outer cell: owner proc index in XTK comm-table || inner cell: position of the IP vertex (for which T-matrix information should be communicated) in the list of all UIPVs with missing T-matrices
         * @param aIPVertIndsToProcs output: outer cell: owner proc index in XTK comm-table || inner cell: index of the UIPV with a missing T-matrix
         * @param aBaseVertexIds output: outer cell: owner proc index in XTK comm-table || inner cell: ID of the base vertex underlying the UIPV
         * @param aUnzippedIpCellIds output: outer cell: owner proc index in XTK comm-table || inner cell: ID of the UIPC the UIPV is attached to
         */
        void
        prepare_requests_for_T_matrices_without_interpolation(
                Vector< mtk::Vertex* > const &      aGhostVerticesWithoutInterpolation,
                Vector< mtk::Cell const * > const & aGhostIpCellConnectedToVertex,
                Vector< Vector< moris_index > >&    aNotOwnedIPVertIndsInNotOwnedList,
                Vector< Vector< moris_index > >&    aIPVertIndsToProcs,
                Vector< Matrix< IdMat > >&          aBaseVertexIds,
                Vector< Matrix< IdMat > >&          aUnzippedIpCellIds );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Find UIPVs based on received information and send T-matrices for this vertices wrt to a given B-spline mesh
         *
         * @param aReceivedBaseVertexIds input: outer cell: owner proc index in XTK comm-table || inner cell: ID of the base vertex underlying the UIPV
         * @param aReceivedUnzippedIpCellIds input: outer cell: owner proc index in XTK comm-table || inner cell: ID of the UIPC the UIPV is attached to
         * @param aTMatrixWeights output: outer cell: owner proc index in XTK comm-table || inner cell: linear array of all T-matrix weights for all UIPVs received
         * @param aTMatrixIds output: outer cell: owner proc index in XTK comm-table || inner cell: linear array of all T-matrix basis IDs for all UIPVs received
         * @param aTMatrixOwners output: outer cell: owner proc index in XTK comm-table || inner cell: linear array of all T-matrix basis owners for all UIPVs received
         * @param aTMatrixOffsets output: outer cell: owner proc index in XTK comm-table || inner cell: start index of the T-matrix info for a given UIPV in the linear arrays above
         * @param aMeshIndex input: B-spline mesh index for which the requested UIPVs need the T-matrix
         */
        void
        prepare_answers_for_T_matrices(
                Vector< Matrix< IdMat > > const & aReceivedBaseVertexIds,
                Vector< Matrix< IdMat > > const & aReceivedUnzippedIpCellIds,
                Vector< Matrix< DDRMat > >&       aTMatrixWeights,
                Vector< Matrix< IdMat > >&        aTMatrixIds,
                Vector< Matrix< IdMat > >&        aTMatrixOwners,
                Vector< Matrix< IndexMat > >&     aTMatrixOffsets,
                const moris_index                 aMeshIndex );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Reconstruct T-matrices from the linear T-matrix arrays received and assign them to the UIPVs with missing T-matrices
         *
         * @param aGhostIpCellConnectedToVertex  unzipped IP cells that the unzipped IP vertices are attached to
         * @param aNotOwnedIPVertIndsInNotOwnedList outer cell: owner proc index in XTK comm-table || inner cell: position of the IP vertex (for which T-matrix information is received) in the list of all UIPVs with missing T-matrices
         * @param aIPVertIndsToProcs outer cell: owner proc index in XTK comm-table || inner cell: index of the UIPV with a missing T-matrix
         * @param aReceivedTMatrixWeights outer cell: owner proc index in XTK comm-table || inner cell: linear array of all T-matrix weights for all UIPVs requested
         * @param aReceivedTMatrixIds outer cell: owner proc index in XTK comm-table || inner cell: linear array of all T-matrix basis IDs for all UIPVs requested
         * @param aReceivedTMatrixOwners outer cell: owner proc index in XTK comm-table || inner cell: linear array of all T-matrix basis owners for all UIPVs requested
         * @param aReceivedTMatrixOffsets outer cell: owner proc index in XTK comm-table || inner cell: start index of the T-matrix info for a given UIPV in the linear arrays above
         * @param aMeshIndex B-spline mesh index for which the T-matrix is received
         */
        void
        handle_requested_T_matrix_answers(
                Vector< mtk::Cell const * > const &     aGhostIpCellConnectedToVertex,
                Vector< Vector< moris_index > > const & aNotOwnedIPVertIndsInNotOwnedList,
                Vector< Vector< moris_index > > const & aIPVertIndsToProcs,
                Vector< Matrix< DDRMat > > const &      aReceivedTMatrixWeights,
                Vector< Matrix< IdMat > > const &       aReceivedTMatrixIds,
                Vector< Matrix< IdMat > > const &       aReceivedTMatrixOwners,
                Vector< Matrix< IndexMat > > const &    aReceivedTMatrixOffsets,
                const moris_index                       aMeshIndex );

        // ----------------------------------------------------------------------------------
        // ----------------------------------------------------------------------------------

        void
        identify_and_setup_aura_vertices_in_ghost_old( Ghost_Setup_Data& aGhostSetupData );

        // ----------------------------------------------------------------------------------

        void
        get_ip_vertices_in_ghost_sets(
                Ghost_Setup_Data&            aGhostSetupData,
                Vector< mtk::Vertex* >&      aGhostVerticesWithInterpolation,
                Vector< mtk::Vertex* >&      aGhostVerticesWithoutInterpolation,
                Vector< mtk::Cell const * >& aGhostIpCellConnectedToVertex );

        // ----------------------------------------------------------------------------------
        /**
         * @brief Using a background vertex id, and a enriched interpolation cell id
         * find the corresponding enriched interpolation vertex
         */
        moris_index
        get_enriched_interpolation_vertex( moris_index const & aBGVertId,
                moris_index const &                            aEnrichedIpCellIndex );

        // ----------------------------------------------------------------------------------

        void
        prepare_ip_cell_id_answers(
                Vector< Matrix< IndexMat > >&             aReceivedEnrCellIds,
                Vector< moris_id >&                       aNewInterpCellIds,
                Vector< Matrix< IndexMat > >&             aEnrCellIds,
                std::unordered_map< moris_id, moris_id >& aBaseEnrIdToIndexInNonTrivialOwned );

        // ----------------------------------------------------------------------------------

        bool
        verify_ghost_communication( Ghost_Setup_Data& aGhostSetupData );

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
         * 2. The owning processor of the leader (first) subphase constructs the ghost facet.
         * 3. Construct from coarse to fine in HMR
         *
         * @param aGhostSetupData
         * @param aFirstSubphase
         * @param aSecondSubphase
         * @param aTrivialFlag output: flag whether the transition between the leader and follower element is trivial
         * @return bool whether Ghost side set should be created or not
         */
        bool
        create_ghost(
                Ghost_Setup_Data&   aGhostSetupData,
                moris_index const & aFirstSubphase,
                moris_index const & aSecondSubphase,
                moris_index&        aTrivialFlag );

        bool
        create_ghost_new(
                moris_index const & aBspMeshListIndex,
                moris_index const & aFirstSpgIndex,
                moris_index const & aSecondSpgIndex,
                moris_index&        aNonTrivialFlag );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Side_Cluster >
        create_follower_side_cluster(
                Ghost_Setup_Data&                       aGhostSetupData,
                Vector< Interpolation_Cell_Unzipped* >& aEnrIpCells,
                uint const &                            aBulkIndex,
                uint const &                            aCellIndex,
                moris_index&                            aCurrentIndex,
                moris_index&                            aCurrentId );

        std::shared_ptr< Side_Cluster >
        create_follower_side_cluster_new(
                Ghost_Setup_Data&                       aGhostSetupData,
                Vector< Interpolation_Cell_Unzipped* >& aEnrIpCells,
                uint const &                            aBsplineMeshListIndex,
                uint const &                            aBulkPhaseIndex,
                uint const &                            aGhostFacetIndexInSet,
                moris_index&                            aCurrentIndex,
                moris_index&                            aCurrentId );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< Side_Cluster >
        create_leader_side_cluster(
                Ghost_Setup_Data&                       aGhostSetupData,
                Vector< Interpolation_Cell_Unzipped* >& aEnrIpCells,
                uint const &                            aBulkIndex,
                uint const &                            aCellIndex,
                Side_Cluster*                           aFollowerSideCluster,
                moris_index&                            aCurrentIndex,
                moris_index&                            aCurrentId );

        std::shared_ptr< Side_Cluster >
        create_leader_side_cluster_new(
                Ghost_Setup_Data&                       aGhostSetupData,
                Vector< Interpolation_Cell_Unzipped* >& aEnrIpCells,
                uint const &                            aBsplineMeshListIndex,
                uint const &                            aBulkPhaseIndex,
                uint const &                            aGhostFacetIndexInSet,
                Side_Cluster*                           aFollowerSideCluster,
                moris_index&                            aCurrentIndex,
                moris_index&                            aCurrentId );

        // ----------------------------------------------------------------------------------

        std::shared_ptr< xtk::Cell_XTK_No_CM >
        create_non_trivial_leader_ig_cell(
                Ghost_Setup_Data& aGhostSetupData,
                uint const &      aBulkIndex,
                uint const &      aCellIndex,
                Side_Cluster*     aLeaderSideCluster,
                Side_Cluster*     aFollowerSideCluster,
                moris_index&      aCurrentIndex,
                moris_index&      aCurrentId );

        std::shared_ptr< xtk::Cell_XTK_No_CM >
        create_non_trivial_leader_ig_cell_new(
                Ghost_Setup_Data& aGhostSetupData,
                uint const &      aBsplineMeshListIndex,
                uint const &      aBulkIndex,
                uint const &      aCellIndex,
                Side_Cluster*     aLeaderSideCluster,
                Side_Cluster*     aFollowerSideCluster,
                moris_index&      aCurrentIndex,
                moris_index&      aCurrentId );

        // ----------------------------------------------------------------------------------

        mtk::Cell*
        get_linear_ig_cell( Ghost_Setup_Data& aGhostSetupData,
                Interpolation_Cell_Unzipped*  aInterpCell,
                moris_index&                  aCurrentIndex,
                moris_index&                  aCurrentId );

        // ----------------------------------------------------------------------------------

        mtk::Cell*
        create_linear_ig_cell( Ghost_Setup_Data&    aGhostSetupData,
                Interpolation_Cell_Unzipped const * aInterpCell,
                moris_index&                        aCurrentIndex,
                moris_index&                        aCurrentId );

        // ----------------------------------------------------------------------------------

        void
        permute_follower_vertices(
                Vector< moris::mtk::Vertex const * > const & aFollowerVertices,
                Vector< moris::mtk::Vertex const * > const & aLeaderVertices,
                Vector< moris::mtk::Vertex const * >&        aPermutedFollowerVertices,
                Vector< moris::mtk::Vertex const * >&        aPermutedAdjMastVertices );

        // ----------------------------------------------------------------------------------

        void
        get_local_coords_on_transition_side( moris_index const & aMySideOrdinal,
                moris_index const &                              aTransitionLoc,
                Vector< Matrix< DDRMat > >&                      aLocCoord );

        // ----------------------------------------------------------------------------------

        moris_index get_side_ordinals_for_non_trivial_leader();
    };

}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_GHOST_STABILIZATION_HPP_ */
