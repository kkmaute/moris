/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Contact_Mesh_Editor.hpp
 *
 */

#ifndef MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP
#define MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP


#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Integrator.hpp"
#include "fn_assert.hpp"

#if MORIS_HAVE_ARBORX == ON
#include "cl_MTK_QuadraturePointMapper_Ray_ArborX.hpp"
using PointMapper = moris::mtk::QuadraturePointMapper_ArborX;
#else
// This will let the code compile without ArborX. But it will throw an error if nonconformal side sets are used.
#include "cl_MTK_QuadraturePointMapper_Ray_Dummy.hpp"
using PointMapper = moris::mtk::QuadraturePointMapper_Ray_Dummy;
#endif

namespace moris::mtk
{
    class Contact_Mesh_Editor
    {
      public:    // member functions
        Contact_Mesh_Editor(
                Integration_Mesh_DataBase_IG                          *aIGMesh,
                Integrator                                            &aIntegrator,
                Vector< Side_Set const * >                            &aCandidateSideSet,
                Vector< std::pair< moris_index, moris_index > > const &aCandidatePairs )
                : mIGMesh( aIGMesh )
                , mIntegrator( std::move( aIntegrator ) )
                , mSideSets( aCandidateSideSet )
                , mCandidatePairs( aCandidatePairs )
                , mPointMapper( PointMapper( aIGMesh, aCandidateSideSet, aCandidatePairs ) )
        {
        }

        void update_nonconformal_side_sets() const;

        void update_displacements( std::unordered_map< moris_index, Vector< real > > const &aNodalDisplacements );

        Vector< Side_Set const * > get_side_sets() const;

        void set_max_ray_length( real aMaxNegativeRayLength, real aMaxPositiveRayLength )
        {
            mMaxNegativeRayLength = aMaxNegativeRayLength;
            mMaxPositiveRayLength = aMaxPositiveRayLength;
        }

      private:
        // types
        /**
         * \brief Holds information about a cluster pair consisting of a source cluster and a target cluster on a target mesh.
         * \details This is used for the conversion from a mapping result to a nonconformal side cluster.
         * The first index is the source cluster index, the second index is the target cluster index and the third index is the target mesh index.
         */
        using ClusterPair = std::tuple< moris_index, moris_index, moris_index >;

        /**
         * \brief Holds information about a cell pair consisting of a source cell and a target cell.
         * \details This is used for the conversion from a mapping result to a nonconformal side cluster.
         * The first index is the source cell index, the second index is the target cell index.
         */
        using CellPair = std::pair< moris_index, moris_index >;

        /**
         * \brief Holds information about the indices in the mapping result that belong to a specific cluster- and cell-pair.
         */
        using ResultIndices = Vector< moris_index >;

        /**
         * \brief Holds information about a set of source- and target cells that are mapped to each other.
         */
        using SetPair = std::pair< moris_index, moris_index >;

        /**
         * \brief This method extracts the cluster and cell pairing from the mapping result and returns it in a structured way.
         * \details The (rather/quite) complex data structure of the return value is required to store the correct indices of the mapping results for
         * each unique pair of clusters (temporarily). Each unique pair of clusters (ClusterPair) is identified by a tuple of the source cluster index,
         * the target cluster index and the target side set index (since the target cluster is not necessarily the same set as the source cluster).
         * As value, a map is stored that maps the source- and target cell pairs (CellPair) to the list of indices to get access to the correct entries
         * of the mapping result (ResultIndices).
         *
         * (SourceCluster, TargetCluster, TargetMesh) -- n -> (SourceCell, TargetCell) -- n -> (Index in MappingResult, e.g. to access integration points)
         */
        static std::map< Contact_Mesh_Editor::CellPair, Contact_Mesh_Editor::ResultIndices > extract_cell_pairing( MappingResult const &aMappingResult, Vector< moris_index > aResultIndices );

        std::unordered_map< moris_index, std::map< std::pair< moris_index, moris_index >, Contact_Mesh_Editor::ResultIndices > > extract_cluster_pairing( MappingResult const &aMappingResult ) const;
        /**
         * \brief Performs the mapping (i.e. ray tracing) from the follower side to the leader side. The follower side will also be called the source side,
         * while the leader side (where the ray hits) will be called the target side. The mapping will be performed for all source sides of the candidate pairs.
         * \return Returns a mapping result for each source side of the candidate pairs.
         */
        Vector< MappingResult > perform_mapping( Matrix< DDRMat > aPointsToMap ) const;

        std::map< SetPair, Vector< Nonconformal_Side_Cluster > >
        convert_mapping_result_to_nonconformal_side_clusters( MappingResult const &aMappingResult ) const;

        std::string get_nonconformal_side_set_name( SetPair const &tSetPair ) const;

        static std::pair< std::string, std::string > get_leaderphase_from_set_name( std::string const &aSideSetName );

        /**
         * @brief The mapper has to map the nodes as well as the integration points. This function returns a combined matrix of the nodal parametric coordinates
         * (e.g. -1.0 and 1.0 for a line) as the first n_node columns and the integration points as the last n_ig columns.
         * @return
         */
        Matrix< DDRMat > get_points_to_map() const;

        /**
         * @brief The mapping does not know if a point is a nodal point or an integration point. It will return a continuous list of results that contains all requested
         * points for each cell in order. This method is used to determine if a result index belongs to a requested nodal coordinate or one of the integration points.
         * E.g. for a line element with 3 Integration points, the list of requested points is [ n1, n2, q1, q2, q3 ]. Thus the results would be ordered (for each source cell c):
         *
         * @verbatim
         *|--------------- c1 ---------------| |-------------- c2 ---------------|
         *[ c1_n1, c1_n2, c1_q1, c1_q2, c1_q3, c2_n1, c2_n2, c2_q1, c2_q2, c2_q3 ]
         *    0      1      2      3      4      5      6      7      8      9
         * @endverbatim
         *
         * This function would return false for indices 0, 1, 5 and 6 and true for 2, 3, 4, 7, 8 and 9.
         * @param aMappingResultColumnIndex
         * @return
         */
        bool is_integration_point_result_index( moris_index aMappingResultColumnIndex ) const;

        Matrix< DDRMat > get_nodal_parametric_coordinates() const;

        moris_index get_integration_point_index( moris_index aResultIndex ) const;

        moris_index get_node_coordinate_index( moris_index aMappingResultColumnIndex ) const;

        IntegrationPointPairs create_integration_point_pairs_from_results( Vector< moris_index > aResultIndices, MappingResult aMappingResult ) const;

        NodalPointPairs create_nodal_point_pairs_from_results( Vector< moris_index > aResultIndices, MappingResult aMappingResult ) const;

        void populate_integration_and_nodal_point_pairs(
                MappingResult const                      &aMappingResult,
                Vector< IntegrationPointPairs >          &aIntegrationPointPairs,
                Vector< NodalPointPairs >                &aNodePointPairs,
                Contact_Mesh_Editor::ResultIndices const &aCellResults ) const;

        Integration_Mesh_DataBase_IG                   *mIGMesh;
        Integrator                                      mIntegrator;
        Vector< Side_Set const * >                      mSideSets;
        Vector< std::pair< moris_index, moris_index > > mCandidatePairs;
        PointMapper                                     mPointMapper;
        real                                            mMaxNegativeRayLength;
        real                                            mMaxPositiveRayLength;
        int                                             mIteration = 0;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP
