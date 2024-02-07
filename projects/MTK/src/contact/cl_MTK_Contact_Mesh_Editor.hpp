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
#include "cl_MTK_QuadraturePointMapper_Ray.hpp"
#include "cl_MTK_Integrator.hpp"
#include "fn_assert.hpp"

namespace moris::mtk
{
    class Contact_Mesh_Editor
    {
      public:    // member functions
        Contact_Mesh_Editor(
                Integration_Mesh_DataBase_IG                          *aIGMesh,
                Interpolation_Mesh_DataBase_IP                        *aIPMesh,
                Integrator                                            &aIntegrator,
                Vector< Side_Set const * >                            &aCandidateSideSet,
                Vector< std::pair< moris_index, moris_index > > const &aCandidatePairs )
                : mIGMesh( aIGMesh )
                , mIPMesh( aIPMesh )
                , mIntegrator( std::move( aIntegrator ) )
                , mSideSets( aCandidateSideSet )
                , mCandidatePairs( aCandidatePairs )
                , mPointMapper( QuadraturePointMapper_Ray( aIGMesh, aCandidateSideSet, aCandidatePairs ) ){};

        void update_nonconformal_side_sets() const;

        void update_displacements( std::map< moris_index, Vector< real > > const &aNodalDisplacements );

        Vector< Side_Set const * > get_side_sets() const;

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
        static std::map< ClusterPair, std::map< CellPair, ResultIndices > > extract_cluster_and_cell_pairing( MappingResult const &aResult );

        /**
         * \brief Performs the mapping (i.e. ray tracing) from the follower side to the leader side. The follower side will also be called the source side, while the leader side (where the ray hits) will be called the target side. The mapping will be performed for all source sides of the candidate pairs.
         * \return Returns a mapping result for each source side of the candidate pairs.
         */
        Vector< MappingResult > perform_mapping() const;

        std::map< SetPair, Vector< Nonconformal_Side_Cluster > >
        convert_mapping_result_to_nonconformal_side_clusters( MappingResult const &aMappingResult ) const;

        void update_ig_mesh_database( const Vector< Nonconformal_Side_Cluster > &aNonconformalSideClusters, std::string const &aSetName, const Matrix< IndexMat > &aSetColor ) const;

        std::string get_nonconformal_side_set_name( SetPair const &tSetPair ) const;

        static std::pair< std::string, std::string > get_leaderphase_from_set_name( std::string const &aSideSetName );

        Integration_Mesh_DataBase_IG                   *mIGMesh;
        Interpolation_Mesh_DataBase_IP                 *mIPMesh;
        Integrator                                      mIntegrator;
        Vector< Side_Set const * >                      mSideSets;
        Vector< std::pair< moris_index, moris_index > > mCandidatePairs;
        QuadraturePointMapper_Ray                       mPointMapper;
        int mIteration = 0;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP
