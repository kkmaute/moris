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
                Integration_Mesh_DataBase_IG                               *aIGMesh,
                Interpolation_Mesh_DataBase_IP                             *aIPMesh,
                const Integrator                                           &aIntegrator,
                moris::Cell< Side_Set * >                                  &aCandidateSideSet,
                moris::Cell< std::pair< moris_index, moris_index > > const &aCandidatePairs )
                : mIGMesh( aIGMesh )
                , mIPMesh( aIPMesh )
                , mIntegrator( aIntegrator )
                , mSideSets( aCandidateSideSet )
                , mCandidatePairs( aCandidatePairs )
                , mPointMapper( QuadraturePointMapper_Ray( aIGMesh, aCandidateSideSet, aCandidatePairs ) ){};

        void update_nonconformal_side_sets();

        void update_displacements( Matrix< DDRMat > &aDisplacements );

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
        using ResultIndices = moris::Cell< moris_index >;

        // methods

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

        moris::Cell< Nonconformal_Side_Cluster > convert_mapping_result_to_nonconformal_side_clusters( moris_index aSourceSideSetIndex, MappingResult aResult );

        // data
        Integration_Mesh_DataBase_IG                        *mIGMesh;
        Interpolation_Mesh_DataBase_IP                      *mIPMesh;
        Integrator                                           mIntegrator;
        moris::Cell< Side_Set * >                            mSideSets;
        moris::Cell< std::pair< moris_index, moris_index > > mCandidatePairs;
        QuadraturePointMapper_Ray                            mPointMapper;
    };
}    // namespace moris::mtk
#endif    // MORIS_CL_MTK_CONTACT_MESH_EDITOR_HPP
