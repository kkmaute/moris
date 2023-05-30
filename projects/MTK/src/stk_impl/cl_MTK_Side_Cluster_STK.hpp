/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster_STK.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_

#include <unordered_map>

#include "cl_MTK_Side_Cluster.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    namespace mtk
    {

        class Side_Cluster_STK : public Side_Cluster
        {
          private:
            bool                                      mTrivial;
            moris::mtk::Cell const                   *mInterpolationCell;
            moris::Cell< moris::mtk::Cell const * >   mIntegrationCells;
            moris::Matrix< moris::IndexMat >          mIntegrationCellSideOrdinals;
            moris::Cell< moris::mtk::Vertex const * > mVerticesInCluster;
            moris::Matrix< moris::DDRMat >            mVertexParamCoords;

            // map from vertex id to local index
            std::unordered_map< moris_index, moris_index > mVertexIdToLocalIndex;

          public:
            //----------------------------------------------------------------

            Side_Cluster_STK();

            //----------------------------------------------------------------

            // trivial constructor
            Side_Cluster_STK(
                    moris::mtk::Cell const                          *aInterpCell,
                    moris::mtk::Cell const                          *aIntegrationCell,
                    moris::Cell< moris::mtk::Vertex const * > const &aVerticesInCluster,
                    moris_index                                      aSideOrdinal );

            //----------------------------------------------------------------
            Side_Cluster_STK(
                    bool                                             aTrivial,
                    moris::mtk::Cell const                          *aInterpolationCell,
                    moris::Cell< moris::mtk::Cell const * > const   &aIntegrationCells,
                    moris::Matrix< moris::IndexMat > const          &aIntegrationCellSideOrdinals,
                    moris::Cell< moris::mtk::Vertex const * > const &aVerticesInCluster,
                    moris::Matrix< moris::DDRMat > const            &aVertexParamCoords );

            //----------------------------------------------------------------

            bool
            is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::mtk::Cell const &
            get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::Cell< moris::mtk::Cell const * > const &
            get_cells_in_side_cluster() const;

            //----------------------------------------------------------------

            moris::Matrix< moris::IndexMat >
            get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris_index
            get_cell_side_ordinal(
                    moris::moris_index      aCellIndexInCluster,
                    const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::Cell< moris::mtk::Vertex const * >
            get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris::Matrix< moris::DDRMat >
            get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris_index
            get_vertex_cluster_index(
                    moris::mtk::Vertex const *aVertex,
                    const mtk::Leader_Follower   aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris_index
            get_vertex_ordinal_on_facet(
                    moris_index               aCellIndexInCluster,
                    moris::mtk::Vertex const *aVertex ) const;

            //----------------------------------------------------------------

            moris::Matrix< moris::DDRMat >
            get_vertex_local_coordinate_wrt_interp_cell(
                    moris::mtk::Vertex const *aVertex,
                    const mtk::Leader_Follower   aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris_index
            get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

            //----------------------------------------------------------------

            moris_index
            get_vertex_cluster_local_index( moris_id aVertexId ) const;

            //----------------------------------------------------------------

            void
            add_vertex_to_map(
                    moris_id    aVertexId,
                    moris_index aVertexLocalIndex );

            //----------------------------------------------------------------
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_ */

