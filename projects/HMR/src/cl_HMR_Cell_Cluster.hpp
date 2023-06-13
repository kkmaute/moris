/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Cell_Cluster.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_CELL_CLUSTER_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_CELL_CLUSTER_HPP_

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_Matrix.hpp"
#include <unordered_map>

namespace moris::hmr
{

    class Cell_Cluster_HMR: public mtk::Cell_Cluster
    {
    private:
        bool                                    mTrivial;
        moris::mtk::Cell const *                mInterpolationCell;
        moris::Cell<moris::mtk::Cell const *>   mPrimaryIntegrationCells;
        moris::Cell<moris::mtk::Cell const *>   mVoidIntegrationCells;
        moris::Cell<moris::mtk::Vertex const *> mVerticesInCluster;
        moris::Matrix<moris::DDRMat>            mVertexParamCoords;

        // map from vertex id to local index
        std::unordered_map<moris_index,moris_index> mVertexIdToLocalIndex;           // FIXME should be ordered map. about 1000 times faster

    public:
        Cell_Cluster_HMR():
            mTrivial(true),
            mInterpolationCell(nullptr),
            mPrimaryIntegrationCells(0,nullptr),
            mVoidIntegrationCells(0,nullptr),
            mVerticesInCluster(0,nullptr),
            mVertexParamCoords(0,0)
        {};

        ~Cell_Cluster_HMR(){}

        //----------------------------------------------------------------
        bool
        is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //##############################################
        // Add and setup of cluster
        //##############################################
        void
        mark_as_nontrivial();

        //----------------------------------------------------------------

        void
        set_interpolation_cell(moris::mtk::Cell const * aInterpCell);

        //----------------------------------------------------------------

        void
        add_primary_integration_cell(moris::mtk::Cell  const * aIntegrationCell);

        //----------------------------------------------------------------

        void
        add_primary_integration_cell(moris::Cell<moris::mtk::Cell  const *> const & aIntegrationCell);

        //----------------------------------------------------------------
        void
        add_void_integration_cell(moris::Cell<moris::mtk::Cell const *> const & aIntegrationCell);

        //----------------------------------------------------------------

        void
        add_vertex_to_cluster(moris::Cell<moris::mtk::Vertex const *> const & aVertex);

        //----------------------------------------------------------------

        void
        add_vertex_local_coordinates_wrt_interp_cell(moris::Matrix<moris::DDRMat> const & aLocalCoords);

        //----------------------------------------------------------------

        //##############################################
        // Required Access Functions
        //##############################################

        moris::Cell<moris::mtk::Cell const *> const &
        get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const;

        //----------------------------------------------------------------

        moris::Cell<moris::mtk::Cell const *> const &
        get_void_cells_in_cluster() const;

        //----------------------------------------------------------------

        moris::mtk::Cell const & get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //----------------------------------------------------------------

        moris::Cell<moris::mtk::Vertex const *>
        get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //----------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const;

        //----------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                                                     const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const;

        //----------------------------------------------------------------

        moris_index
        get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //----------------------------------------------------------------
        //##############################################
        // Private access functions
        //##############################################
    private:
        moris_index
        get_vertex_cluster_local_index(moris_id aVertexId) const;

        //----------------------------------------------------------------

        void
        add_vertex_to_map(moris_id aVertexId,
                          moris_index aVertexLocalIndex);
    };
}

#endif /* PROJECTS_HMR_SRC_CL_HMR_CELL_CLUSTER_HPP_ */

