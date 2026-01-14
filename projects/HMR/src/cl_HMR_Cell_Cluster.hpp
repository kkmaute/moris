/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Cell_Cluster.hpp
 *
 */

#pragma once

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
        Vector<moris::mtk::Cell const *>   mPrimaryIntegrationCells;
        Vector<moris::mtk::Cell const *>   mVoidIntegrationCells;
        Vector<moris::mtk::Vertex const *> mVerticesInCluster;
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

        ~Cell_Cluster_HMR() override{}

        //----------------------------------------------------------------
        bool
        is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        //##############################################
        // Add and setup of cluster
        //##############################################
        void
        mark_as_nontrivial() override;

        //----------------------------------------------------------------

        void
        set_interpolation_cell(moris::mtk::Cell const * aInterpCell) override;

        //----------------------------------------------------------------

        void
        add_primary_integration_cell(moris::mtk::Cell  const * aIntegrationCell);

        //----------------------------------------------------------------

        void
        add_primary_integration_cell(Vector<moris::mtk::Cell  const *> const & aIntegrationCell) override;

        //----------------------------------------------------------------
        void
        add_void_integration_cell(Vector<moris::mtk::Cell const *> const & aIntegrationCell) override;

        //----------------------------------------------------------------

        void
        add_vertex_to_cluster(Vector<moris::mtk::Vertex const *> const & aVertex);

        //----------------------------------------------------------------

        void
        add_vertex_local_coordinates_wrt_interp_cell(moris::Matrix<moris::DDRMat> const & aLocalCoords);

        //----------------------------------------------------------------

        //##############################################
        // Required Access Functions
        //##############################################

        Vector<moris::mtk::Cell const *> const &
        get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const override;

        //----------------------------------------------------------------

        Vector<moris::mtk::Cell const *> const &
        get_void_cells_in_cluster() const override;

        //----------------------------------------------------------------

        moris::mtk::Cell const & get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        //----------------------------------------------------------------

        Vector<moris::mtk::Vertex const *>
        get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        //----------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const override;

        //----------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                                                     const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const override;

        //----------------------------------------------------------------

        moris_index
        get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

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
