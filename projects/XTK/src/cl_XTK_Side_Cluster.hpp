/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Side_Cluster.hpp
 *
 */
#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_SIDE_CLUSTER_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_SIDE_CLUSTER_HPP_

#include "cl_MTK_Side_Cluster.hpp"
#include <unordered_map>

// forward declare the mtk::Cluster_Group
namespace moris::mtk
{
    class Cluster_Group;
}

using namespace moris;

namespace moris::xtk
{
    class Interpolation_Cell_Unzipped;
    class Cell_Cluster;

    struct IG_Vertex_Group;

    class Side_Cluster : public mtk::Side_Cluster
    {
        //---------------------------------------------------------------------------------------

        friend class Enriched_Integration_Mesh;
        friend class Ghost_Stabilization;

        //---------------------------------------------------------------------------------------

      protected:
        bool                                          mTrivial;
        Interpolation_Cell_Unzipped const *           mInterpolationCell;
        Vector< moris::mtk::Cell const * >            mIntegrationCells;
        Matrix< IndexMat >                            mIntegrationCellSideOrdinals;
        Vector< moris::mtk::Vertex const * >          mVerticesInCluster;
        std::shared_ptr< IG_Vertex_Group >            mVertexGroup;
        Vector< Matrix< DDRMat > >                    mVertexLocalCoords;
        Matrix< DDRMat >                              mVertexLocalCoordsMat;  /*FIXME: get rid of mVertexLocalCoords*/
        xtk::Cell_Cluster const *                     mAssociatedCellCluster; /* Associated cell cluster (needed for volume computations in Nitsche).*/
        Vector< std::weak_ptr< mtk::Cluster_Group > > mClusterGroups;

        //---------------------------------------------------------------------------------------

      public:
        //---------------------------------------------------------------------------------------

        Side_Cluster();

        //---------------------------------------------------------------------------------------

        ~Side_Cluster(){};

        //---------------------------------------------------------------------------------------

        bool is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::mtk::Cell const &
        get_interpolation_cell(
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        Vector< mtk::Cell const * > const &
        get_cells_in_side_cluster() const;

        //---------------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_cell_side_ordinals(
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris_index
        get_cell_side_ordinal(
                moris::moris_index         aCellIndexInCluster,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        Vector< moris::mtk::Vertex const * >
        get_vertices_in_cluster(
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        Matrix< DDRMat >
        get_vertices_local_coordinates_wrt_interp_cell(
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::moris_index
        get_vertex_cluster_index(
                const moris::mtk::Vertex*  aVertex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris_index
        get_vertex_ordinal_on_facet(
                moris_index                aCellIndexInCluster,
                moris::mtk::Vertex const * aVertex ) const;

        //---------------------------------------------------------------------------------------

        Matrix< DDRMat >
        get_vertex_local_coordinate_wrt_interp_cell(
                const moris::mtk::Vertex*  aVertex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris_index
        get_dim_of_param_coord(
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        Matrix< DDRMat >
        compute_cluster_ig_cell_measures(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_measure_derivative(
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const;

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_side_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const;

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_side_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const;

        //---------------------------------------------------------------------------------------

        void
        set_ig_vertex_group( std::shared_ptr< IG_Vertex_Group > aVertexGroup );

        //---------------------------------------------------------------------------------------

        bool
        has_cluster_group( const moris_index aDiscretizationMeshIndex ) const override;

        //---------------------------------------------------------------------------------------

        std::shared_ptr< mtk::Cluster_Group >
        get_cluster_group( const moris_index aDiscretizationMeshIndex ) const override;

        //---------------------------------------------------------------------------------------

        void
        set_cluster_group(
                const moris_index                     aDiscretizationMeshIndex,
                std::shared_ptr< mtk::Cluster_Group > aClusterGroupPtr ) override;

        //---------------------------------------------------------------------------------------

        // memory
        size_t
        capacity();

        //---------------------------------------------------------------------------------------

        void print_vertex_map() const;

        //---------------------------------------------------------------------------------------

    };    // class Side_Cluster

    //---------------------------------------------------------------------------------------

}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_SIDE_CLUSTER_HPP_ */
