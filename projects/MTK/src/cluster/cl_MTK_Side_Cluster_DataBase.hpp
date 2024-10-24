/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster_DataBase.hpp
 *
 */
#ifndef SRC_cl_MTK_Side_Cluster_DataBase
#define SRC_cl_MTK_Side_Cluster_DataBase

#include "cl_MTK_Side_Cluster.hpp"

namespace moris::mtk
{
    class Cell_Cluster;
    class Mesh;

    class Side_Cluster_DataBase : public mtk::Side_Cluster
    {
      private:
        // FIXME: old data member that needs to be deleted
        Vector< moris::mtk::Cell const * > mPrimaryIntegrationCells;
        Vector< moris::mtk::Cell const * > mVoidIntegrationCells;

        moris_index mSideClusterIndex;
        mtk::Mesh*  mMesh = nullptr;

      public:
        // ----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Side_Cluster_DataBase
         *
         */
        Side_Cluster_DataBase() = default;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Side_Cluster_DataBase object
         *
         * @param aSideClusterIndex
         * @param aMesh
         */

        Side_Cluster_DataBase(
                moris_index aSideClusterIndex,
                mtk::Mesh*  aMesh );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Destroy the Side_Cluster_DataBase object
         *
         */

        ~Side_Cluster_DataBase() override = default;

        // ----------------------------------------------------------------------------------

        /*!
         * Indicates there is a 1 to 1 relationship between
         * integration cell and interpolation cells in this cluster
         */

        bool
        is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /*!
         * Get interpolation cell interpolating into this side cluster
         * @param[in] aIsLeader  Leader or Follower Selector Enum (for Double side clusters only)
         * @return Interpolation cell related to cluster
         */

        moris::mtk::Cell const &
        get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /*!
         * Get all integration cells in this side cluster
         */

        Vector< mtk::Cell const * > const &
        get_cells_in_side_cluster() const override;

        // ----------------------------------------------------------------------------------

        /*!
         * @param[in] aIsLeader  Leader or Follower Selector Enum (for Double side clusters only)
         * @return all integration cell side ordinals
         */

        moris::Matrix< moris::IndexMat >
        get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /*!
         * Single side ordinal version of above
         * @param[in] aCellIndexInCluster Integration cell cluster index
         * @return single integration cell side ordinal
         *
         */
        moris_index
        get_cell_side_ordinal(
                moris::moris_index         aCellIndexInCluster,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /*!
         * @return all the vertices in this cluster
         */

        Vector< moris::mtk::Vertex const * >
        get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------
        // Local Coordinate Access
        // ----------------------------------------------------------------------------------

        /*!
         * Access the full array of local coordinates
         * @param[in] aIsLeader - Leader or Follower Selector Enum (for Double side clusters only)
         * @return All vertex local coordinates
         */

        moris::Matrix< moris::DDRMat >
        get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /*!
         * Access a single local coordinate of a vertex
         * @param[in] aIsLeader - Leader or Follower Selector Enum (for Double side clusters only)
         * @return Single vertex local coordinates
         */
        moris::Matrix< moris::DDRMat >
        get_vertex_local_coordinate_wrt_interp_cell(
                moris::mtk::Vertex const * aVertex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /*!
         * @return Size of the xsi vector in this side cluster
         */

        moris_index
        get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /**
         * @brief populate the data that is being returned by refernce
         *
         */

        void
        set_outward_data();

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the vertex cluster index
         *
         * @param aVertex
         * @param aIsLeader
         * @return moris_index
         */

        moris_index
        get_vertex_cluster_index(
                const Vertex*              aVertex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const override;

        // ----------------------------------------------------------------------------------

        /**
         * @brief Get the vertex ordinal on facet
         *
         * @param aCellIndexInCluster
         * @param aVertex
         * @return moris_index
         */

        moris_index
        get_vertex_ordinal_on_facet(
                moris_index                aCellIndexInCluster,
                moris::mtk::Vertex const * aVertex ) const override;

        //---------------------------------------------------------------------------------------

        /**
         * @brief compute cluster cell measure
         *
         * @param aPrimaryOrVoid
         * @param aIsLeader
         * @return moris::real
         */

        moris::real
        compute_cluster_cell_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override;

        //---------------------------------------------------------------------------------------

        /**
         * @brief volume of the ig cells in matrix
         *
         * @param aPrimaryOrVoid
         * @param aIsLeader
         * @return Matrix< DDRMat >
         */

        Matrix< DDRMat >
        compute_cluster_ig_cell_measures(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override;

        //---------------------------------------------------------------------------------------

        /**
         * @brief compute cluster cell measure derivative
         *
         * @param aPerturbedVertexCoords
         * @param aDirection
         * @param aPrimaryOrVoid
         * @param aIsLeader
         * @return moris::real
         */

        moris::real
        compute_cluster_cell_measure_derivative(
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const override;

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( false, "mtk::Side_Cluster_DataBase::compute_cluster_group_cell_measure() - Not implemented in this child class." );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( false, "mtk::Side_Cluster_DataBase::compute_cluster_group_cell_measure_derivative() - Not implemented in this child class." );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_side_measure(
                const moris_index          aDiscretizationMeshIndex,
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster_DataBase::compute_cluster_group_cell_side_measure() - Not implemented in this child class." );
            return 0.0;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cluster_group_cell_side_measure_derivative(
                const moris_index          aDiscretizationMeshIndex,
                const Matrix< DDRMat >&    aPerturbedVertexCoords,
                uint                       aDirection,
                const mtk::Primary_Void    aPrimaryOrVoid,
                const mtk::Leader_Follower aIsLeader ) const override
        {
            MORIS_ERROR( false, "mtk::Cell_Cluster_DataBase::compute_cluster_group_cell_side_measure_derivative() - Not implemented in this child class." );
            return 0.0;
        }

        //---------------------------------------------------------------------------------------

        /**
         * @brief memory usage of the side cluster
         *
         * @return size_t
         */

        size_t
        capacity();

        //---------------------------------------------------------------------------------------

        /**
         * @brief
         *
         * @param aNewIndex
         */

        void
        update_cluster_index( moris_index aNewIndex );
    };

}    // namespace moris::mtk

#endif /* cl_MTK_Side_Cluster_DataBase.hpp */
