/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Cluster_Proxy.hpp
 *
 */

#ifndef PROJECTS_MTK_TEST_CL_MTK_CELL_CLUSTER_PROXY_HPP_
#define PROJECTS_MTK_TEST_CL_MTK_CELL_CLUSTER_PROXY_HPP_

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
    class Cell_Cluster_Proxy : public mtk::Cell_Cluster
    {
      public:
        bool                                      mTrivial;
        moris::mtk::Cell                         *mInterpolationCell;
        Vector< moris::mtk::Cell const * >   mPrimaryIntegrationCells;
        Vector< moris::mtk::Cell const * >   mVoidIntegrationCells;
        Vector< moris::mtk::Vertex const * > mVerticesInCluster;
        moris::Matrix< moris::DDRMat >            mVertexParamCoords;

      public:
        Cell_Cluster_Proxy(){};

        //---------------------------------------------------------------------------------------

        virtual ~Cell_Cluster_Proxy(){};

        //---------------------------------------------------------------------------------------

        bool
        is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            return mTrivial;
        }

        //---------------------------------------------------------------------------------------

        Vector< moris::mtk::Cell const * > const &
        get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            return mPrimaryIntegrationCells;
        }

        //---------------------------------------------------------------------------------------

        Vector< moris::mtk::Cell const * > const &
        get_void_cells_in_cluster() const
        {
            return mVoidIntegrationCells;
        }

        //---------------------------------------------------------------------------------------

        moris::mtk::Cell const &
        get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            return *mInterpolationCell;
        }

        //---------------------------------------------------------------------------------------

        Vector< moris::mtk::Vertex const * >
        get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            return mVerticesInCluster;
        }

        //---------------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            return mVertexParamCoords;
        }

        //---------------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        get_vertex_local_coordinate_wrt_interp_cell(
                moris::mtk::Vertex const  *aVertex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            MORIS_ERROR( 0, "get_vertex_local_coordinate_wrt_interp_cell not implemented in proxy cell cluster" );
            return moris::Matrix< moris::DDRMat >( 0, 0 );
        }

        //---------------------------------------------------------------------------------------

        moris_index
        get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
        {
            return mVertexParamCoords.n_cols();
        }

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const
        {
            moris::real tVolume = 0.0;

            Vector< moris::mtk::Cell const * > const *tCells = nullptr;

            if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
            {
                tCells = &this->get_primary_cells_in_cluster();
            }
            else
            {
                tCells = &this->get_void_cells_in_cluster();
            }

            for ( auto iC = tCells->cbegin(); iC < tCells->cend(); iC++ )
            {
                tVolume += ( *iC )->compute_cell_measure();
            }

            return tVolume;
        }

        //---------------------------------------------------------------------------------------

        Matrix< DDRMat >
        compute_cluster_ig_cell_measures(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const
        {
            Vector< moris::mtk::Cell const * > const *tCells = nullptr;

            if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY )
            {
                tCells = &this->get_primary_cells_in_cluster();
            }
            else
            {
                tCells = &this->get_void_cells_in_cluster();
            }

            Matrix< DDRMat > tMeasureVec( tCells->size(), 1 );

            for ( uint iC = 0; iC < tCells->size(); ++iC )
            {
                tMeasureVec( iC ) = ( *tCells )( iC )->compute_cell_measure();
            }

            return tMeasureVec;
        }

        //---------------------------------------------------------------------------------------

        moris::real
        compute_cluster_cell_side_measure(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const
        {
            MORIS_ERROR( 0, "compute_cluster_cell_side_measure only valid on side clusters" );
            return 0;
        }

        //---------------------------------------------------------------------------------------

        Matrix< DDRMat >
        compute_cluster_ig_cell_side_measures(
                const mtk::Primary_Void    aPrimaryOrVoid = mtk::Primary_Void::PRIMARY,
                const mtk::Leader_Follower aIsLeader      = mtk::Leader_Follower::LEADER ) const
        {
            MORIS_ERROR( 0, "compute_cluster_cell_side_measure only valid on side clusters" );
            return { { 0 } };
        }
    };
}    // namespace moris
#endif /* PROJECTS_MTK_TEST_CL_MTK_CELL_CLUSTER_PROXY_HPP_ */
