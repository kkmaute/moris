/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Cell_Cluster_Visualization.cpp
 *
 */

#include "cl_VIS_Cell_Cluster_Visualization.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"

namespace moris::vis
{
    //----------------------------------------------------------------
    bool
    Cell_Cluster_Visualization::is_trivial( const mtk::Leader_Follower aIsLeader ) const
    {
        return mTrivial;
    }

    // ##############################################
    //  Add and setup of cluster
    // ##############################################
    //----------------------------------------------------------------
    void
    Cell_Cluster_Visualization::mark_as_nontrivial()
    {
        mTrivial = false;
    }
    //----------------------------------------------------------------

    void
    Cell_Cluster_Visualization::set_interpolation_cell( moris::mtk::Cell const *aInterpCell )
    {
        MORIS_ASSERT( mInterpolationCell == nullptr, "Interpolation Cell already set" );
        mInterpolationCell = aInterpCell;
    }

    //----------------------------------------------------------------
    void
    Cell_Cluster_Visualization::add_primary_integration_cell( moris::mtk::Cell const *aIntegrationCell )
    {
        mPrimaryIntegrationCells.push_back( aIntegrationCell );
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster_Visualization::add_primary_integration_cell( Vector< moris::mtk::Cell const * > const &aIntegrationCell )
    {
        mPrimaryIntegrationCells.append( aIntegrationCell );
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster_Visualization::add_void_integration_cell( Vector< moris::mtk::Cell const * > const &aIntegrationCell )
    {

        mVoidIntegrationCells.append( aIntegrationCell );
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster_Visualization::add_vertex_to_cluster( Vector< moris::mtk::Vertex const * > const &aVertex )
    {
        // add vertices to map
        moris_index tIndex = mVerticesInCluster.size();

        // add vertices to map
        for ( moris::uint i = 0; i < aVertex.size(); i++ )
        {
            this->add_vertex_to_map( aVertex( i )->get_id(), tIndex );
            tIndex++;
        }

        mVerticesInCluster.append( aVertex );
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster_Visualization::add_vertex_local_coordinates_wrt_interp_cell( moris::Matrix< moris::DDRMat > const &aLocalCoords )
    {
        MORIS_ASSERT( aLocalCoords.n_rows() == mVerticesInCluster.size(),
                "Local coordinates need to match the number of vertices in the cluster" );

        mVertexParamCoords = aLocalCoords.copy();
    }

    //----------------------------------------------------------------

    Vector< moris::mtk::Cell const * > const &
    Cell_Cluster_Visualization::get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return mPrimaryIntegrationCells;
    }

    //----------------------------------------------------------------

    Vector< moris::mtk::Cell const * > const &
    Cell_Cluster_Visualization::get_void_cells_in_cluster() const
    {
        return mVoidIntegrationCells;
    }

    //----------------------------------------------------------------

    moris::mtk::Cell const &
    Cell_Cluster_Visualization::get_interpolation_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        return *mInterpolationCell;
    }

    //----------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Cell_Cluster_Visualization::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return mVerticesInCluster;
    }

    //----------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell_Cluster_Visualization::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        return mVertexParamCoords;
    }

    //----------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell_Cluster_Visualization::get_vertex_local_coordinate_wrt_interp_cell(
            moris::mtk::Vertex const  *aVertex,
            const mtk::Leader_Follower aIsLeader ) const
    {
        moris_index tLocalVertIndex = this->get_vertex_cluster_local_index( aVertex->get_id() );

        MORIS_ASSERT( tLocalVertIndex < (moris_index)mVertexParamCoords.n_rows(),
                "Vertex local side cluster index out of bounds. This could be cause by not adding parametric coordinates" );

        return mVertexParamCoords.get_row( tLocalVertIndex );
    }

    //----------------------------------------------------------------

    moris_index
    Cell_Cluster_Visualization::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader ) const
    {
        return mVertexParamCoords.n_cols();
    }

    //----------------------------------------------------------------

    moris_index
    Cell_Cluster_Visualization::get_vertex_cluster_local_index( moris_id aVertexId ) const
    {
        auto tIter = mVertexIdToLocalIndex.find( aVertexId );

        MORIS_ERROR( tIter != mVertexIdToLocalIndex.end(), "Vertex not found in side cluster" );

        return tIter->second;
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster_Visualization::add_vertex_to_map( moris_id aVertexId,
            moris_index                                     aVertexLocalIndex )
    {
        MORIS_ERROR( mVertexIdToLocalIndex.find( aVertexId ) == mVertexIdToLocalIndex.end(),
                "Trying to add vertex already found in side cluster" );

        mVertexIdToLocalIndex[ aVertexId ] = aVertexLocalIndex;
    }

    //----------------------------------------------------------------
}    // namespace moris::vis
