/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Side_Cluster_Visualization.cpp
 *
 */

#include "cl_VIS_Side_Cluster_Visualization.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Cluster.hpp"

namespace moris
{
    namespace vis
    {
        //----------------------------------------------------------------
        bool
        Side_Cluster_Visualization::is_trivial( const mtk::Leader_Follower aIsLeader ) const
        {
            return mTrivial;
        }

        // ##############################################
        //  Add and setup of cluster
        // ##############################################

        //----------------------------------------------------------------
        void
        Side_Cluster_Visualization::mark_as_nontrivial()
        {
            mTrivial = false;
        }
        //----------------------------------------------------------------

        void
        Side_Cluster_Visualization::set_interpolation_cell( moris::mtk::Cell const *aInterpCell )
        {
            MORIS_ASSERT( mInterpolationCell == nullptr,
                    "Side_Cluster_Visualization::set_interpolation_cell() - Interpolation Cell already set" );
            mInterpolationCell = aInterpCell;
        }

        //----------------------------------------------------------------
        void
        Side_Cluster_Visualization::add_primary_integration_cell( moris::mtk::Cell const *aIntegrationCell )
        {
            mPrimaryIntegrationCells.push_back( aIntegrationCell );
        }

        //----------------------------------------------------------------

        void
        Side_Cluster_Visualization::add_primary_integration_cell( moris::Cell< moris::mtk::Cell const * > const &aIntegrationCell )
        {
            mPrimaryIntegrationCells.append( aIntegrationCell );
        }

        //----------------------------------------------------------------

        void
        Side_Cluster_Visualization::add_integration_cell_side_ordinals( moris::Matrix< moris::IndexMat > aIntegrationCellSideOrdinals )
        {

            mIntegrationCellSideOrdinals = aIntegrationCellSideOrdinals.copy();
        }

        //----------------------------------------------------------------

        void
        Side_Cluster_Visualization::add_integration_cell_side_ordinals( moris_index aSideOrdinal )
        {

            mIntegrationCellSideOrdinals = { { aSideOrdinal } };
        }

        //----------------------------------------------------------------

        void
        Side_Cluster_Visualization::add_vertex_to_cluster( moris::Cell< moris::mtk::Vertex const * > const &aVertex )
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
        Side_Cluster_Visualization::add_vertex_local_coordinates_wrt_interp_cell( moris::Matrix< moris::DDRMat > const &aLocalCoords )
        {
            MORIS_ASSERT( aLocalCoords.n_rows() == mVerticesInCluster.size(), "Local coordinates need to match the number of vertices in the cluster" );
            mVertexParamCoords = aLocalCoords.copy();
        }

        //----------------------------------------------------------------

        // ##############################################
        //  Required Access Functions
        // ##############################################

        moris::mtk::Cell const &
        Side_Cluster_Visualization::get_interpolation_cell( const mtk::Leader_Follower aIsLeader ) const
        {
            return *mInterpolationCell;
        }

        //----------------------------------------------------------------

        moris::Cell< moris::mtk::Cell const * > const &
        Side_Cluster_Visualization::get_cells_in_side_cluster() const
        {
            return mPrimaryIntegrationCells;
        }

        //----------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        Side_Cluster_Visualization::get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const
        {
            return mIntegrationCellSideOrdinals;
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_Visualization::get_cell_side_ordinal(
                moris::moris_index         aCellIndexInCluster,
                const mtk::Leader_Follower aIsLeader ) const
        {
            return mIntegrationCellSideOrdinals( aCellIndexInCluster );
        }

        //----------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex const * >
        Side_Cluster_Visualization::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
        {
            return mVerticesInCluster;
        }

        //----------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        Side_Cluster_Visualization::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const
        {
            return mVertexParamCoords;
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_Visualization::get_vertex_cluster_index(
                const mtk::Vertex         *aVertex,
                const mtk::Leader_Follower aIsLeader ) const
        {
            moris_id tVertexId = aVertex->get_id();
            return this->get_vertex_cluster_local_index( tVertexId );
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_Visualization::get_vertex_ordinal_on_facet(
                moris_index               aCellIndexInCluster,
                moris::mtk::Vertex const *aVertex ) const
        {
            // get the side ordinal of the facet
            moris_index tSideOrd = mIntegrationCellSideOrdinals( aCellIndexInCluster );

            // get the vertices 
            moris::Cell< mtk::Vertex const * > tVerticesOnSide = 
                    mPrimaryIntegrationCells( aCellIndexInCluster )->get_vertices_on_side_ordinal( tSideOrd );
            
            // iterate through vertices and see if the ids match
            for ( moris_index iVertex = 0; iVertex < (moris_index)tVerticesOnSide.size(); iVertex++ )
            {
                moris_index tVertexIndex = tVerticesOnSide( iVertex )->get_index();
                if ( tVertexIndex == aVertex->get_index() )
                {
                    return iVertex;
                }
            }

            // if the code makes it here, the vertex has not been found
            MORIS_ERROR( false, "Side_Cluster_Visualization::get_vertex_ordinal_on_facet() - Vertex not found on facet." );
            return MORIS_INDEX_MAX;
        }

        //----------------------------------------------------------------

        moris::Matrix< moris::DDRMat >
        Side_Cluster_Visualization::get_vertex_local_coordinate_wrt_interp_cell(
                moris::mtk::Vertex const  *aVertex,
                const mtk::Leader_Follower aIsLeader ) const
        {
            moris_index tLocalVertIndex = this->get_vertex_cluster_local_index( aVertex->get_id() );

            MORIS_ASSERT( tLocalVertIndex < (moris_index)mVertexParamCoords.n_rows(),
                    "Side_Cluster_Visualization::get_vertex_local_coordinate_wrt_interp_cell() - "
                    "Vertex local side cluster index out of bounds. This could be cause by not adding parametric coordinates" );

            return mVertexParamCoords.get_row( tLocalVertIndex );
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_Visualization::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader ) const
        {
            return mVertexParamCoords.n_cols();
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_Visualization::get_vertex_cluster_local_index( moris_id aVertexId ) const
        {
            MORIS_ASSERT( mVertexIdToLocalIndex.size() > 0,
                    "Side_Cluster_Visualization::get_vertex_cluster_local_index() - "
                    "mVertexIdToLocalIndex has not been initialized." );

            auto tIter = mVertexIdToLocalIndex.find( aVertexId );

            MORIS_ERROR( tIter != mVertexIdToLocalIndex.end(),
                    "Side_Cluster_Visualization::get_vertex_cluster_local_index() - "
                    "Vertex not found in side cluster." );

            return tIter->second;
        }

        //----------------------------------------------------------------

        void
        Side_Cluster_Visualization::add_vertex_to_map(
                moris_id    aVertexId,
                moris_index aVertexLocalIndex )
        {
            MORIS_ERROR( mVertexIdToLocalIndex.find( aVertexId ) == mVertexIdToLocalIndex.end(), "Trying to add vertex already found in side cluster" );
            mVertexIdToLocalIndex[ aVertexId ] = aVertexLocalIndex;
        }

        //----------------------------------------------------------------

    }    // namespace vis
}    // namespace moris
