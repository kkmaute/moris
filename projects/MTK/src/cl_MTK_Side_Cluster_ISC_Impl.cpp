/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Cluster_ISC_Impl.cpp
 *
 */

#include "cl_MTK_Side_Cluster_ISC_Impl.hpp"

namespace moris
{
    namespace mtk
    {
        //----------------------------------------------------------------
        Side_Cluster_ISC::Side_Cluster_ISC():
                    mTrivial(true),
                    mInterpolationCell(nullptr),
                    mIntegrationCells(0,nullptr),
                    mIntegrationCellSideOrdinals(0,0),
                    mVerticesInCluster(0,nullptr)
                    //mVertexParamCoords(0,0)
                    {}
        //----------------------------------------------------------------

        Side_Cluster_ISC::Side_Cluster_ISC( moris::mtk::Cell const *                        aInterpCell,
                moris::mtk::Cell const *                        aIntegrationCell,
                Vector<moris::mtk::Vertex const *> const & aVerticesInCluster,
                moris_index aSideOrdinal):
                    mTrivial(true),
                    mInterpolationCell(aInterpCell),
                    mIntegrationCells(0,nullptr),
                    mIntegrationCellSideOrdinals({{aSideOrdinal}}),
                    mVerticesInCluster(aVerticesInCluster)
                    //mVertexParamCoords(0,0)
        {
            mIntegrationCells.push_back(aIntegrationCell);

            // add to map
            for(moris::uint i = 0; i <aVerticesInCluster.size(); i++)
            {
                this->add_vertex_to_map(aVerticesInCluster(i)->get_id(),i);
            }

        }

        //----------------------------------------------------------------

        Side_Cluster_ISC::Side_Cluster_ISC(
                bool                                            aTrivial,
                moris::mtk::Cell const *                       aInterpolationCell,
                Vector<moris::mtk::Cell const *>   const & aIntegrationCells,
                moris::Matrix<moris::IndexMat>          const & aIntegrationCellSideOrdinals,
                Vector<moris::mtk::Vertex const *> const & aVerticesInCluster,
                moris::Matrix<moris::DDRMat>            const & aVertexParamCoords):
                    mTrivial(aTrivial),
                    mInterpolationCell(aInterpolationCell),
                    mIntegrationCells(aIntegrationCells),
                    mIntegrationCellSideOrdinals(aIntegrationCellSideOrdinals),
                    mVerticesInCluster(aVerticesInCluster),
                    mVertexParamCoords(aVertexParamCoords)
        {
            MORIS_ERROR(aVerticesInCluster.size() == aVertexParamCoords.n_rows(),"Dimension mismatch between parametric coordinates provided and vertices provided (one row in parametric coordinates for each vertex in cluster).");

            // add vertices to map
            for(moris::uint i = 0; i <aVerticesInCluster.size(); i++)
            {
                this->add_vertex_to_map(aVerticesInCluster(i)->get_id(),i);
            }
        }

        //----------------------------------------------------------------

        bool
        Side_Cluster_ISC::is_trivial( const mtk::Leader_Follower aIsLeader ) const
        {
            return mTrivial;
        }

        //----------------------------------------------------------------

        moris::mtk::Cell const &
        Side_Cluster_ISC::get_interpolation_cell(const mtk::Leader_Follower aIsLeader ) const
        {
            return *mInterpolationCell;
        }

        //----------------------------------------------------------------

        Vector<moris::mtk::Cell const *> const &
        Side_Cluster_ISC::get_cells_in_side_cluster() const
        {
            return mIntegrationCells;
        }

        //----------------------------------------------------------------

        moris::Matrix<moris::IndexMat>
        Side_Cluster_ISC::get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const
        {
            return mIntegrationCellSideOrdinals;
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_ISC::get_cell_side_ordinal(moris::moris_index aCellIndexInCluster,
                const mtk::Leader_Follower aIsLeader  ) const
        {
            MORIS_ASSERT(aCellIndexInCluster<(moris_index)mIntegrationCellSideOrdinals.numel(),"Cell index in cluster out of bounds");

            return mIntegrationCellSideOrdinals(aCellIndexInCluster);
        }

        //----------------------------------------------------------------

        Vector<moris::mtk::Vertex const *>
        Side_Cluster_ISC::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
        {
            return mVerticesInCluster;
        }

        //----------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Side_Cluster_ISC::get_vertices_local_coordinates_wrt_interp_cell(const mtk::Leader_Follower aIsLeader ) const
        {
            return mVertexParamCoords;
            //return moris::Matrix<moris::DDRMat>();;
        }

        //----------------------------------------------------------------

        moris::moris_index
        Side_Cluster_ISC::get_vertex_cluster_index( moris::mtk::Vertex const * aVertex,
                const mtk::Leader_Follower    aIsLeader ) const
        {

            return this->get_vertex_cluster_local_index(aVertex->get_id());
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_ISC::get_vertex_ordinal_on_facet( moris_index aCellIndexInCluster, moris::mtk::Vertex const * aVertex ) const
        {
            moris_index tSideOrd = mIntegrationCellSideOrdinals(aCellIndexInCluster);

            Vector<moris::mtk::Vertex const *> tVerticesOnSide = mIntegrationCells(aCellIndexInCluster)->get_vertices_on_side_ordinal(tSideOrd);

            // iterate through vertices and see if the ids match
            for(moris::moris_index i = 0; i < (moris_index)tVerticesOnSide.size(); i++)
            {

                if(tVerticesOnSide(i)->get_id() == aVertex->get_id())
                {
                    return i;
                }
            }
            return MORIS_INDEX_MAX;
        }

        //----------------------------------------------------------------

        moris::Matrix<moris::DDRMat>
        Side_Cluster_ISC::get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                const mtk::Leader_Follower aIsLeader ) const
        {
           moris_index tLocalVertIndex = this->get_vertex_cluster_local_index(aVertex->get_id());

            MORIS_ASSERT( tLocalVertIndex < (moris_index)mVertexParamCoords.n_rows(),"Vertex local side cluster index out of bounds. This could be cause by not adding parametric coordinates");

            return mVertexParamCoords.get_row(tLocalVertIndex);
            //return moris::Matrix<moris::DDRMat>();
        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_ISC::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader ) const
        {
            return this->get_vertices_local_coordinates_wrt_interp_cell(aIsLeader).n_cols();

        }

        //----------------------------------------------------------------

        moris_index
        Side_Cluster_ISC::get_vertex_cluster_local_index(moris_id aVertexId) const
        {

            auto tIter = mVertexIdToLocalIndex.find(aVertexId);

            MORIS_ERROR(tIter != mVertexIdToLocalIndex.end(),"Vertex not found in side cluster");

            return tIter->second;
        }

        //----------------------------------------------------------------

        void
        Side_Cluster_ISC::add_vertex_to_map(moris_id aVertexId,
                moris_index aVertexLocalIndex)
        {
            MORIS_ERROR(mVertexIdToLocalIndex.find(aVertexId) == mVertexIdToLocalIndex.end(),"Trying to add vertex already found in side cluster");
            mVertexIdToLocalIndex[aVertexId] = aVertexLocalIndex;
        }

    } /*namespace mtk */
} /*namespace moris */

