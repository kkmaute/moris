/*
 * cl_MTK_Side_Cluster_STK.hpp
 *
 *  Created on: May 13, 2019
 *      Author: doble
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

class Side_Cluster_STK: public Side_Cluster
{
private:

    bool                                    mTrivial;
    moris::mtk::Cell const *                mInterpolationCell;
    moris::Cell<moris::mtk::Cell const *>   mIntegrationCells;
    moris::Matrix<moris::IndexMat>          mIntegrationCellSideOrdinals;
    moris::Cell<moris::mtk::Vertex const *> mVerticesInCluster;
    moris::Matrix<moris::DDRMat>            mVertexParamCoords;

    // map from vertex id to local index
    std::unordered_map<moris_index,moris_index> mVertexIdToLocalIndex;


public:
    Side_Cluster_STK():
        mTrivial(true),
        mInterpolationCell(nullptr),
        mIntegrationCells(0,nullptr),
        mIntegrationCellSideOrdinals(0,0),
        mVerticesInCluster(0,nullptr),
        mVertexParamCoords(0,0)
    {};


    // trivial constructor
    Side_Cluster_STK( moris::mtk::Cell const * aInterpCell,
                      moris::mtk::Cell const * aIntegrationCell,
                      moris_index aSideOrdinal):
        mTrivial(true),
        mInterpolationCell(aInterpCell),
        mIntegrationCells(0,nullptr),
        mIntegrationCellSideOrdinals({{aSideOrdinal}}),
        mVerticesInCluster(0,nullptr),
        mVertexParamCoords(0,0)
    {
        mIntegrationCells.push_back(aIntegrationCell);
    };

    Side_Cluster_STK(bool aTrivial,
                     moris::mtk::Cell const *                        aInterpolationCell,
                     moris::Cell<moris::mtk::Cell const *>   const & aIntegrationCells,
                     moris::Matrix<moris::IndexMat>          const & aIntegrationCellSideOrdinals,
                     moris::Cell<moris::mtk::Vertex const *> const & aVerticesInCluster,
                     moris::Matrix<moris::DDRMat> const & aVertexParamCoords):
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
    };


    bool
    is_trivial() const
    {
        return mTrivial;
    }

    moris::mtk::Cell const &
    get_interpolation_cell() const
    {
        return *mInterpolationCell;
    }

    moris::Cell<moris::mtk::Cell const *> const &
    get_cells_in_side_cluster() const
    {
        return mIntegrationCells;
    }

    moris::Matrix<moris::IndexMat>
    get_cell_side_ordinals() const
    {
        return mIntegrationCellSideOrdinals;
    }

    moris_index
    get_cell_side_ordinal(moris::moris_index aCellIndexInCluster) const
    {
        MORIS_ASSERT(aCellIndexInCluster<(moris_index)mIntegrationCellSideOrdinals.numel(),"Cell index in cluster out of bounds");

        return mIntegrationCellSideOrdinals(aCellIndexInCluster);
    }


    moris::Cell<moris::mtk::Vertex const *> const &
    get_vertices_in_cluster() const
    {
        MORIS_ERROR(!mTrivial,"Accessing vertices in a trivial cluster is not allowed");

        return mVerticesInCluster;
    }

    moris::Matrix<moris::DDRMat> const &
    get_vertices_local_coordinates_wrt_interp_cell() const
    {
        MORIS_ERROR(!mTrivial,"Accessing local coordinates on a trivial side cluster is not allowed");

        return mVertexParamCoords;
    }

    moris::Matrix<moris::DDRMat>
    get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex ) const
    {
        MORIS_ERROR(!mTrivial,"Accessing local coordinates on a trivial side cluster is not allowed");

        moris_index tLocalVertIndex = this->get_vertex_cluster_local_index(aVertex->get_id());

        MORIS_ASSERT( tLocalVertIndex < (moris_index)mVertexParamCoords.n_rows(),"Vertex local side cluster index out of bounds. This could be cause by not adding parametric coordinates");

        return mVertexParamCoords.get_row(tLocalVertIndex);
    }

    //##############################################
    // Size Access
    // (Pure Virtual)
    //##############################################
    /*!
     * Size of the xsi vector in this side cluster
     */
    moris_index
    get_dim_of_param_coord() const
    {
        MORIS_ERROR(!mTrivial,"Accessing size of local coordinates on a trivial side cluster is not allowed");
        return mVertexParamCoords.n_cols();
    }

    moris_index
    get_vertex_cluster_local_index(moris_id aVertexId) const
    {


        auto tIter = mVertexIdToLocalIndex.find(aVertexId);

        MORIS_ERROR(tIter != mVertexIdToLocalIndex.end(),"Vertex not found in side cluster");

        return tIter->second;
    }

    void
    add_vertex_to_map(moris_id aVertexId,
                      moris_index aVertexLocalIndex)
    {
        MORIS_ERROR(mVertexIdToLocalIndex.find(aVertexId) == mVertexIdToLocalIndex.end(),"Trying to add vertex already found in side cluster");
        mVertexIdToLocalIndex[aVertexId] = aVertexLocalIndex;
    }

};
}
}



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_ */
