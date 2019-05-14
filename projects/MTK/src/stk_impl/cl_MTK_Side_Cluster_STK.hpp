/*
 * cl_MTK_Side_Cluster_STK.hpp
 *
 *  Created on: May 13, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_

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
    Side_Cluster_STK( moris::mtk::Cell const * aInterpCell):
        mTrivial(true),
        mInterpolationCell(aInterpCell),
        mIntegrationCells(0,nullptr),
        mIntegrationCellSideOrdinals(0,0),
        mVerticesInCluster(0,nullptr),
        mVertexParamCoords(0,0)
    {};

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
    {};

    void
    set_interpolation_cell(moris::mtk::Cell const * aInterpCell)
    {
        MORIS_ASSERT(mInterpolationCell == nullptr,"Interpolation Cell already set");
        mInterpolationCell = aInterpCell;
    }

    void
    mark_as_non_trivial()
    {
        mTrivial = false;
    }

    void
    add_integration_cell(moris::mtk::Cell  const * aIntegrationCell)
    {
        mIntegrationCells.push_back( aIntegrationCell );
    }

    void
    add_integration_cells(moris::Cell<moris::mtk::Cell  const *> const & aIntegrationCell)
    {
        mIntegrationCells.append( aIntegrationCell );
    }

    void
    add_integration_cell_side_ordinals(moris::Matrix<moris::IndexMat> const & aIntegrationCellSideOrds)
    {
        mIntegrationCellSideOrdinals = aIntegrationCellSideOrds.copy();
    }

    void
    add_vertex_to_cluster(moris::Cell<moris::mtk::Vertex const *> const & aVertex)
    {
        mVerticesInCluster.append(aVertex);
    }

    void
    add_vertex_local_coordinates_wrt_interp_cell(moris::Matrix<moris::DDRMat> const & aLocalCoords)
    {
        MORIS_ASSERT(aLocalCoords.n_rows() == mVerticesInCluster.size(),"Local coordinates need to match the number of vertices in the cluster");
        mVertexParamCoords = aLocalCoords.copy();
    }


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

    moris::Cell<moris::mtk::Vertex const *> const &
    get_vertices_in_cluster() const
    {
        return mVerticesInCluster;
    }

    moris::Matrix<moris::DDRMat> const &
    get_vertices_local_coordinates_wrt_interp_cell() const
    {
        return mVertexParamCoords;
    }
};
}
}



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_SIDE_CLUSTER_STK_HPP_ */
