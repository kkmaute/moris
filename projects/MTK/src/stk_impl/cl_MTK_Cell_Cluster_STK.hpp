/*
 * cl_MTK_Cell_Cluster_STK.hpp
 *
 *  Created on: Apr 29, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_CLUSTER_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_CLUSTER_STK_HPP_

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace mtk
{

class Cell_Cluster_STK: public Cell_Cluster
{
private:
    moris::mtk::Cell const *                mInterpolationCell;
    moris::Cell<moris::mtk::Cell const *>   mPrimaryIntegrationCells;
    moris::Cell<moris::mtk::Cell const *>   mVoidIntegrationCells;
    moris::Cell<moris::mtk::Vertex const *> mVerticesInCluster;
    moris::Matrix<moris::DDRMat>            mVertexParamCoords;
public:
    Cell_Cluster_STK():
        mInterpolationCell(nullptr),
        mPrimaryIntegrationCells(0,nullptr),
        mVoidIntegrationCells(0,nullptr),
        mVerticesInCluster(0,nullptr),
        mVertexParamCoords(0,0)
    {};

    void
    set_interpolation_cell(moris::mtk::Cell const * aInterpCell)
    {
        MORIS_ASSERT(mInterpolationCell == nullptr,"Interpolation Cell already set");
        mInterpolationCell = aInterpCell;
    }

    void
    add_primary_integration_cell(moris::mtk::Cell  const * aIntegrationCell)
    {
        mPrimaryIntegrationCells.push_back( aIntegrationCell );
    }

    void
    add_primary_integration_cell(moris::Cell<moris::mtk::Cell  const *> const & aIntegrationCell)
    {
        mPrimaryIntegrationCells.append( aIntegrationCell );
    }

    void
    add_void_integration_cell(moris::Cell<moris::mtk::Cell const *> const & aIntegrationCell)
    {
        mVoidIntegrationCells.append( aIntegrationCell );
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


    moris::Cell<moris::mtk::Cell const *> const &
    get_primary_cells_in_cluster() const
    {
        return mPrimaryIntegrationCells;
    }

    moris::Cell<moris::mtk::Cell const *> const &
    get_void_cells_in_cluster() const
    {
        return mVoidIntegrationCells;
    }

    moris::mtk::Cell const &
    get_interpolation_cell() const
    {
        return *mInterpolationCell;
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


#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_CLUSTER_STK_HPP_ */
