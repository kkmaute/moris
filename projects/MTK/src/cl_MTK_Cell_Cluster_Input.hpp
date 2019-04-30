/*
 * cl_MTK_Cell_Cluster_Input.hpp
 *
 *  Created on: Apr 29, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_INPUT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_INPUT_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Matrix.hpp"
#include <unordered_map>
namespace moris
{

namespace mtk
{
///////////////////////////////////
// STRUC FOR CELL CLUSTER INPUT  //
///////////////////////////////////

struct Cell_Cluster_Input
{
    Cell_Cluster_Input()
    {

    }

    void
    add_cluster_data(mtk::Cell*               aInterpCell,
                     moris::Matrix<IndexMat>* aPrimaryCellIds,
                     moris::Matrix<IndexMat>* aVoidCellIds,
                     moris::Matrix<IndexMat>* aVerticesInCluster,
                     moris::Matrix<DDRMat>*   aLocalCoordsRelativeToInterpCell)
    {
        MORIS_ASSERT(aVerticesInCluster->numel() == aLocalCoordsRelativeToInterpCell->n_rows(),"Number of vertices in the cluster must match the number of rows in local coord mat");

        moris::uint tIndex = mInterpolationCells.size();

        // add to map
        MORIS_ASSERT(mInterpCellIndexToIndex.find(aInterpCell->get_id())==mInterpCellIndexToIndex.end(),"Trying to add cluster on an interpolation cell which already has a cluster.");
        mInterpCellIndexToIndex[aInterpCell->get_id()] = tIndex;

        mInterpolationCells.push_back(aInterpCell);
        mPrimaryCellIds.push_back(aPrimaryCellIds);
        mVoidCellsIds.push_back(aVoidCellIds);
        mVerticesIdsInCluster.push_back(aVerticesInCluster);
        mLocalCoordsRelativeToInterpCell.push_back(aLocalCoordsRelativeToInterpCell);
    }

    moris_index
    get_cluster_index(moris_id aInterpCellId)
    {
        auto tIter = mInterpCellIndexToIndex.find(aInterpCellId);

        if(tIter == mInterpCellIndexToIndex.end())
        {
            return MORIS_INDEX_MAX;
        }

        else
        {
            return tIter->second;
        }
    }

    moris::uint
    get_num_cell_clusters()
    {
        return mInterpolationCells.size();
    }

    mtk::Cell*
    get_interp_cell(moris_index aClusterIndex)
    {
        MORIS_ASSERT(aClusterIndex< (moris_index)get_num_cell_clusters()," Cell cluster index out of bounds");
        return mInterpolationCells(aClusterIndex);
    }

    moris::Matrix<IndexMat>*
    get_primary_cell_ids(moris_index aClusterIndex)
    {
        MORIS_ASSERT(aClusterIndex< (moris_index)get_num_cell_clusters()," Cell cluster index out of bounds");
        return mPrimaryCellIds(aClusterIndex);
    }

    moris::Matrix<IndexMat>*
    get_void_cell_ids(moris_index aClusterIndex)
    {
        MORIS_ASSERT(aClusterIndex< (moris_index)get_num_cell_clusters()," Cell cluster index out of bounds");
        return mVoidCellsIds(aClusterIndex);
    }

    moris::Matrix<IndexMat>*
    get_vertex_in_cluster_ids(moris_index aClusterIndex)
    {
        MORIS_ASSERT(aClusterIndex<(moris_index) get_num_cell_clusters()," Cell cluster index out of bounds");
        return mVerticesIdsInCluster(aClusterIndex);
    }

    moris::Matrix<DDRMat>*
    get_vertex_local_coords_wrt_interpolation_cell(moris_index aClusterIndex)
    {
        MORIS_ASSERT(aClusterIndex< (moris_index)get_num_cell_clusters()," Cell cluster index out of bounds");
        return mLocalCoordsRelativeToInterpCell(aClusterIndex);
    }


private:
    moris::Cell<mtk::Cell*>               mInterpolationCells;
    moris::Cell<moris::Matrix<IndexMat>*> mPrimaryCellIds;
    moris::Cell<moris::Matrix<IndexMat>*> mVoidCellsIds;
    moris::Cell<moris::Matrix<IndexMat>*> mVerticesIdsInCluster;
    moris::Cell<moris::Matrix<DDRMat>*>   mLocalCoordsRelativeToInterpCell;

    std::unordered_map<moris_index, moris_index> mInterpCellIndexToIndex;


};

}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_INPUT_HPP_ */
