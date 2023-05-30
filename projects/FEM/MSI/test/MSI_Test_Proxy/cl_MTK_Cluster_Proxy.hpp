/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cluster_Proxy.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_CLUSTER_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_CLUSTER_STK_HPP_

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_Matrix.hpp"
#include <unordered_map>

namespace moris
{
namespace mtk
{

class Cluster_Proxy: public Cluster
{
private:
    bool                                    mTrivial;
    moris::mtk::Cell const *                mInterpolationCell;
    moris::Cell<moris::mtk::Cell const *>   mPrimaryIntegrationCells;
    moris::Cell<moris::mtk::Cell const *>   mVoidIntegrationCells;
    moris::Cell<moris::mtk::Vertex const *> mVerticesInCluster;
    moris::Matrix<moris::DDRMat>            mVertexParamCoords;

    // map from vertex id to local index
    std::unordered_map<moris_index,moris_index> mVertexIdToLocalIndex;   // FIXME should be ordered map. about 1000 times faster

public:
    Cluster_Proxy( moris::Cell<moris::mtk::Cell const *> aPrimaryCells,
                   moris::Cell<moris::mtk::Cell const *> aVoidCells,
                   moris::Matrix<moris::DDRMat>          aLocalCoords   ) : mTrivial(true),
                      mInterpolationCell(nullptr),
                      mPrimaryIntegrationCells(aPrimaryCells),
                      mVoidIntegrationCells(aVoidCells),
                      mVerticesInCluster(0,nullptr),
                      mVertexParamCoords(aLocalCoords)
    {};

    //----------------------------------------------------------------

    ~Cluster_Proxy(){};

    //----------------------------------------------------------------
    bool
    is_trivial( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
    { return false; }

    //##############################################
    // Add and setup of cluster
    //##############################################
    void
    mark_as_nontrivial(){};

    //----------------------------------------------------------------

    void
    set_interpolation_cell(moris::mtk::Cell const * aInterpCell){};

    //----------------------------------------------------------------

    void
    add_primary_integration_cell(moris::mtk::Cell  const * aIntegrationCell){};

    //----------------------------------------------------------------

    void
    add_vertex_to_cluster(moris::Cell<moris::mtk::Vertex const *> const & aVertex){};

    //----------------------------------------------------------------

    void
    add_vertex_local_coordinates_wrt_interp_cell(moris::Matrix<moris::DDRMat> const & aLocalCoords){};

    //----------------------------------------------------------------

    //##############################################
    // Required Access Functions
    //##############################################

    moris::Cell<moris::mtk::Cell const *> const &
    get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const
	{
    	return mPrimaryIntegrationCells;
	};

    //----------------------------------------------------------------

    moris::Cell<moris::mtk::Cell const *> const &
    get_void_cells_in_cluster() const
	{
    	return mVoidIntegrationCells;
	};

    //----------------------------------------------------------------

    moris::mtk::Cell const &
    get_interpolation_cell( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
    {     	MORIS_ERROR( false, "not implemented for proxy");
    return *mInterpolationCell; };

    //----------------------------------------------------------------

    moris::Cell<moris::mtk::Vertex const *> const &
    get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
	{     	MORIS_ERROR( false, "not implemented for proxy");
	return mVerticesInCluster;};

    //----------------------------------------------------------------

    moris::Matrix<moris::DDRMat> const &
    get_vertices_local_coordinates_wrt_interp_cell(const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const
	{
        return mVertexParamCoords;
	};

    //----------------------------------------------------------------

    moris::Matrix<moris::DDRMat>
    get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
            const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER) const
    {
    	MORIS_ERROR( false, "not implemented for proxy");
    	return moris::Matrix<moris::DDRMat>( 0,0 );};

    //----------------------------------------------------------------

    moris_index
    get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER ) const
    {
        MORIS_ERROR( false, "not implemented for proxy");
        return 0; };

};
}
}

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_CLUSTER_STK_HPP_ */

