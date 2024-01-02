/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_XTK_Cell_Cluster.cpp  
 * 
 */

#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_MTK_Cluster_Group.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "fn_stringify_matrix.hpp"

// namespace moris
// {
namespace xtk
{
    //----------------------------------------------------------------
    
    Cell_Cluster::Cell_Cluster()
            : mTrivial( true )
            , mVoid( false )
            , mInvalid( false )
            , mInterpolationCell( nullptr )
            , mChildMesh( nullptr )
            , mPrimaryIntegrationCells( 0, nullptr )
            , mVoidIntegrationCells( 0, nullptr )
            , mVerticesInCluster( 0, nullptr )
            , mClusterGroups( 0 )
    {
    }

    Cell_Cluster::Cell_Cluster( bool aOnlyForVisualization )
            : mTrivial( true )
            , mVoid( false )
            , mInvalid( false )
            , mInterpolationCell( nullptr )
            , mChildMesh( nullptr )
            , mPrimaryIntegrationCells( 0, nullptr )
            , mVoidIntegrationCells( 0, nullptr )
            , mVerticesInCluster( 0, nullptr )
            , mClusterGroups( 0 )
    {
        mOnlyForVis = aOnlyForVisualization;
    }

    //----------------------------------------------------------------

    Cell_Cluster::~Cell_Cluster(){}

    //----------------------------------------------------------------

    bool
    Cell_Cluster::is_trivial( const mtk::Leader_Follower aIsLeader ) const
    {
        return mTrivial;
    }

    //----------------------------------------------------------------

    bool
    Cell_Cluster::is_full() const
    {
        // Cluster is full if it is trivial and not void
        return ( mTrivial && !mVoid );
    }

    //----------------------------------------------------------------

    bool
    Cell_Cluster::is_void() const
    {
        return mVoid;
    }

    //----------------------------------------------------------------

    bool
    Cell_Cluster::is_invalid() const
    {
        return mVoid;
    }

    //----------------------------------------------------------------

    Vector<moris::mtk::Cell const *> const &
    Cell_Cluster::get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return mPrimaryIntegrationCells;
    }

    //----------------------------------------------------------------

    Vector<moris::mtk::Cell const *> const &
    Cell_Cluster::get_void_cells_in_cluster() const
    {
        return mVoidIntegrationCells;
    }

    //----------------------------------------------------------------

    moris::mtk::Cell const &
    Cell_Cluster::get_interpolation_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        return *mInterpolationCell;
    }

    //----------------------------------------------------------------

    Vector<moris::mtk::Vertex const *>
    Cell_Cluster::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
            return mVerticesInCluster;
    }

    //----------------------------------------------------------------

    moris::Matrix<moris::DDRMat>
    Cell_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader )  const
    {
        if(!mTrivial)
        {
            return mLocalCoords;
        }
        else
        {
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const * tCellInfo = mInterpolationCell->get_cell_info();

            // local coordinate matrix
            Matrix<DDRMat> tXi;

            // get the local coordinates on the side ordinal
            tCellInfo->get_loc_coords_of_cell(tXi);

            return tXi;
        }
    }

    //----------------------------------------------------------------

    moris::Matrix<moris::DDRMat>
    Cell_Cluster::get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
            const mtk::Leader_Follower aIsLeader ) const
    {
        MORIS_ERROR(!mTrivial,"Accessing local coordinates on a trivial cell cluster is not allowed");
        return *mVertexGroup->get_vertex_local_coords(aVertex->get_index());
    }

    //----------------------------------------------------------------

    moris_index
    Cell_Cluster::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader) const
    {
        return this->get_vertices_local_coordinates_wrt_interp_cell(aIsLeader).n_cols();
    }

    moris::Matrix<moris::DDRMat>                  
    Cell_Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aPrimaryCellClusterIndex) const
    {
        if(mTrivial)
        {
            return mLocalCoords;
        }
        else
        {
             // MORIS_ERROR(!this->is_trivial(),"get_primary_cell_local_coords_on_side_wrt_interp_cell on trivial cluster is not allowed");
             MORIS_ASSERT(aPrimaryCellClusterIndex < (moris_index)this->get_num_primary_cells(),"Integration Cell Cluster index out of bounds");

            // get the integration cell of interest
            moris::mtk::Cell const * tIntegrationCell = this->get_primary_cells_in_cluster()(aPrimaryCellClusterIndex);

             // get the vertex pointers on the side - for the bulk this is all vertices on the integration cell
            Vector<moris::mtk::Vertex *> tVerticesOnCell = tIntegrationCell->get_vertex_pointers();

            // allocate output (n_node x dim_xsi)
            moris::Matrix<moris::DDRMat> tVertexParamCoords( tVerticesOnCell.size(), this->get_dim_of_param_coord());

            // iterate through vertices and collect local coordinates
            for(moris::uint i = 0; i < tVerticesOnCell.size(); i++)
            {
                 tVertexParamCoords.get_row(i) = this->get_vertex_local_coordinate_wrt_interp_cell(tVerticesOnCell(i)).get_row(0);
            }


            return tVertexParamCoords;
        }
    }

    //----------------------------------------------------------------
    Interpolation_Cell_Unzipped const *
    Cell_Cluster::get_xtk_interpolation_cell() const
    {
        return mInterpolationCell;
    }

    //----------------------------------------------------------------

    Matrix< IndexMat > Cell_Cluster::get_hanging_nodes(  ) const
    {
        MORIS_ERROR(0,"FIXME");
        return mChildMesh->get_hanging_nodes(  );
    }

    //----------------------------------------------------------------
 
    size_t
    Cell_Cluster::capacity()
    {   
        size_t tTotalSize = 0;
        tTotalSize += sizeof(mTrivial);
        tTotalSize += sizeof(mInterpolationCell);
        tTotalSize += sizeof(mChildMesh);
        tTotalSize += mPrimaryIntegrationCells.capacity();
        tTotalSize += mVoidIntegrationCells.capacity();
        tTotalSize += mVerticesInCluster.capacity();
        return tTotalSize;
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::set_primary_integration_cell_group( std::shared_ptr<IG_Cell_Group> aPrimaryIgCells )
    {
        mPrimaryIgCellGroup = { aPrimaryIgCells };

        mPrimaryIntegrationCells.resize(aPrimaryIgCells->mIgCellGroup.size());

        for(moris::uint i = 0; i < aPrimaryIgCells->mIgCellGroup.size(); i++)
        {
            mPrimaryIntegrationCells(i) = aPrimaryIgCells->mIgCellGroup(i);
        }
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::set_primary_integration_cell_groups( Vector< std::shared_ptr< IG_Cell_Group > > aPrimaryIgCells )
    {
        // store IG cell groups with cell cluster
        mPrimaryIgCellGroup = aPrimaryIgCells;

        // count total number of IG cells in all groups passed into function
        moris::uint tCount = 0;
        for( moris::uint iCellGroup = 0; iCellGroup < aPrimaryIgCells.size(); iCellGroup++ )
        {
            tCount = tCount + aPrimaryIgCells( iCellGroup )->mIgCellGroup.size();
        }

        // initialize list of IG cells
        mPrimaryIntegrationCells.resize( tCount );

        // reset counter to track position in list
        tCount = 0;

        // store IG cells in list
        for( moris::uint iCellGroup = 0; iCellGroup < aPrimaryIgCells.size(); iCellGroup++ )
        {
            for( moris::uint jCellInGroup = 0; jCellInGroup < aPrimaryIgCells( iCellGroup )->mIgCellGroup.size(); jCellInGroup++ )
            {
                mVoidIntegrationCells( tCount ) = aPrimaryIgCells( iCellGroup )->mIgCellGroup( jCellInGroup );
                tCount++;
            }
        }
    }

    //----------------------------------------------------------------

    void
    Cell_Cluster::set_void_integration_cell_groups(Vector<std::shared_ptr<IG_Cell_Group>> & aVoidIgCells)
    {
        mVoidIgCellGroup = aVoidIgCells;

        moris::uint tCount = 0;
        for(moris::uint i = 0; i < aVoidIgCells.size(); i++)
        {
            tCount = tCount + aVoidIgCells(i)->mIgCellGroup.size();
        }


        mVoidIntegrationCells.resize(tCount);

        tCount = 0;

        for(moris::uint i = 0; i < aVoidIgCells.size(); i++)
        {
            for(moris::uint j = 0; j < aVoidIgCells(i)->mIgCellGroup.size(); j++)
            {
                mVoidIntegrationCells(tCount++) = aVoidIgCells(i)->mIgCellGroup(j);
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    Cell_Cluster::set_ig_vertex_group(std::shared_ptr<IG_Vertex_Group> aVertexGroup)
    {
        mVertexGroup = aVertexGroup;

        mVerticesInCluster.resize(mVertexGroup->size());
        mLocalCoords.resize(mVertexGroup->size(),mVertexGroup->get_vertex_local_coords(aVertexGroup->get_vertex(0)->get_index())->n_cols());

        for(moris::uint i = 0; i < mVertexGroup->size(); i++)
        {
            mVerticesInCluster(i) = mVertexGroup->get_vertex(i);
            mLocalCoords.set_row(i, *mVertexGroup->get_vertex_local_coords(mVerticesInCluster(i)->get_index()));
        }
    }


    //----------------------------------------------------------------

    bool
    Cell_Cluster::has_cluster_group( const moris_index aDiscretizationMeshIndex ) const
    {
        // if the list of B-spline meshes for which the cluster groups have been set isn't large enough ...
        if( mClusterGroups.size() < (uint) aDiscretizationMeshIndex + 1 )
        {
            // ... the cluster group has not been set yet
            return false;
        }

        // if the entry exists but is empty
        else if( mClusterGroups( aDiscretizationMeshIndex ).expired() )
        {
            return false;
        }

        // if the entry both exists and is not empty, the cluster has a group
        else
        {
            return true;
        }
    }

    //----------------------------------------------------------------

    std::shared_ptr< mtk::Cluster_Group >
    Cell_Cluster::get_cluster_group( const moris_index aDiscretizationMeshIndex ) const
    {
        // check that the cluster group exists and is set
        MORIS_ASSERT( this->has_cluster_group( aDiscretizationMeshIndex ),
            "xtk::Cell_Cluster::get_cluster_group() - Cluster group is not set or does not exist." );

        // return the pointer to the cluster group
        return mClusterGroups( aDiscretizationMeshIndex ).lock();
    }

    //------------------------------------------------------------------------------

    void
    Cell_Cluster::set_cluster_group( 
            const moris_index aDiscretizationMeshIndex,
            std::shared_ptr< mtk::Cluster_Group > aClusterGroupPtr )
    {
        // check that the cluster group is set to the correct B-spline list index
        MORIS_ASSERT( aClusterGroupPtr->get_discretization_mesh_index_for_cluster_group() == aDiscretizationMeshIndex,
            "xtk::Cell_Cluster::set_cluster_group() - Index which the cluster group lives on is not the list index it gets set to on the cluster." );

        // check if the list of cluster groups is big enough to accommodate the cluster groups for each B-spline mesh 
        if( mClusterGroups.size() < (uint) aDiscretizationMeshIndex + 1 )
        {
            // ... if not increase the size
            mClusterGroups.resize( aDiscretizationMeshIndex + 1 );
        }
        
        // store pointer to the cluster group associated with this B-spline mesh
        mClusterGroups( aDiscretizationMeshIndex ) = aClusterGroupPtr;
    }

    //------------------------------------------------------------------------------

    std::shared_ptr<IG_Vertex_Group>
    Cell_Cluster::get_ig_vertex_group()
    {
        return mVertexGroup;
    }
    
    //------------------------------------------------------------------------------

    moris::real
    Cell_Cluster::compute_cluster_group_cell_measure(
            const moris_index aDiscretizationMeshIndex,
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // only compute // FIXME: ghost clusters for visualization should also have their cluster volumes assigned
        if( !mOnlyForVis ) 
        {
            // check that the cluster group exists and is set
            MORIS_ASSERT( this->has_cluster_group( aDiscretizationMeshIndex ),
                "xtk::Cell_Cluster::compute_cluster_group_cell_measure() - Cluster group is not set or does not exist." );

            // compute the group measure and return it
            return mClusterGroups( aDiscretizationMeshIndex ).lock()->compute_cluster_group_cell_measure( aPrimaryOrVoid, aIsLeader );
        }
        else // cluster groups are not defined on clusters that are only for visualization purposes (e.g. for ghost visualization)
        {        
            return -1.0;
        }
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell_Cluster::compute_cluster_group_cell_measure_derivative(
            const moris_index       aDiscretizationMeshIndex,
            const Matrix< DDRMat >& aPerturbedVertexCoords,
            uint aDirection,
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Leader_Follower aIsLeader ) const
    {
        // check that the cluster group exists and is set
        MORIS_ASSERT( this->has_cluster_group( aDiscretizationMeshIndex ),
            "xtk::Cell_Cluster::compute_cluster_group_cell_measure_derivative() - Cluster group is not set or does not exist." );

        // compute the group measure derivative and return it
        return mClusterGroups( aDiscretizationMeshIndex ).lock()->compute_cluster_group_cell_measure_derivative( aPerturbedVertexCoords, aDirection, aPrimaryOrVoid, aIsLeader );
    }

    //------------------------------------------------------------------------------
 
} // namespace xtk
// } // namespace moris

