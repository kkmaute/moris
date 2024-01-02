/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Cluster_DataBase.cpp
 *
 */

#include "cl_MTK_Cell_Cluster_DataBase.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "fn_TOL_Capacities.hpp"

namespace moris::mtk
{

    //------------------------------------------------------------------------------

    Cell_Cluster_DataBase::Cell_Cluster_DataBase( moris_index aCellClusterIndex,
        mtk::Mesh*                                            aMesh ) :
        mCellClusterIndex( aCellClusterIndex ),
        mMesh( aMesh )
    {
        // populate the data that is returned by reference
        this->set_outward_data();
    }

    //------------------------------------------------------------------------------

    bool
    Cell_Cluster_DataBase::is_trivial( const mtk::Leader_Follower aIsLeader ) const
    {
        return mMesh->cluster_is_trivial( ClusterType::CELL, mCellClusterIndex );
    }

    //------------------------------------------------------------------------------

    Vector< moris::mtk::Cell const* > const&
    Cell_Cluster_DataBase::get_primary_cells_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        return mPrimaryIntegrationCells;
    }

    //------------------------------------------------------------------------------

    Vector< moris::mtk::Cell const* > const&
    Cell_Cluster_DataBase::get_void_cells_in_cluster() const
    {
        return mVoidIntegrationCells;
    }

    //------------------------------------------------------------------------------

    moris::mtk::Cell const&
    Cell_Cluster_DataBase::get_interpolation_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        // defer the call to the mesh
        moris::mtk::Cell* tInterpolationCell = mMesh->get_ip_cell_in_cluster( ClusterType::CELL, mCellClusterIndex );

        // return the output
        return *tInterpolationCell;
    }

    //------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const* >
    Cell_Cluster_DataBase::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        // extract the vertices and number of them
        Vertex* const* tVertices  = mMesh->get_vertices_in_cluster( ClusterType::CELL, mCellClusterIndex );
        uint           tNumVertex = mMesh->get_num_vertices_in_cluster( ClusterType::CELL, mCellClusterIndex );

        // initialize the ouput
        Vector< moris::mtk::Vertex const* > tVerticesInCluster;

        // insert vertex pointers in the cell
        tVerticesInCluster.insert( 0, tVertices, tVertices + tNumVertex );

        // return the output
        return tVerticesInCluster;
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell_Cluster_DataBase::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        // determine if the cluster is trivial, coordinates for trivial is hard-coded and should not be stored
        bool tTrivial = mMesh->cluster_is_trivial( ClusterType::CELL, mCellClusterIndex );

        if ( !tTrivial )
        {
            // extract the necessary information coordinate matrix pointer, beginning of the numer, number of vertices  and spatial dim
            //  these data are to extract the right rows and columns of the coordinate matrix
            Matrix< DDRMat >* tVertexCoordsPtr = mMesh->get_local_coord_matrix_ptr( ClusterType::CELL, mCellClusterIndex );
            uint              tRowNumber       = mMesh->get_row_number_local_coords_matrix( ClusterType::CELL, mCellClusterIndex );
            uint              tNumVertices     = mMesh->get_num_vertices_in_cluster( ClusterType::CELL, mCellClusterIndex );
            uint              tSpatialDim      = mMesh->get_spatial_dim();

            // return the function
            return tVertexCoordsPtr->operator()( { tRowNumber, tRowNumber + tNumVertices - 1 }, { 0, tSpatialDim - 1 } );
        }
        else
        {
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const* tCellInfo = this->get_interpolation_cell().get_cell_info();

            // local coordinate matrix
            Matrix< DDRMat > tXi;

            // get the local coordinates on the side ordinal
            tCellInfo->get_loc_coords_of_cell( tXi );

            return tXi;
        }
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell_Cluster_DataBase::get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const* aVertex,
        const mtk::Leader_Follower                                                                   aIsLeader ) const
    {
        // get list of vertices and number of them
        Vertex* const* tVertices  = mMesh->get_vertices_in_cluster( ClusterType::CELL, mCellClusterIndex );
        uint           tNumVertex = mMesh->get_num_vertices_in_cluster( ClusterType::CELL, mCellClusterIndex );

        // lambda function to iterate through vertices and find the vertex ordinal
        auto tVertexOrdinalFinder = [aVertex]( mtk::Vertex* aVertices ) { return aVertices->get_index() == aVertex->get_index(); };

        // define the iterator to find the  ordinal
        auto itr = std::find_if( tVertices, tVertices + tNumVertex, tVertexOrdinalFinder );

        // check of the vertex exists
        if ( itr == tVertices + tNumVertex )
        {
            MORIS_ERROR( 0, " tried to access the coordinate of the vertex that does not belong the cluster" );
        }

        // convert iterator to ordinal
        moris_index tVertexOrdinal = std::distance( tVertices, itr );

        // ask for the coordinate matrix pointer and the start index of the coordinate block
        Matrix< DDRMat >* tVertexCoordsPtr = mMesh->get_local_coord_matrix_ptr( ClusterType::CELL, mCellClusterIndex );
        uint              tRowNumber       = mMesh->get_row_number_local_coords_matrix( ClusterType::CELL, mCellClusterIndex );

        // return the specific row
        return tVertexCoordsPtr->get_row( tRowNumber + tVertexOrdinal );
    }

    //------------------------------------------------------------------------------

    moris_index
    Cell_Cluster_DataBase::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader ) const
    {
        return (moris_index)mMesh->get_spatial_dim();
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell_Cluster_DataBase::get_primary_cell_local_coords_on_side_wrt_interp_cell( moris::moris_index aPrimaryCellClusterIndex ) const
    {
        // determine if the cluster is trivial
        bool tTrivial = mMesh->cluster_is_trivial( ClusterType::CELL, mCellClusterIndex );

        if ( tTrivial )
        {
            return this->get_vertices_local_coordinates_wrt_interp_cell();
        }
        else
        {
            // This part if the default!
            //  MORIS_ERROR(!this->is_trivial(),"get_primary_cell_local_coords_on_side_wrt_interp_cell on trivial cluster is not allowed");
            MORIS_ASSERT( aPrimaryCellClusterIndex < (moris_index)this->get_num_primary_cells(), "Integration Cell Cluster index out of bounds" );

            // get the integration cell of interest
            moris::mtk::Cell const* tIntegrationCell = this->get_primary_cells_in_cluster()( aPrimaryCellClusterIndex );

            // get the vertex pointers on the side
            Vector< moris::mtk::Vertex* > tVerticesOnCell = tIntegrationCell->get_vertex_pointers();

            // allocate output (n_node x dim_xsi)
            moris::Matrix< moris::DDRMat > tVertexParamCoords( tVerticesOnCell.size(), this->get_dim_of_param_coord() );

            // iterate through vertices and collect local coordinates
            for ( moris::uint i = 0; i < tVerticesOnCell.size(); i++ )
            {
                tVertexParamCoords.get_row( i ) = this->get_vertex_local_coordinate_wrt_interp_cell( tVerticesOnCell( i ) ).get_row( 0 );
            }

            return tVertexParamCoords;
        }
    }

    //------------------------------------------------------------------------------

    void
    Cell_Cluster_DataBase::set_outward_data()
    {
        //clear the data in the cells
        mPrimaryIntegrationCells.clear();
        mVoidIntegrationCells.clear();

        // get primary cell data
        mtk::Cell* const* tPrimaryCells    = mMesh->get_ig_cells_in_cluster( ClusterType::CELL, mtk::Primary_Void::PRIMARY, mCellClusterIndex );
        uint              tNumPrimaryCells = mMesh->get_num_cells_in_cluster( ClusterType::CELL, mtk::Primary_Void::PRIMARY, mCellClusterIndex );

        // populate the primary cell
        mPrimaryIntegrationCells.reserve( tNumPrimaryCells );
        mPrimaryIntegrationCells.insert( 0, tPrimaryCells, tPrimaryCells + tNumPrimaryCells );

        // get primary cell data
        mtk::Cell* const* tVoidCells    = mMesh->get_ig_cells_in_cluster( ClusterType::CELL, mtk::Primary_Void::VOID, mCellClusterIndex );
        uint              tNumVoidCells = mMesh->get_num_cells_in_cluster( ClusterType::CELL, mtk::Primary_Void::VOID, mCellClusterIndex );

        // populate the primary cell
        mVoidIntegrationCells.reserve( tNumVoidCells );
        mVoidIntegrationCells.insert( 0, tVoidCells, tVoidCells + tNumVoidCells );
    }

    //------------------------------------------------------------------------------

    size_t
    Cell_Cluster_DataBase::capacity()
    {
        // capacity of cell cluster
        size_t tCapacity = 0;

        // add up the member data size
        tCapacity += sizeof( mCellClusterIndex );
        tCapacity += sizeof( mMesh );
        tCapacity += mPrimaryIntegrationCells.capacity() * ( ( sizeof( void* ) ) + 1 );
        tCapacity += mVoidIntegrationCells.capacity() * ( ( sizeof( void* ) ) + 1 );

        return tCapacity;
    }

}// namespace moris::mtk

