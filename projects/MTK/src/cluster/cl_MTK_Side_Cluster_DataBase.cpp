/*
 * Copyright (c) 2022 University of Colorado 
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details. 
 * 
 * ------------------------------------------------------------------------------------ 
 * 
 * cl_MTK_Side_Cluster_DataBase.cpp  
 * 
 */

#include "cl_MTK_Side_Cluster_DataBase.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "fn_TOL_Capacities.hpp"

namespace moris::mtk
{

    // ----------------------------------------------------------------------------------

    Side_Cluster_DataBase::Side_Cluster_DataBase( moris_index aSideClusterIndex,
        mtk::Mesh*                                            aMesh ) :
        mSideClusterIndex( aSideClusterIndex ),
        mMesh( aMesh )
    {
        this->set_outward_data();
    }

    // ----------------------------------------------------------------------------------

    bool
    Side_Cluster_DataBase::is_trivial( const mtk::Leader_Follower aIsLeader ) const
    {
        return mMesh->cluster_is_trivial( ClusterType::SIDE, mSideClusterIndex );
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Cell const&
    Side_Cluster_DataBase::get_interpolation_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        moris::mtk::Cell* tInterpolationCell = mMesh->get_ip_cell_in_cluster( ClusterType::SIDE, mSideClusterIndex );
        return *tInterpolationCell;
    }

    // ----------------------------------------------------------------------------------

    Vector< mtk::Cell const* > const&
    Side_Cluster_DataBase::get_cells_in_side_cluster() const
    {
        return mPrimaryIntegrationCells;
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::IndexMat >
    Side_Cluster_DataBase::get_cell_side_ordinals( const mtk::Leader_Follower aIsLeader ) const
    {
        // get the first side ord position and number of cells
        moris_index* tSideOrdinalPtr  = mMesh->get_side_ordinals_in_cluster( ClusterType::SIDE, mSideClusterIndex );
        uint         tNumPrimaryCells = mMesh->get_num_cells_in_cluster( ClusterType::SIDE, mtk::Primary_Void::PRIMARY, mSideClusterIndex );

        // advance constrcut the matrix
        moris::Matrix< moris::IndexMat > tIntegrationCellSideOrdinals( tSideOrdinalPtr, 1, tNumPrimaryCells, false, true );

        return tIntegrationCellSideOrdinals;
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Side_Cluster_DataBase::get_cell_side_ordinal( moris::moris_index aCellIndexInCluster,
        const mtk::Leader_Follower                                      aIsLeader ) const
    {
        // get the first side ordinal position
        moris_index* tSideOrdinalPtr = mMesh->get_side_ordinals_in_cluster( ClusterType::SIDE, mSideClusterIndex );

        return *( tSideOrdinalPtr + aCellIndexInCluster );
    }

    // ----------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const* >
    Side_Cluster_DataBase::get_vertices_in_cluster( const mtk::Leader_Follower aIsLeader ) const
    {
        // initialize the output
        Vector< moris::mtk::Vertex const* > tVerticesInCluster;

        // get vertex pointer position and number of vertices
        Vertex* const* tVertices  = mMesh->get_vertices_in_cluster( ClusterType::SIDE, mSideClusterIndex );
        uint           tNumVertex = mMesh->get_num_vertices_in_cluster( ClusterType::SIDE, mSideClusterIndex );

        // Insert pointers in the output cell
        tVerticesInCluster.insert( 0, tVertices, tVertices + tNumVertex );

        return tVerticesInCluster;
    }


    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Side_Cluster_DataBase::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Leader_Follower aIsLeader ) const
    {
        // determine if the cluster is trivial
        bool tTrivial = mMesh->cluster_is_trivial( ClusterType::SIDE, mSideClusterIndex );

        // if it is not trivial
        if ( !tTrivial )
        {
            // get the matrix pointer and row index
            Matrix< DDRMat >* tVertexCoordsPtr = mMesh->get_local_coord_matrix_ptr( ClusterType::SIDE, mSideClusterIndex );
            uint              tRowNumber       = mMesh->get_row_number_local_coords_matrix( ClusterType::SIDE, mSideClusterIndex );

            // if it is a not a ghost
            if ( !mMesh->is_secondary_cluster( mSideClusterIndex ) )
            {
                // get number of vertices in the cell cluster
                uint tNumVertices = mMesh->get_associated_cell_cluster( mSideClusterIndex )->get_num_vertices_in_cluster();
                uint tSpatialDim  = mMesh->get_spatial_dim();

                return tVertexCoordsPtr->operator()( { tRowNumber, tRowNumber + tNumVertices - 1 }, { 0, tSpatialDim - 1 } );
            }

            // if it is a ghost
            else
            {
                uint tNumVertices = this->get_num_vertices_in_cluster();
                uint tSpatialDim  = mMesh->get_spatial_dim();

                return tVertexCoordsPtr->operator()( { tRowNumber, tRowNumber + tNumVertices - 1 }, { 0, tSpatialDim - 1 } );
            }
        }
        else
        {
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const* tCellInfo = this->get_interpolation_cell().get_cell_info();

            moris_index* tSideOrdinalPtr = mMesh->get_side_ordinals_in_cluster( ClusterType::SIDE, mSideClusterIndex );

            // side ordinal on interpolation cell
            moris::uint tSideOrd = (uint)( *tSideOrdinalPtr );

            // local coordinate matrix
            Matrix< DDRMat > tXi;

            // get the local coordinates on the side ordinal
            tCellInfo->get_loc_coord_on_side_ordinal( tSideOrd, tXi );

            return tXi;
        }
    }

    // ----------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Side_Cluster_DataBase::get_vertex_local_coordinate_wrt_interp_cell( 
        moris::mtk::Vertex const* aVertex,
        const mtk::Leader_Follower   aIsLeader ) const
    {
        // determine if the cluster is trivial
        bool tTrivial = mMesh->cluster_is_trivial( ClusterType::SIDE, mSideClusterIndex );

        // if it is not trivial
        if ( !tTrivial )
        {
            // if it is a not a ghost
             if ( !mMesh->is_secondary_cluster( mSideClusterIndex ) )
            {
                mtk::Cell_Cluster const* tAssociatedCellCluster = mMesh->get_associated_cell_cluster( mSideClusterIndex );

                return tAssociatedCellCluster->get_vertex_local_coordinate_wrt_interp_cell( aVertex, aIsLeader );
            }

            // if it is a ghost
            else
            {
                Vertex* const* tVertices  = mMesh->get_vertices_in_cluster( ClusterType::SIDE, mSideClusterIndex );
                uint           tNumVertex = mMesh->get_num_vertices_in_cluster( ClusterType::SIDE, mSideClusterIndex );

                auto tVertexOrdinalFinder = [aVertex]( mtk::Vertex* aVertices ) { return aVertices->get_index() == aVertex->get_index(); };

                auto itr = std::find_if( tVertices, tVertices + tNumVertex, tVertexOrdinalFinder );

                if ( itr == tVertices + tNumVertex )
                {
                    MORIS_ERROR( 0, " tried to access the coordinate of the vertex that does not belong the cluster" );
                }

                moris_index tVertexOrdinal = std::distance( tVertices, itr );

                Matrix< DDRMat >* tVertexCoordsPtr = mMesh->get_local_coord_matrix_ptr( ClusterType::SIDE, mSideClusterIndex );
                uint              tRowNumber       = mMesh->get_row_number_local_coords_matrix( ClusterType::SIDE, mSideClusterIndex );

                return tVertexCoordsPtr->get_row( tRowNumber + tVertexOrdinal );
            }
        }

        // if it is trivial
        else
        {
            mtk::Cell* const* tPrimaryCells = mMesh->get_ig_cells_in_cluster( ClusterType::SIDE, mtk::Primary_Void::PRIMARY, mSideClusterIndex );

            // Cluster index is also the ordinal in trivial cluster
            moris_index tVertexOrdinal = tPrimaryCells[0]->get_vertex_ordinal_wrt_cell( aVertex->get_index() );

            // std::cout<<"XTK Ord = "<<tVertexOrdinal<<std::endl;
            //  get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const* tCellInfo = this->get_interpolation_cell().get_cell_info();

            // get the local coordinates on the side ordinal
            Matrix< DDRMat > tXi = tCellInfo->get_vertex_loc_coord( tVertexOrdinal );

            return tXi;
        }
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Side_Cluster_DataBase::get_dim_of_param_coord( const mtk::Leader_Follower aIsLeader ) const
    {
        return (moris_index)mMesh->get_spatial_dim();
    }

    // ----------------------------------------------------------------------------------

    void
    Side_Cluster_DataBase::set_outward_data()
    {
        mPrimaryIntegrationCells.clear();

        mtk::Cell* const* tPrimaryCells    = mMesh->get_ig_cells_in_cluster( ClusterType::SIDE, mtk::Primary_Void::PRIMARY, mSideClusterIndex );
        uint              tNumPrimaryCells = mMesh->get_num_cells_in_cluster( ClusterType::SIDE, mtk::Primary_Void::PRIMARY, mSideClusterIndex );

        mPrimaryIntegrationCells.reserve( tNumPrimaryCells );
        mPrimaryIntegrationCells.insert( 0, tPrimaryCells, tPrimaryCells + tNumPrimaryCells );

        // mVoidIntegrationCells.resize( mNumVoidCells );
        //  mVoidIntegrationCells.insert( 0, mVoidCells, mVoidCells + mNumVoidCells );
    }

    // ----------------------------------------------------------------------------------

    moris_index
    Side_Cluster_DataBase::get_vertex_cluster_index( const Vertex* aVertex,
        const mtk::Leader_Follower                                    aIsLeader ) const
    {
        auto tVertexOrdinalFinder = [aVertex]( mtk::Vertex* aVertices ) { return aVertices->get_index() == aVertex->get_index(); };

        Vertex* const* tVertices = mMesh->get_vertices_in_cluster( ClusterType::SIDE, mSideClusterIndex );

        uint tNumVertex = mMesh->get_num_vertices_in_cluster( ClusterType::SIDE, mSideClusterIndex );

        auto itr = std::find_if( tVertices, tVertices + tNumVertex, tVertexOrdinalFinder );

        moris_index tVertexOrdinal = std::distance( tVertices, itr );

        return tVertexOrdinal;
    }

    //-----------------------------------------------
    moris_index
    Side_Cluster_DataBase::get_vertex_ordinal_on_facet(
        moris_index               aCellIndexInCluster,
        moris::mtk::Vertex const* aVertex ) const
    {
        mtk::Cell* const* tPrimaryCells = mMesh->get_ig_cells_in_cluster( ClusterType::SIDE, mtk::Primary_Void::PRIMARY, mSideClusterIndex );

        moris_index* tSideOrdinalPtr = mMesh->get_side_ordinals_in_cluster( ClusterType::SIDE, mSideClusterIndex );

        moris_index tSideOrd = tSideOrdinalPtr[aCellIndexInCluster];

        Vector< moris::mtk::Vertex const* > tVerticesOnSide = tPrimaryCells[aCellIndexInCluster]->get_vertices_on_side_ordinal( tSideOrd );

        // iterate through vertices and see if the ids match
        for ( moris::moris_index i = 0; i < (moris_index)tVerticesOnSide.size(); i++ )
        {
            if ( tVerticesOnSide( i )->get_id() == aVertex->get_id() )
            {
                return i;
            }
        }

        // print debug information if search for vertex on side cluster fails
        std::cout << "Side_Cluster_DataBase::get_vertex_ordinal_on_facet() - " << 
            "Looking for vertex with ID " << aVertex->get_id() << " in Dbl-Side Cluster with Vertices [ ";
        for ( moris::moris_index i = 0; i < (moris_index)tVerticesOnSide.size(); i++ )
        {
            std::cout << tVerticesOnSide( i )->get_id() << " ";
        }
        std::cout << "]" << std::endl;

        // throw error and return default
        MORIS_ERROR( false, "Side_Cluster_DataBase::get_vertex_ordinal_on_facet() - Vertex not found on facet." );
        return MORIS_INDEX_MAX;
    }


    //----------------------------------------------------------------

    moris::real
    Side_Cluster_DataBase::compute_cluster_cell_measure(
        const mtk::Primary_Void aPrimaryOrVoid,
        const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY || aPrimaryOrVoid == mtk::Primary_Void::VOID )
        {
            mtk::Cell_Cluster const* tAssociatedCellCluster = mMesh->get_associated_cell_cluster( mSideClusterIndex );

            return tAssociatedCellCluster->compute_cluster_cell_measure( aPrimaryOrVoid, aIsLeader );
        }
        else
        {
            return this->get_interpolation_cell().compute_cell_measure();
        }
    }

    //----------------------------------------------------------------

    Matrix< DDRMat >
    Side_Cluster_DataBase::compute_cluster_ig_cell_measures(
        const mtk::Primary_Void aPrimaryOrVoid,
        const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY || aPrimaryOrVoid == mtk::Primary_Void::VOID )
        {
            mtk::Cell_Cluster const* tAssociatedCellCluster = mMesh->get_associated_cell_cluster( mSideClusterIndex );

            return tAssociatedCellCluster->compute_cluster_ig_cell_measures( aPrimaryOrVoid, aIsLeader );
        }
        else
        {
            return { { this->get_interpolation_cell().compute_cell_measure() } };
        }
    }

    //----------------------------------------------------------------

    moris::real
    Side_Cluster_DataBase::compute_cluster_cell_measure_derivative(
        const Matrix< DDRMat >& aPerturbedVertexCoords,
        uint                    aDirection,
        const mtk::Primary_Void aPrimaryOrVoid,
        const mtk::Leader_Follower aIsLeader ) const
    {
        if ( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY || aPrimaryOrVoid == mtk::Primary_Void::VOID )
        {
            mtk::Cell_Cluster const* tAssociatedCellCluster = mMesh->get_associated_cell_cluster( mSideClusterIndex );

            return tAssociatedCellCluster->compute_cluster_cell_measure_derivative(
                aPerturbedVertexCoords,
                aDirection,
                aPrimaryOrVoid,
                aIsLeader );
        }
        else
        {
            MORIS_LOG( "mInterpolationCell->compute_cell_measure_derivative() is set to zero" );
            return 0.0;
        }
    }

    //----------------------------------------------------------------

    size_t
    Side_Cluster_DataBase::capacity()
    {
        size_t tCapacity = 0;

        tCapacity += sizeof( mSideClusterIndex );
        tCapacity += sizeof( mMesh );
        tCapacity += mPrimaryIntegrationCells.capacity() * ( ( sizeof( void* ) ) + 1 );
        tCapacity += mVoidIntegrationCells.capacity() * ( ( sizeof( void* ) ) + 1 );

        return tCapacity;
    }

    //----------------------------------------------------------------
    void
    Side_Cluster_DataBase::update_cluster_index( moris_index aNewIndex )
    {
        mSideClusterIndex = aNewIndex;
    }
}// namespace moris::mtk
