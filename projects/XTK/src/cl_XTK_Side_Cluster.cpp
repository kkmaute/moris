/*
 * cl_XTK_Side_Cluster.cpp
 *
 *  Created on: Jul 22, 2019
 *      Author: doble
 */

#include "cl_XTK_Side_Cluster.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "assert.hpp"
#include "fn_TOL_Capacities.hpp"



namespace xtk
{
    Side_Cluster::Side_Cluster():
                    mTrivial(true),
                    mInterpolationCell(nullptr),
                    mChildMesh(nullptr),
                    mIntegrationCells(0,nullptr),
                    mIntegrationCellSideOrdinals(0,0),
                    mVerticesInCluster(0,nullptr),
                    mVertexLocalCoords(0)
    {}
    
    //----------------------------------------------------------------
    
    bool
    Side_Cluster::is_trivial( const mtk::Master_Slave aIsMaster ) const
    {
        return mTrivial;
    }
    
    //----------------------------------------------------------------
    
    moris::mtk::Cell const &
    Side_Cluster::get_interpolation_cell(const mtk::Master_Slave aIsMaster ) const
    {
        return *mInterpolationCell;
    }
    
    //----------------------------------------------------------------
    
    moris::Cell<moris::mtk::Cell const *> const &
    Side_Cluster::get_cells_in_side_cluster() const
    {
        return mIntegrationCells;
    }
    
    //----------------------------------------------------------------
    
    moris::Matrix<moris::IndexMat>
    Side_Cluster::get_cell_side_ordinals( const mtk::Master_Slave aIsMaster ) const
    {
        return mIntegrationCellSideOrdinals;
    }
    
    //----------------------------------------------------------------
    
    moris_index
    Side_Cluster::get_cell_side_ordinal(
            moris::moris_index      aCellIndexInCluster,
            const mtk::Master_Slave aIsMaster  ) const
    {
        MORIS_ASSERT(aCellIndexInCluster<(moris_index)mIntegrationCellSideOrdinals.numel(),"Cell index in cluster out of bounds");

        return mIntegrationCellSideOrdinals(aCellIndexInCluster);
    }
    
    //----------------------------------------------------------------
    
    moris::Cell<moris::mtk::Vertex const *>
    Side_Cluster::get_vertices_in_cluster( const mtk::Master_Slave aIsMaster ) const
    {
        return mVerticesInCluster;
    }
    
    //----------------------------------------------------------------
    
    moris::Matrix<moris::DDRMat>
    Side_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster )  const
    {
        if(!mTrivial)
        {
            return mVertexLocalCoordsMat;
        }
        else
        {
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const * tCellInfo = mInterpolationCell->get_cell_info();

            // side ordinal on interpolation cell
            moris::uint tSideOrd = (uint) mIntegrationCellSideOrdinals(0);

            // local coordinate matrix
            Matrix<DDRMat> tXi;

            // get the local coordinates on the side ordinal
            tCellInfo->get_loc_coord_on_side_ordinal(tSideOrd,tXi);

            return tXi;
        }
    }
    
    //----------------------------------------------------------------
    
    moris::Matrix<moris::DDRMat>
    Side_Cluster::get_vertex_local_coordinate_wrt_interp_cell(
            moris::mtk::Vertex const * aVertex,
            const mtk::Master_Slave    aIsMaster ) const
    {
        if(!mTrivial)
        {
            return *mVertexGroup->get_vertex_local_coords(aVertex->get_index());
        }
        else
        {
            // Cluster index is also the ordinal in trivial cluster
            moris_index tVertexOrdinal = mIntegrationCells(0)->get_vertex_ordinal_wrt_cell(aVertex->get_index());

            //std::cout<<"XTK Ord = "<<tVertexOrdinal<<std::endl;
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const * tCellInfo = mInterpolationCell->get_cell_info();

            // get the local coordinates on the side ordinal
            Matrix<DDRMat> tXi  = tCellInfo->get_vertex_loc_coord(tVertexOrdinal);

            return tXi;
        }
    }
    
    //----------------------------------------------------------------

    moris::moris_index
    Side_Cluster::get_vertex_cluster_index( moris::mtk::Vertex const * aVertex,
                                            const mtk::Master_Slave aIsMaster) const
    {
        return mVertexGroup->get_vertex_group_ordinal(aVertex->get_index());
    }
    
    
    //----------------------------------------------------------------
    
    moris_index
    Side_Cluster::get_vertex_ordinal_on_facet( 
            moris_index                aCellIndexInCluster, 
            moris::mtk::Vertex const * aVertex ) const
    {
        moris_index tSideOrd = mIntegrationCellSideOrdinals(aCellIndexInCluster);


        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = mIntegrationCells(aCellIndexInCluster)->get_vertices_on_side_ordinal(tSideOrd);
        // iterate through vertices and see if the ids match
        for(moris::moris_index i = 0; i < (moris_index)tVerticesOnSide.size(); i++)
        {
            
            if(tVerticesOnSide(i)->get_id() == aVertex->get_id())
            {
                return i;
            }
        }

        MORIS_ERROR(0,"Vertex not found on facet.");
        return MORIS_INDEX_MAX;
    }
    
    //----------------------------------------------------------------
    
    moris_index
    Side_Cluster::get_dim_of_param_coord( const mtk::Master_Slave aIsMaster) const
    {
        return this->get_vertices_local_coordinates_wrt_interp_cell(aIsMaster).n_cols();
    }

    //----------------------------------------------------------------
    
    moris::real
    Side_Cluster::compute_cluster_cell_measure(
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Master_Slave aIsMaster) const
    {
        if( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY ||  aPrimaryOrVoid == mtk::Primary_Void::VOID )
        {
            MORIS_ASSERT( mAssociatedCellCluster,
                    "Side_Cluster::compute_cluster_cell_measure - Associated cell cluster not set.\n");

            return mAssociatedCellCluster->compute_cluster_cell_measure(aPrimaryOrVoid,aIsMaster);
        }
        else
        {
            MORIS_ASSERT( mInterpolationCell,
                    "Side_Cluster::compute_cluster_cell_measure - Interpolation cell not set.\n");

            return mInterpolationCell->compute_cell_measure();
        }
    }

    //----------------------------------------------------------------

    Matrix<DDRMat>
    Side_Cluster::compute_cluster_ig_cell_measures(
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Master_Slave aIsMaster) const
    {
        if( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY ||  aPrimaryOrVoid == mtk::Primary_Void::VOID )
        {
            MORIS_ASSERT( mAssociatedCellCluster,
                    "Side_Cluster::compute_cluster_ig_cell_measures - Associated cell cluster not set.\n");

            return mAssociatedCellCluster->compute_cluster_ig_cell_measures(aPrimaryOrVoid,aIsMaster);
        }
        else
        {
            MORIS_ASSERT( mInterpolationCell,
                    "Side_Cluster::compute_cluster_ig_cell_measures - Interpolation cell not set.\n");

            return {{mInterpolationCell->compute_cell_measure()}};
        }
    }
    
    //----------------------------------------------------------------
    
    moris::real
    Side_Cluster::compute_cluster_cell_measure_derivative(
            const Matrix< DDRMat > & aPerturbedVertexCoords,
            uint aDirection,
            const mtk::Primary_Void aPrimaryOrVoid,
            const mtk::Master_Slave aIsMaster) const
    {
        if( aPrimaryOrVoid == mtk::Primary_Void::PRIMARY ||  aPrimaryOrVoid == mtk::Primary_Void::VOID )
        {
            return mAssociatedCellCluster->compute_cluster_cell_measure_derivative(
                    aPerturbedVertexCoords,
                    aDirection,
                    aPrimaryOrVoid,
                    aIsMaster);
        }
        else
        {
            MORIS_LOG( "mInterpolationCell->compute_cell_measure_derivative() is set to zero" );
            return 0.0;
        }
    }

    void
    Side_Cluster::set_ig_vertex_group(std::shared_ptr<IG_Vertex_Group> aVertexGroup)
    {
        mVertexGroup = aVertexGroup;

        mVerticesInCluster.resize(mVertexGroup->size());
        mVertexLocalCoordsMat.resize(mVertexGroup->size(),mVertexGroup->get_vertex_local_coords(aVertexGroup->get_vertex(0)->get_index())->n_cols());

        for(moris::uint i = 0; i < mVertexGroup->size(); i++)
        {
            mVerticesInCluster(i) = mVertexGroup->get_vertex(i);
            mVertexLocalCoordsMat.set_row(i, *mVertexGroup->get_vertex_local_coords(mVerticesInCluster(i)->get_index()));
        }
    }

    //----------------------------------------------------------------

    
    //----------------------------------------------------------------

    size_t
    Side_Cluster::capacity()
    {

        size_t tTotalSize = 0;
        tTotalSize += sizeof(mTrivial);
        tTotalSize += sizeof(mInterpolationCell);
        tTotalSize += sizeof(mChildMesh);
        tTotalSize += mIntegrationCells.capacity();
        tTotalSize += mIntegrationCellSideOrdinals.capacity();

        tTotalSize += mVerticesInCluster.capacity();
        tTotalSize += moris::internal_capacity(mVertexLocalCoords);
        tTotalSize += sizeof(mAssociatedCellCluster);

        return tTotalSize;
    }

    
    //----------------------------------------------------------------
    

    //----------------------------------------------------------------
}


