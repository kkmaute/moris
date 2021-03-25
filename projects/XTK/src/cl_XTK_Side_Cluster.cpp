/*
 * cl_XTK_Side_Cluster.cpp
 *
 *  Created on: Jul 22, 2019
 *      Author: doble
 */

#include "cl_XTK_Side_Cluster.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
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
        if(mChildMesh != nullptr)
        {
            return mChildMesh->get_parametric_coordinates();
        }
        else if(mTrivial)
        {
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const * tCellInfo = mInterpolationCell->get_connectivity();

            // side ordinal on interpolation cell
            moris::uint tSideOrd = (uint) mIntegrationCellSideOrdinals(0);

            // local coordinate matrix
            Matrix<DDRMat> tXi;

            // get the local coordinates on the side ordinal
            tCellInfo->get_loc_coord_on_side_ordinal(tSideOrd,tXi);

            return tXi;
        }
        else
        {
            return mVertexLocalCoordsMat;
        }
    }
    
    //----------------------------------------------------------------
    
    moris::Matrix<moris::DDRMat>
    Side_Cluster::get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
            const mtk::Master_Slave aIsMaster ) const
    {
        if(mChildMesh != nullptr)
        {
            return mChildMesh->get_parametric_coordinates(aVertex->get_index());
        }
        else if(mTrivial)
        {
            // Cluster index is also the ordinal in trivial cluster
            moris_index tVertexOrdinal = mIntegrationCells(0)->get_vertex_ordinal_wrt_cell(aVertex->get_index());

            //std::cout<<"XTK Ord = "<<tVertexOrdinal<<std::endl;
            // get the interpolation cell's connectivity information
            moris::mtk::Cell_Info const * tCellInfo = mInterpolationCell->get_connectivity();

            // get the local coordinates on the side ordinal
            Matrix<DDRMat> tXi  = tCellInfo->get_vertex_loc_coord(tVertexOrdinal);

            return tXi;
        }
        else
        {
            return mVertexLocalCoords(this->get_vertex_cluster_index(aVertex));
        }
    }
    
    //----------------------------------------------------------------

    moris::moris_index
    Side_Cluster::get_vertex_cluster_index( moris::mtk::Vertex const * aVertex,
                                            const mtk::Master_Slave aIsMaster) const
    {

        // if(mTrivial ||  mChildMesh == nullptr)
        // {
            auto tIter = mVertexIdToLocalIndex.find(aVertex->get_id());

            MORIS_ERROR(tIter != mVertexIdToLocalIndex.end(),"Vertex not found in side cluster");

            return tIter->second;
        // }

        // else
        // {
            // return mChildMesh->get_cm_local_node_index(aVertex->get_index());
        // }
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
            MORIS_ERROR( mVertexIdToLocalIndex.find(aVertex->get_id()) != mVertexIdToLocalIndex.end(),"Vertex not found in side cluster");

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
            return mAssociatedCellCluster->compute_cluster_cell_measure(aPrimaryOrVoid,aIsMaster);
        }
        else
        {
            return mInterpolationCell->compute_cell_measure();
        }
    }
    
    //----------------------------------------------------------------
    
    void
    Side_Cluster::print_vertex_map() const
    {
        for(auto i :mVertexIdToLocalIndex)
        {
            std::cout<<"Vertex Id: "<<i.first<<" | Cluster Index: "<<i.second<< std::endl;
        }
    }
    
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

    void
    Side_Cluster::finalize_setup()
    {
        moris::Cell<moris::mtk::Vertex const *> tVerticesInCluster = this->get_vertices_in_cluster();


        for(moris::uint i = 0; i <tVerticesInCluster.size(); i++)
        {
            this->add_vertex_to_map(tVerticesInCluster(i)->get_id(),i);
        }

        if( mChildMesh == nullptr && !mTrivial )
        {
            for(moris::uint i = 0; i <tVerticesInCluster.size(); i++)
            {
                if( i == 0)
                {
                    mVertexLocalCoordsMat.resize(tVerticesInCluster.size(),mVertexLocalCoords(i).n_cols());
                }

                mVertexLocalCoordsMat.get_row(i) = mVertexLocalCoords(i).matrix_data();
            }
        }

    }
    
    //----------------------------------------------------------------
    
    void
    Side_Cluster::add_vertex_to_map(
            moris_id    aVertexId,
            moris_index aVertexLocalIndex)
    {
        MORIS_ERROR(mVertexIdToLocalIndex.find(aVertexId) == mVertexIdToLocalIndex.end(),"Trying to add vertex already found in side cluster");
        mVertexIdToLocalIndex[aVertexId] = aVertexLocalIndex;
    }

    //----------------------------------------------------------------
}


