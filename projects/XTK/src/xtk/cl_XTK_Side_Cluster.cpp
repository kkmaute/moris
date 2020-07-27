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
    Side_Cluster::get_cell_side_ordinal(moris::moris_index aCellIndexInCluster,
            const mtk::Master_Slave aIsMaster  ) const
    {
        MORIS_ASSERT(aCellIndexInCluster<(moris_index)mIntegrationCellSideOrdinals.numel(),"Cell index in cluster out of bounds");

        return mIntegrationCellSideOrdinals(aCellIndexInCluster);
    }
    //----------------------------------------------------------------
    moris::Cell<moris::mtk::Vertex const *> const &
    Side_Cluster::get_vertices_in_cluster( const mtk::Master_Slave aIsMaster ) const
    {
        return mVerticesInCluster;
    }
    //----------------------------------------------------------------
    moris::Matrix<moris::DDRMat> const &
    Side_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster )  const
    {
        MORIS_ERROR(!mTrivial,"Accessing local coordinates on a trivial cell cluster is not allowed");

        if(mChildMesh != nullptr)
        {
            return mChildMesh->get_parametric_coordinates();
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
        MORIS_ERROR(!mTrivial,"Accessing local coordinates on a trivial cell cluster is not allowed");

        if(mChildMesh != nullptr)
        {
            return mChildMesh->get_parametric_coordinates(aVertex->get_index());
        }
        else
        {
            return mVertexLocalCoords(this->get_vertex_cluster_index(aVertex));
        }
    }
    //----------------------------------------------------------------

    moris::moris_index
    Side_Cluster::get_vertex_cluster_index( moris::mtk::Vertex const * aVertex ) const
    {

        if(mTrivial ||  mChildMesh == nullptr)
        {
            auto tIter = mVertexIdToLocalIndex.find(aVertex->get_id());

            MORIS_ERROR(tIter != mVertexIdToLocalIndex.end(),"Vertex not found in side cluster");

            return tIter->second;
        }

        else
        {
            return mChildMesh->get_cm_local_node_index(aVertex->get_index());
        }
    }
    //----------------------------------------------------------------
    moris_index
    Side_Cluster::get_vertex_ordinal_on_facet( moris_index aCellIndexInCluster, moris::mtk::Vertex const * aVertex ) const
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
        return MORIS_INDEX_MAX;
    }
    //----------------------------------------------------------------
    moris_index
    Side_Cluster::get_dim_of_param_coord( const mtk::Master_Slave aIsMaster) const
    {
        return this->get_vertices_local_coordinates_wrt_interp_cell(aIsMaster).n_cols();
    }

    moris::real
    Side_Cluster::compute_cluster_cell_measure(const mtk::Primary_Void aPrimaryOrVoid,
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
    void
    Side_Cluster::finalize_setup()
    {
        moris::Cell<moris::mtk::Vertex const *> const & tVerticesInCluster = this->get_vertices_in_cluster();

        // add to map if trivial otherwise the child mesh takes care of this
        if(mTrivial || mChildMesh == nullptr)
        {
            for(moris::uint i = 0; i <tVerticesInCluster.size(); i++)
            {
                this->add_vertex_to_map(tVerticesInCluster(i)->get_id(),i);
            }
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
    Side_Cluster::add_vertex_to_map(moris_id aVertexId,
            moris_index aVertexLocalIndex)
    {
        MORIS_ERROR(mVertexIdToLocalIndex.find(aVertexId) == mVertexIdToLocalIndex.end(),"Trying to add vertex already found in side cluster");
        mVertexIdToLocalIndex[aVertexId] = aVertexLocalIndex;
    }

    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    //----------------------------------------------------------------

}


