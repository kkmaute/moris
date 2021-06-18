/*
 * cl_XTK_Cell_Cluster.cpp
 *
 *  Created on: Jul 22, 2019
 *      Author: doble
 */

#include "cl_XTK_Cell_Cluster.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Child_Mesh.hpp"
namespace xtk
{
    //----------------------------------------------------------------
    Cell_Cluster::Cell_Cluster():
            mTrivial(true),
            mInterpolationCell(nullptr),
            mChildMesh(nullptr),
            mPrimaryIntegrationCells(0,nullptr),
            mVoidIntegrationCells(0,nullptr),
            mVerticesInCluster(0,nullptr)
    {

    }
    //----------------------------------------------------------------
    Cell_Cluster::~Cell_Cluster()
    {

    }
    //----------------------------------------------------------------
    bool
    Cell_Cluster::is_trivial( const mtk::Master_Slave aIsMaster ) const
    {
        return mTrivial;
    }
    //----------------------------------------------------------------
    moris::Cell<moris::mtk::Cell const *> const &
    Cell_Cluster::get_primary_cells_in_cluster( const mtk::Master_Slave aIsMaster ) const
    {
        return mPrimaryIntegrationCells;
    }
    //----------------------------------------------------------------
    moris::Cell<moris::mtk::Cell const *> const &
    Cell_Cluster::get_void_cells_in_cluster() const
    {
        return mVoidIntegrationCells;
    }
    //----------------------------------------------------------------
    moris::mtk::Cell const &
    Cell_Cluster::get_interpolation_cell( const mtk::Master_Slave aIsMaster ) const
    {
        return *mInterpolationCell;
    }
    //----------------------------------------------------------------
    moris::Cell<moris::mtk::Vertex const *>
    Cell_Cluster::get_vertices_in_cluster( const mtk::Master_Slave aIsMaster ) const
    {
            return mVerticesInCluster;
    }
    //----------------------------------------------------------------
    moris::Matrix<moris::DDRMat>
    Cell_Cluster::get_vertices_local_coordinates_wrt_interp_cell( const mtk::Master_Slave aIsMaster )  const
    {
        if(!mTrivial)
        {
            return mChildMesh->get_parametric_coordinates();
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
            const mtk::Master_Slave aIsMaster ) const
    {
        MORIS_ERROR(!mTrivial,"Accessing local coordinates on a trivial cell cluster is not allowed");

        return mChildMesh->get_parametric_coordinates(aVertex->get_index());
    }
    //----------------------------------------------------------------
    moris_index
    Cell_Cluster::get_dim_of_param_coord( const mtk::Master_Slave aIsMaster) const
    {
        return this->get_vertices_local_coordinates_wrt_interp_cell(aIsMaster).n_cols();
    }

    moris::Matrix<moris::DDRMat>                  
    Cell_Cluster::get_primary_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aPrimaryCellClusterIndex) const
    {
        if(mTrivial)
        {
            return this->get_vertices_local_coordinates_wrt_interp_cell();
        }
        else
        {
             // MORIS_ERROR(!this->is_trivial(),"get_primary_cell_local_coords_on_side_wrt_interp_cell on trivial cluster is not allowed");
             MORIS_ASSERT(aPrimaryCellClusterIndex < (moris_index)this->get_num_primary_cells(),"Integration Cell Cluster index out of bounds");

            // get the integration cell of interest
            moris::mtk::Cell const * tIntegrationCell = this->get_primary_cells_in_cluster()(aPrimaryCellClusterIndex);

             // get the vertex pointers on the side
            moris::Cell<moris::mtk::Vertex *> tVerticesOnCell = tIntegrationCell->get_vertex_pointers();

            // allocate output (nnode x dim_xsi)
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
    Child_Mesh const *
    Cell_Cluster::get_xtk_child_mesh() const
    {
        MORIS_ASSERT(!this->is_trivial(),"Trivial xtk cell clusters do not have a child mesh.");

        return mChildMesh;
    }

    //----------------------------------------------------------------

    Matrix< IndexMat > Cell_Cluster::get_hanging_nodes(  ) const
    {
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
 
}


