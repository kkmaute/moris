/*
 * cl_MTK_Cluster.hpp
 *
 *  Created on: Apr 26, 2019
 *      Author: Schmidt
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
namespace moris
{
namespace mtk
{
class Cluster
{
private:
    moris::Cell<moris::mtk::Cell const *>   mDummCellCell;

public:
    Cluster(){};

    //##############################################
    // Characteristic functions
    //##############################################

    virtual bool is_trivial( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;

    //##############################################
    // Cell/Vertex Access
    //##############################################

    virtual
    moris_index
    get_vertex_cluster_index( const Vertex * aVertex,
                              const moris::uint aSide = 0 ) const
    {
        MORIS_ERROR(false, "get_vertex_cluster_index(): not implemented for this cluster type");
        return 0;
    }

    virtual
    moris::Cell<moris::mtk::Cell const *> const &
    get_primary_cells_in_cluster( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;

    virtual
    moris::Cell<moris::mtk::Cell const *> const &
    get_void_cells_in_cluster() const
    {
        MORIS_ERROR(false, "get_void_cells_in_cluster(): not implemented for this cluster type");
        return mDummCellCell;
    }

    virtual
    moris::mtk::Cell const &
    get_interpolation_cell( const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER ) const = 0;

    virtual
    moris::Matrix<moris::IndexMat>
    get_cell_side_ordinals( const moris::uint aSide = 0 ) const
    {
        MORIS_ERROR(false, "get_interpolation_cell(): not implemented for this cluster type");
        return moris::Matrix<moris::IndexMat>(0,0);
    };

    virtual
    moris_index
    get_cell_side_ordinal(moris::moris_index aCellIndexInCluster,
                          const moris::uint aSide = 0) const
    {
        MORIS_ERROR(false, "get_cell_side_ordinal(): not implemented for this cluster type");
        return 0;
    };

    virtual
    moris::Cell<moris::mtk::Vertex const *> const &
    get_vertices_in_cluster( const moris::uint aSide = 0 ) const = 0;

    virtual
    moris::mtk::Vertex const *
    get_left_vertex_pair(moris::mtk::Vertex const * aLeftVertex) const
    {
        MORIS_ERROR(false, "get_left_vertex_pair(): not implemented for this cluster type");
        return nullptr;
    }

    //##############################################
    // Local Coordinate Access
    // (Pure Virtual)
    //##############################################
    virtual
    moris::Matrix<moris::DDRMat> const &
    get_vertices_local_coordinates_wrt_interp_cell( const moris::uint aSide = 0 ) const = 0;

    /*
     * Access a single local coordinate of a vertex
     */
    virtual
    moris::Matrix<moris::DDRMat>
    get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex,
                                                 const moris::uint aSide = 0 ) const = 0;

    virtual
    moris::Matrix<moris::DDRMat>
    get_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aLeftClusterLocalIndex,
                                                  const moris::uint aSide = 0 ) const
    {
        MORIS_ERROR(false, "get_cell_local_coords_on_side_wrt_interp_cell(): not implemented for this cluster type");
        return moris::Matrix<moris::DDRMat>(0,0);
    }

    //##############################################
    // Size Access
    // (Pure Virtual)
    //##############################################
    /*!
     * Size of the xsi vector in this side cluster
     */
    virtual
    moris_index
    get_dim_of_param_coord( const moris::uint aSide = 0) const = 0;

    // ---------------------------------------------
    // EVERYTHING BELOW THIS LINE HAS A DEFAULT
    // IMPLEMENTATION
    // ---------------------------------------------

    //##############################################
    // Cell/Vertex Index Access
    //##############################################
    virtual
    moris::Matrix<moris::IndexMat>
    get_primary_cell_indices_in_cluster() const
    {
        MORIS_ERROR(false, "get_primary_cell_indices_in_cluster(): not implemented for this cluster type");
        return moris::Matrix<moris::IndexMat>(0,0);
    }

    virtual
    moris::Matrix<moris::IndexMat>
    get_void_cell_indices_in_cluster() const
    {
        MORIS_ERROR(false, "get_void_cell_indices_in_cluster(): not implemented for this cluster type");
        return moris::Matrix<moris::IndexMat>(0,0);
    }

    virtual
    moris::moris_index
    get_interpolation_cell_index() const
    {
        MORIS_ERROR(false, "get_interpolation_cell_index(): not implemented for this cluster type");
        return moris::moris_index(0);
    }

    virtual
    moris::Matrix<moris::IndexMat>
    get_vertex_indices_in_cluster() const
    {
        MORIS_ERROR(false, "get_vertex_indices_in_cluster(): not implemented for this cluster type");
        return moris::Matrix<moris::IndexMat>(0,0);
    }

    //##############################################
    // Cell/Vertex Id Access
    //##############################################
    virtual
    moris::Matrix<moris::IdMat>
    get_primary_cell_ids_in_cluster() const
    {
        MORIS_ERROR(false, "get_primary_cell_ids_in_cluster(): not implemented for this cluster type");
        return moris::Matrix<moris::IndexMat>(0,0);
    }

    virtual
    moris::Matrix<moris::IdMat>
    get_void_cell_ids_in_cluster() const
    {
        MORIS_ERROR(false, "get_void_cell_ids_in_cluster(): not implemented for this cluster type");
        return moris::Matrix<moris::IdMat>(0,0);
    }

    virtual
    moris::moris_id
    get_interpolation_cell_id() const
    {
        MORIS_ERROR(false, "get_interpolation_cell_id(): not implemented for this cluster type");
        return 0;
    }

    virtual
    moris::Matrix<moris::IdMat>
    get_vertex_ids_in_cluster() const
    {
        MORIS_ERROR(false, "get_vertex_ids_in_cluster(): not implemented for this cluster type");
        return moris::Matrix<moris::IdMat>(0,0);
    }

    //##############################################
    // Local Coordinate access
    //##############################################

    /*!
     * Access a primary integration cells parametric coordinates relative to the interpolation cell
     * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
     */
    virtual
    moris::Matrix<moris::DDRMat>
    get_primary_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aPrimaryCellClusterIndex) const
    {
        MORIS_ERROR(false, "get_primary_cell_local_coords_on_side_wrt_interp_cell(): not implemented for this cluster type");
        return moris::Matrix<moris::DDRMat>(0,0);
    }

    /*!
     * Access a void integration cells parametric coordinates relative to the interpolation cell
     * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
     */
    virtual
    moris::Matrix<moris::DDRMat>
    get_void_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aVoidCellClusterIndex) const
    {
        MORIS_ERROR(false, "get_void_cell_local_coords_on_side_wrt_interp_cell(): not implemented for this cluster type");
        return moris::Matrix<moris::DDRMat>(0,0);
    }

    //##############################################
    // Size Access
    //##############################################
    virtual
    moris::uint
    get_num_primary_cells() const
    {
        MORIS_ERROR(false, "get_num_primary_cells(): not implemented for this cluster type");
        return 0;
    }

    virtual
    moris::uint
    get_num_void_cells() const
    {
        MORIS_ERROR(false, "get_num_void_cells(): not implemented for this cluster type");
        return 0;
    }

    virtual moris::uint get_num_vertices_in_cluster() const
    {
        MORIS_ERROR(false, "get_num_vertices_in_cluster(): not implemented for this cluster type");
        return 0;
    }


};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_CLUSTER_HPP_ */
