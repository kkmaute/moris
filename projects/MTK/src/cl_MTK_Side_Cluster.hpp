/*
 * cl_MTK_Side_Cluster.hpp
 *
 *  Created on: May 9, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_HPP_

#include "typedefs.hpp"
#include "fn_assert.hpp"
#include "cl_Matrix.hpp"
#include "cl_Cell.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Cell.hpp"


namespace moris
{
namespace mtk
{

class Side_Cluster
{
public:
    /*
     * Default constructor
     */
    Side_Cluster(){}

    //##############################################
    // Cell Side Ordinals/Vertex Access
    // (Pure Virtual)
    //##############################################
    /*!
     * Indicates there is a 1 to 1 relationship between
     * integration cell and interpolation cells in this cluster
     */
    virtual
    bool
    is_trivial() const  = 0;

    /*!
     * Get interpolation cell interpolating into this side cluster
     */
    virtual
    moris::mtk::Cell const &
    get_interpolation_cell() const = 0;

    /*!
     * Get all integration cells in this side cluster
     */
    virtual
    moris::Cell<mtk::Cell const *> const &
    get_cells_in_side_cluster() const = 0;

    /*!
     * Return all integration cell side ordinals in cluster
     */
    virtual
    moris::Matrix<moris::IndexMat>
    get_cell_side_ordinals() const  = 0;

    /*!
     * Single side ordinal version of above
     */
    virtual
    moris_index
    get_cell_side_ordinal(moris::moris_index aCellIndexInCluster) const = 0;


    /*!
     * Returns all the vertices in this cluster
     */
    virtual
    moris::Cell<moris::mtk::Vertex const *> const &
    get_vertices_in_cluster() const = 0;

    //##############################################
    // Local Coordinate Access
    // (Pure Virtual)
    //##############################################

    /*
     * Access the full array of local coordinates
     */
    virtual
    moris::Matrix<moris::DDRMat> const &
    get_vertices_local_coordinates_wrt_interp_cell() const = 0;

    /*
     * Access a single local coordinate of a vertex
     */
    virtual
    moris::Matrix<moris::DDRMat>
    get_vertex_local_coordinate_wrt_interp_cell( moris::mtk::Vertex const * aVertex ) const  = 0;

    //##############################################
    // Size Access
    // (Pure Virtual)
    //##############################################
    /*!
     * Size of the xsi vector in this side cluster
     */
    virtual
    moris_index
    get_dim_of_param_coord() const  = 0;

    //----------------------------------------------------------------
    // EVERYTHING BELOW THIS LINE HAS A DEFAULT IMPLEMENTATION
    //----------------------------------------------------------------

    virtual
    moris::Matrix<moris::IndexMat>
    get_cell_indices_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumCells = this->get_num_sides_in_cluster();

        // cell access
        moris::Cell<moris::mtk::Cell const *> const & tCells = this->get_cells_in_side_cluster();

        // initialize output
        moris::Matrix<moris::IndexMat> tCellIndices(1,tNumCells);

        // get cell indices and store
        for(moris::uint i = 0 ; i < tNumCells; i++)
        {
            tCellIndices(i) = tCells(i)->get_index();
        }

        return tCellIndices;
    }

    // ---------------------------------------------

    virtual
    moris::moris_index
    get_interpolation_cell_index() const
    {
        return get_interpolation_cell().get_index();
    }

    virtual
    moris::Matrix<moris::IndexMat>
    get_vertex_indices_in_cluster() const
    {
        // number of cells in cluster
         moris::uint tNumVertices = this->get_num_vertices_in_cluster();

         // cell access
         moris::Cell<moris::mtk::Vertex const *> const & tVertices = this->get_vertices_in_cluster();

         // initialize output
         moris::Matrix<moris::IndexMat> tVertexIndices(1,tNumVertices);

         // get cell indices and store
         for(moris::uint i = 0 ; i < tNumVertices; i++)
         {
             tVertexIndices(i) = tVertices(i)->get_index();
         }

         return tVertexIndices;
    }

    // ---------------------------------------------



    //##############################################
    // Cell/Vertex Id Access
    //##############################################
    virtual
    moris::Matrix<moris::IdMat>
    get_cell_ids_in_cluster() const
    {
        // number of cells in cluster
        moris::uint tNumCells = this->get_num_sides_in_cluster();

        // cell access
        moris::Cell<moris::mtk::Cell const *> const & tCells = this->get_cells_in_side_cluster();

        // initialize output
        moris::Matrix<moris::IdMat> tCellIds(1,tNumCells);

        // get cell indices and store
        for(moris::uint i = 0 ; i < tNumCells; i++)
        {
            tCellIds(i) = tCells(i)->get_id();
        }

        return tCellIds;
    }

    virtual
    moris::Matrix<moris::IndexMat>
    get_vertex_ids_in_cluster() const
    {
        // number of cells in cluster
         moris::uint tNumVertices = this->get_num_vertices_in_cluster();

         // cell access
         moris::Cell<moris::mtk::Vertex const *> const & tVertices = this->get_vertices_in_cluster();

         // initialize output
         moris::Matrix<moris::IndexMat> tVertexIds(1,tNumVertices);

         // get cell indices and store
         for(moris::uint i = 0 ; i < tNumVertices; i++)
         {
             tVertexIds(i) = tVertices(i)->get_id();
         }

         return tVertexIds;
    }

    //##############################################
    // Local Coordinate access
    //##############################################

    /*!
     * Access an integration cells parametric coordinates on a side
     * @param[in] - Local integration cell index with respect to the cluster (not proc local index)
     */
    virtual
    moris::Matrix<moris::DDRMat>
    get_cell_local_coords_on_side_wrt_interp_cell(moris::moris_index aClusterLocalIndex) const
    {
        MORIS_ASSERT(aClusterLocalIndex < (moris_index)this->get_num_sides_in_cluster(),"Integration Cell Cluster index out of bounds");

        // get side ordinal of interest
        moris_index tSideOrdinal = this->get_cell_side_ordinal(aClusterLocalIndex);

        // get the integration cell of interest
        moris::mtk::Cell const * tIntegrationCell = this->get_cells_in_side_cluster()(aClusterLocalIndex);

        // get the vertex pointers on the side
        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = tIntegrationCell->get_vertices_on_side_ordinal(tSideOrdinal);

        // allocate output (nnode x dim_xsi)
        moris::Matrix<moris::DDRMat> tVertexParamCoords( tVerticesOnSide.size(), this->get_dim_of_param_coord());

        // iterate through vertices and collect local coordinates
        for(moris::uint i = 0; i < tVerticesOnSide.size(); i++)
        {
            tVertexParamCoords.get_row(i) = this->get_vertex_local_coordinate_wrt_interp_cell(tVerticesOnSide(i)).get_row(0);
        }

        return tVertexParamCoords;
    }

    //##############################################
    // Size Access
    //##############################################
    virtual
    moris::uint
    get_num_sides_in_cluster() const
    {
        return this->get_cells_in_side_cluster().size();
    }

    virtual
    moris::uint
    get_num_vertices_in_cluster() const
    {
        return this->get_vertices_in_cluster().size();
    }

};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_HPP_ */
