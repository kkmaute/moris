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
    virtual
    bool
    is_trivial() const  = 0;

    virtual
    moris::mtk::Cell const &
    get_interpolation_cell() const = 0;

    virtual
    moris::Cell<mtk::Cell const *> const &
    get_cells_in_side_cluster() const = 0;

    virtual
    moris::Matrix<moris::IndexMat>
    get_cell_side_ordinals() const  = 0;

    virtual
    moris::Cell<moris::mtk::Vertex const *> const &
    get_vertices_in_cluster() const = 0;

    //##############################################
    // Local Coordinate Access
    // (Pure Virtual)
    //##############################################
    virtual
    moris::Matrix<moris::DDRMat> const &
    get_vertices_local_coordinates_wrt_interp_cell() const = 0;

//    virtual
//    moris::Matrix<moris::DDRMat>
//    get_vertex_local_coordinates_wrt_interp_cell(moris::mtk::Vertex const * aVertex) const = 0;


    // ---------------------------------------------
    // EVERYTHING BELOW THIS LINE HAS A DEFAULT
    // IMPLEMENTATION
    // ---------------------------------------------

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
