/*
 * cl_MTK_Cell.hpp
 *
 *  Created on: Jul 23, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_CELL_HPP_
#define SRC_MESH_CL_MTK_CELL_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_Matrix.hpp"
#include "fn_isrow.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_MTK_Cell_Info_Hex8.hpp"

//------------------------------------------------------------------------------
namespace moris
{
namespace mtk
{
//------------------------------------------------------------------------------
/**
 * \brief the mtk::Cell class provides the cell information that is
 * provided by the mesh.
 */

class Cell
{
    //------------------------------------------------------------------------------
public:
    //------------------------------------------------------------------------------

    /**
     * trivial constructor
     */
    Cell(){};

    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    virtual
    ~Cell(){};

    //------------------------------------------------------------------------------

    /**
     * returns the domain wide id of the cell
     *
     * @return moris_id ID
     */
    virtual moris_id
    get_id() const = 0;

    //------------------------------------------------------------------------------

    /**
     * returns the local index of the cell
     *
     * @return moris_index ID
     */
    virtual moris_index
    get_index() const = 0;


    //------------------------------------------------------------------------------

    /**
     * returns the proc id of the owner of this cell
     * ( this information is needed for STK )
     */
    virtual moris_id
    get_owner() const = 0;

    //------------------------------------------------------------------------------

    /*!
     * Returns the level that this cell is on. For most meshes this returns 0. However,
     * for HMR this is not trivial
     */
    virtual uint
    get_level() const;

    //------------------------------------------------------------------------------

    /**
     * fills a moris::cell with pointers to connected vertices
     */
    virtual moris::Cell< Vertex* >
    get_vertex_pointers() const = 0;

    //------------------------------------------------------------------------------

    /**
     * tells how many vertices are connected to this cell
     */
    virtual uint
    get_number_of_vertices() const;

    //------------------------------------------------------------------------------

    /**
     * returns a Mat with IDs of connected vertices
     */
    virtual Matrix< IdMat >
    get_vertex_ids() const;

    //------------------------------------------------------------------------------

    /**
     * returns a Mat with indices of connected vertices
     */
    virtual Matrix< IndexMat >
    get_vertex_inds() const;

    //------------------------------------------------------------------------------

    /**
     * returns a Mat of dimension
     * < number of vertices * number of dimensions >
     */
    virtual Matrix< DDRMat >
    get_vertex_coords() const = 0;

    //------------------------------------------------------------------------------

    virtual
    moris::Cell<mtk::Vertex_Interpolation*>
    get_vertex_interpolations( const uint aOrder ) const;

    //------------------------------------------------------------------------------

    /*!
     * get vertices on side ordinal.
     * This functions is needed for side clustering
     */
    virtual
    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const;

    //------------------------------------------------------------------------------

    /*!
     * get vertices on side ordinal that define the geometry (i.e. the corner nodes)
     * This functions is needed for side clustering
     */
    virtual
    moris::Cell<moris::mtk::Vertex const *>
    get_geometric_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const;

    //------------------------------------------------------------------------------

    /*!
     * Get vertex coordinates on side ordinal
     */

    virtual
    moris::Matrix<moris::DDRMat>
    get_cell_physical_coords_on_side_ordinal(moris::moris_index aSideOrdinal) const;

    //------------------------------------------------------------------------------

    /*!
     * get vertices on side ordinal.
     * This functions is needed for side clustering
     */
    moris::Matrix< IndexMat >
    get_vertices_ind_on_side_ordinal(moris::moris_index aSideOrdinal) const;

    //------------------------------------------------------------------------------

    /**
     * returns an enum that defines the geometry type of the element
     */
    virtual Geometry_Type
    get_geometry_type() const = 0;

    //------------------------------------------------------------------------------

    /*!
     * Compute facet normal
     */
    virtual
    moris::Matrix<moris::DDRMat>
    compute_outward_side_normal(moris::moris_index aSideOrdinal) const;

    //------------------------------------------------------------------------------

    /*
     * Volume in 3D, Surface Area in 2D
     */
    virtual
    moris::real
    compute_cell_measure() const;

    //------------------------------------------------------------------------------

    /*
     * Surface Area on side of cell in 3D, line length on side in 2D
     */
    virtual
    moris::real
    compute_cell_side_measure(moris_index const & aCellSideOrd) const;

    //------------------------------------------------------------------------------

    /**
     * returns the order of the element
     */
    virtual Interpolation_Order
    get_interpolation_order() const = 0;

    //------------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------

#endif /* SRC_MESH_CL_MTK_CELL_HPP_ */
