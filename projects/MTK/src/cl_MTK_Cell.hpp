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

    /**
     * fills a moris::cell with pointers to connected vertices
     */
    //FIXME: SDF's Triangle_Vertex causes this to not be able to return a reference.
    virtual moris::Cell< Vertex* >
    get_vertex_pointers() const = 0;

    /**
     * tells how many vertices are connected to this cell
     */
    virtual uint
    get_number_of_vertices() const
    {
        return this->get_vertex_pointers().size();
    };
    //------------------------------------------------------------------------------

    /**
     * returns a Mat with IDs of connected vertices
     */
    virtual Matrix< IdMat >
    get_vertex_ids() const
    {
        uint tNumVertices = this->get_number_of_vertices();
        moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

        Matrix< IdMat > tVertexIds(1, tNumVertices);
        for(uint i = 0; i<tNumVertices; i++)
        {
            tVertexIds(i) = tVertices(i)->get_id();
        }
        return tVertexIds;
    }

    //------------------------------------------------------------------------------

    /**
     * returns a Mat with indices of connected vertices
     */
    virtual Matrix< IndexMat >
    get_vertex_inds() const
    {
        uint tNumVertices = this->get_number_of_vertices();
        moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

        Matrix< IdMat > tVertexInds(1, tNumVertices);
        for(uint i = 0; i<tNumVertices; i++)
        {
            tVertexInds(i) = tVertices(i)->get_index();
        }
        return tVertexInds;
    }

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
    get_vertex_interpolations( const uint aOrder ) const
    {
        uint tNumVerts = this->get_number_of_vertices();
        moris::Cell< mtk::Vertex* > tVertexPointers = this->get_vertex_pointers();
        moris::Cell<mtk::Vertex_Interpolation*> tVertexInterp(tNumVerts);

        for(moris::uint i = 0; i < tNumVerts; i++)
        {
            tVertexInterp(i) =tVertexPointers(i)->get_interpolation(aOrder);
        }

        return tVertexInterp;

    }
    //------------------------------------------------------------------------------

    /*!
     * get vertices on side ordinal.
     * This functions is needed for side clustering
     */
    virtual
    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {
        MORIS_ERROR(0,"get_vertices_on_side_ordinal has no default implementation");
        return  moris::Cell<moris::mtk::Vertex const *>(0);
    }

    /*!
     * Get vertex coordinates on side ordinal
     */

    virtual
    moris::Matrix<moris::DDRMat>
    get_cell_physical_coords_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {

        // FIXME: Add assert to check side ordinal

        // get the vertex pointers on the side
        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide = this->get_vertices_on_side_ordinal(aSideOrdinal);

        // allocate output coords (note we do not know the spatial dimension at this time)
        moris::Matrix<moris::DDRMat> tVertexPhysCoords(0,0);

        // iterate through vertices and collect local coordinates
        for(moris::uint i = 0; i < tVerticesOnSide.size(); i++)
        {
            moris::Matrix<moris::DDRMat> tVertexCoord = tVerticesOnSide(i)->get_coords();

            if( i == 0 )
            {
                MORIS_ASSERT(isrow(tVertexCoord),"Default implementation assumes row based coordinates");
                tVertexPhysCoords.resize(tVerticesOnSide.size(), tVertexCoord.numel());
            }

            tVertexPhysCoords.get_row(i) = tVertexCoord.get_row(0);
        }

        return tVertexPhysCoords;
    }

    //------------------------------------------------------------------------------

    /*!
     * get vertices on side ordinal.
     * This functions is needed for side clustering
     */
    moris::Matrix< IndexMat >
    get_vertices_ind_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {
        moris::Cell<moris::mtk::Vertex const *> tVertices = this->get_vertices_on_side_ordinal(aSideOrdinal);

        uint tNumVertices = tVertices.size();

        Matrix< IndexMat > tVertexInd( 1, tNumVertices );

        for(uint i = 0; i < tNumVertices; i++ )
        {
            tVertexInd( 0, i ) = tVertices( i )->get_index();
        }
        return  tVertexInd;
    }

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
    compute_outward_side_normal(moris::moris_index aSideOrdinal) const
    {
        MORIS_ERROR(0,"compute_outward_side_normal has no default implementation");
        return  moris::Matrix<moris::DDRMat>(0,0);
    }

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
