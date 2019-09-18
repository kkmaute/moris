/*
 * cl_MTK_Cell_XTK_Impl.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_

#include "cl_MTK_Cell.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Vertex.hpp" //MTK/src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_MTK_Tet4_Cell_Info.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"

namespace xtk
{
class Background_Mesh;
}

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

class Cell_XTK: public Cell
{
private:

    //------------------------------------------------------------------------------
public:
    //------------------------------------------------------------------------------

    /**
     * trivial constructor
     */
    Cell_XTK(){};

    Cell_XTK(moris::moris_id       aElementId,
             moris::moris_index    aElementIndex,
             moris::moris_index    aElementOwner,
             moris::moris_index    aCMElementIndex,
             xtk::Child_Mesh*      aChildMeshPtr,
             xtk::Background_Mesh* aBackgroundMeshPtr);
    //------------------------------------------------------------------------------

    /**
     * Destructor. Must be virtual.
     */
    ~Cell_XTK(){};

    //------------------------------------------------------------------------------

    /**
     * returns the domain wide id of the cell
     *
     * @return moris_id ID
     */
    moris_id
    get_id() const
    {
        return mElementId;
    }

    //------------------------------------------------------------------------------

    /**
     * returns the local index of the cell
     *
     * @return moris_index ID
     */
    moris_index
    get_index() const
    {
        return mElementIndex;
    }

    //------------------------------------------------------------------------------

    /**
     * tells how many vertices are connected to this cell
     */
    uint
    get_number_of_vertices() const
    {
        return mChildMeshPtr->get_element_to_node().n_cols();
    }

    //------------------------------------------------------------------------------

    /**
     * returns the proc id of the owner of this cell
     * ( this information is needed for STK )
     */
    moris_id
    get_owner() const
    {
        return mElementOwner;
    }

    //------------------------------------------------------------------------------

    /**
     * fills a moris::cell with pointers to connected vertices
     */
    //FIXME: SDF's Triangle_Vertex causes this to not be able to return a reference.
    moris::Cell< Vertex* >
    get_vertex_pointers() const;

    //------------------------------------------------------------------------------

    /**
     * returns a Mat with IDs of connected vertices
     */
    Matrix< IdMat >
    get_vertex_ids() const
    {
        return mChildMeshPtr->get_element_to_node_glob_ids(mCMElementIndex);
    }

    //------------------------------------------------------------------------------

    /**
     * returns a Mat with indices of connected vertices
     */
    Matrix< IndexMat >
    get_vertex_inds() const
    {
        return mChildMeshPtr->get_element_to_node().get_row(mCMElementIndex);
    }

    //------------------------------------------------------------------------------

    /**
     * returns a Mat of dimension
     * < number of vertices * number of dimensions >
     */
    Matrix< DDRMat >
    get_vertex_coords() const;

    //------------------------------------------------------------------------------

    /**
     * returns an enum that defines the geometry type of the element
     */
    Geometry_Type
    get_geometry_type() const
    {
        return mChildMeshPtr->get_child_geometry_type();
    }

    //------------------------------------------------------------------------------

    /**
     * returns the order of the element
     */
    virtual Interpolation_Order
    get_interpolation_order() const
    {
        return mChildMeshPtr->get_child_interpolation_order();
    }

    //------------------------------------------------------------------------------

    moris::Cell<moris::mtk::Vertex const *>
    get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {

        moris::Cell< Vertex* > tVertices = this->get_vertex_pointers();

        moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = mChildMeshPtr->get_cell_info()->get_node_to_facet_map(aSideOrdinal);

        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(tNodeOrdsOnSide.numel());
        for(moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++)
        {
            tVerticesOnSide(i) = tVertices(tNodeOrdsOnSide(i));
        }

        return tVerticesOnSide;
    }

    //------------------------------------------------------------------------------

    moris::Matrix<moris::DDRMat>
    compute_outward_side_normal(moris::moris_index aSideOrdinal) const
    {
        // get the vertex coordinates
        moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

        // Get vector along these edges
        moris::Matrix<moris::DDRMat> tEdge0Vector(tVertexCoords.numel(),1);
        moris::Matrix<moris::DDRMat> tEdge1Vector(tVertexCoords.numel(),1);

        // Get the nodes which need to be used to compute normal
        moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = mChildMeshPtr->get_cell_info()->get_node_map_outward_normal(aSideOrdinal);

        // Get vector along these edges
        tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,0)));
        tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1,1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0,1)));

        // Take the cross product to get the normal
        Matrix<DDRMat> tOutwardNormal = moris::cross(tEdge0Vector,tEdge1Vector);

        // Normalize
        Matrix<DDRMat> tUnitOutwardNormal = tOutwardNormal / moris::norm(tOutwardNormal);


        return tUnitOutwardNormal;

    }

    //------------------------------------------------------------------------------

private:
    moris::moris_id       mElementId;
    moris::moris_index    mElementIndex;
    moris::moris_index    mElementOwner;
    moris::moris_index    mCMElementIndex;    /* Needed to access connectivity (verts) */

    xtk::Child_Mesh*      mChildMeshPtr;      /* Needed to access connectivity (verts) */
    xtk::Background_Mesh* mBackgroundMeshPtr; /* Needed to access coordinates */

    //------------------------------------------------------------------------------
};

//------------------------------------------------------------------------------
} /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------



#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_ */
