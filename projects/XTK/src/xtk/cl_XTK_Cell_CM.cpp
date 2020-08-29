/*
 * cl_MTK_Cell_XTK_CM_Impl.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: doble
 */

#include "cl_XTK_Cell_CM.hpp"
#include "cl_XTK_Background_Mesh.hpp"

namespace xtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------
    Cell_XTK_CM::Cell_XTK_CM(){};

    Cell_XTK_CM::Cell_XTK_CM(moris::moris_id aElementId,
                             moris::moris_index aElementIndex,
                             moris::moris_index aElementOwner,
                             moris::moris_index aCMElementIndex,
                             xtk::Child_Mesh *aChildMeshPtr,
                             xtk::Background_Mesh *aBackgroundMeshPtr) 
                             : mElementId(aElementId),
                               mElementIndex(aElementIndex),
                               mElementOwner(aElementOwner),
                               mCMElementIndex(aCMElementIndex),
                               mChildMeshPtr(aChildMeshPtr),
                               mBackgroundMeshPtr(aBackgroundMeshPtr)
    {
    }

    moris::Cell<mtk::Vertex *>
    Cell_XTK_CM::get_vertex_pointers() const
    {
        Matrix<IndexMat> tVertexIndices = this->get_vertex_inds();
        moris::Cell<mtk::Vertex *> tVertices(tVertexIndices.numel());

        for (moris::uint i = 0; i < tVertices.size(); i++)
        {
            tVertices(i) = &mBackgroundMeshPtr->get_mtk_vertex(tVertexIndices(i));
        }
        return tVertices;
    }

    // ----------------------------------------------------------------------------------

    Cell_XTK_CM::~Cell_XTK_CM(){};

    // ----------------------------------------------------------------------------------
    moris_id
    Cell_XTK_CM::get_id() const
    {
        return mElementId;
    }
    // ----------------------------------------------------------------------------------
    moris_index
    Cell_XTK_CM::get_index() const
    {
        return mElementIndex;
    }
    // ----------------------------------------------------------------------------------
    uint
    Cell_XTK_CM::get_number_of_vertices() const
    {
        return mChildMeshPtr->get_element_to_node().n_cols();
    }

    // ----------------------------------------------------------------------------------
    moris_id
    Cell_XTK_CM::get_owner() const
    {
        return mElementOwner;
    }

    // ----------------------------------------------------------------------------------

    Matrix<IdMat>
    Cell_XTK_CM::get_vertex_ids() const
    {
        return mChildMeshPtr->get_element_to_node_glob_ids(mCMElementIndex);
    }

    // ----------------------------------------------------------------------------------
    Matrix<IndexMat>
    Cell_XTK_CM::get_vertex_inds() const
    {
        return mChildMeshPtr->get_element_to_node().get_row(mCMElementIndex);
    }

    // ----------------------------------------------------------------------------------
    mtk::Geometry_Type
    Cell_XTK_CM::get_geometry_type() const
    {
        return mChildMeshPtr->get_child_geometry_type();
    }
    // ----------------------------------------------------------------------------------
    mtk::Interpolation_Order
    Cell_XTK_CM::get_interpolation_order() const
    {
        return mChildMeshPtr->get_child_interpolation_order();
    }
    // ----------------------------------------------------------------------------------
    moris::Cell<moris::mtk::Vertex const *>
    Cell_XTK_CM::get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {

        moris::Cell<mtk::Vertex *> tVertices = this->get_vertex_pointers();

        moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = mChildMeshPtr->get_cell_info()->get_node_to_facet_map(aSideOrdinal);

        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(tNodeOrdsOnSide.numel());
        for (moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++)
        {
            tVerticesOnSide(i) = tVertices(tNodeOrdsOnSide(i));
        }

        return tVerticesOnSide;
    }
    // ----------------------------------------------------------------------------------
    moris::Matrix<moris::DDRMat>
    Cell_XTK_CM::compute_outward_side_normal(moris::moris_index aSideOrdinal) const
    {
        // get the vertex coordinates
        moris::Matrix<moris::DDRMat> tVertexCoords = this->get_vertex_coords();

        // Get vector along these edges
        moris::Matrix<moris::DDRMat> tEdge0Vector(tVertexCoords.numel(), 1);
        moris::Matrix<moris::DDRMat> tEdge1Vector(tVertexCoords.numel(), 1);

        // Get the nodes which need to be used to compute normal
        moris::Matrix<moris::IndexMat> tEdgeNodesForNormal = mChildMeshPtr->get_cell_info()->get_node_map_outward_normal(aSideOrdinal);

        // Get vector along these edges
        tEdge0Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1, 0)) - tVertexCoords.get_row(tEdgeNodesForNormal(0, 0)));
        tEdge1Vector = moris::linalg_internal::trans(tVertexCoords.get_row(tEdgeNodesForNormal(1, 1)) - tVertexCoords.get_row(tEdgeNodesForNormal(0, 1)));

        // Take the cross product to get the normal
        Matrix<DDRMat> tOutwardNormal = moris::cross(tEdge0Vector, tEdge1Vector);

        // Normalize
        Matrix<DDRMat> tUnitOutwardNormal = tOutwardNormal / moris::norm(tOutwardNormal);

        return tUnitOutwardNormal;
    }
    // ----------------------------------------------------------------------------------
    moris::real
    Cell_XTK_CM::compute_cell_measure() const
    {
        return mChildMeshPtr->get_cell_info()->compute_cell_size(this);
    }
    // ----------------------------------------------------------------------------------
    moris::real
    Cell_XTK_CM::compute_cell_side_measure(moris_index const &aSideOrdinal) const
    {
        return mChildMeshPtr->get_cell_info()->compute_cell_side_size(this, aSideOrdinal);
    }

    // ----------------------------------------------------------------------------------
    // Cell get functions
    // ----------------------------------------------------------------------------------
    Matrix<DDRMat>
    Cell_XTK_CM::get_vertex_coords() const
    {
        return mBackgroundMeshPtr->get_selected_node_coordinates_loc_inds(this->get_vertex_inds());
    }

    // ----------------------------------------------------------------------------------

    size_t
    Cell_XTK_CM::capacity()
    {
        size_t tTotal = 0;
        tTotal +=  sizeof(mElementId);
        tTotal +=  sizeof(mElementIndex);
        tTotal +=  sizeof(mElementOwner);
        tTotal +=  sizeof(mCMElementIndex);
        tTotal +=  sizeof(mChildMeshPtr);
        tTotal +=  sizeof(mBackgroundMeshPtr);
        return tTotal;
    }

} // namespace xtk
