/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_CM.cpp
 *
 */

#include "cl_XTK_Cell_CM.hpp"
#include "cl_XTK_Background_Mesh.hpp"

namespace xtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------
    Cell_XTK_CM::Cell_XTK_CM(){}

    Cell_XTK_CM::Cell_XTK_CM(moris::moris_id aElementId,
                             moris::moris_index aElementIndex,
                             moris::moris_index aElementOwner,
                             moris::moris_index aCMElementIndex,
                             xtk::Child_Mesh *aChildMeshPtr,
                             xtk::Background_Mesh *aBackgroundMeshPtr)
                             : Cell(aElementId,aElementIndex,aElementOwner, aChildMeshPtr->get_cell_info_sp()),
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

    Cell_XTK_CM::~Cell_XTK_CM(){}

    uint
    Cell_XTK_CM::get_number_of_vertices() const
    {
        return mChildMeshPtr->get_element_to_node().n_cols();
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
        tTotal +=  sizeof(mCMElementIndex);
        tTotal +=  sizeof(mChildMeshPtr);
        tTotal +=  sizeof(mBackgroundMeshPtr);
        return tTotal;
    }

} // namespace xtk

