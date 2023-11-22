/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_CM.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_

#include "cl_MTK_Cell.hpp"
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Cell.hpp"     //MRS/CNT/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK/src
#include "cl_MTK_Enums.hpp"     //MTK/src
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"

namespace xtk
{
    class Background_Mesh;
}

//------------------------------------------------------------------------------
namespace xtk
{
    //------------------------------------------------------------------------------

    /**
     * \brief This is the XTK cell implementation when there is a child mesh to use.
     */
    class Cell_XTK_CM : public mtk::Cell
    {

      private:
        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Cell_XTK_CM();

        Cell_XTK_CM(
                moris::moris_id       aElementId,
                moris::moris_index    aElementIndex,
                moris::moris_index    aElementOwner,
                moris::moris_index    aCMElementIndex,
                xtk::Child_Mesh      *aChildMeshPtr,
                xtk::Background_Mesh *aBackgroundMeshPtr );
        //------------------------------------------------------------------------------

        /**
         * Destructor
         */
        ~Cell_XTK_CM();

        //------------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this cell
         */
        uint
        get_number_of_vertices() const;

        //------------------------------------------------------------------------------

        /**
         * fills a moris::cell with pointers to connected vertices
         */
        // FIXME: SDF's Facet_Vertex causes this to not be able to return a reference.
        moris::Cell< mtk::Vertex * >
        get_vertex_pointers() const;

        //------------------------------------------------------------------------------

        void
        remove_vertex_pointer( moris_index aIndex );

        //------------------------------------------------------------------------------

        /**
         * returns a Mat with IDs of connected vertices
         */
        Matrix< IdMat >
        get_vertex_ids() const;

        //------------------------------------------------------------------------------

        /**
         * returns a Mat with indices of connected vertices
         */
        Matrix< IndexMat >
        get_vertex_inds() const;

        //------------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        Matrix< DDRMat >
        get_vertex_coords() const;

        //------------------------------------------------------------------------------
        /*!
         * @brief capacity of the cell
         */
        size_t
        capacity();
        //------------------------------------------------------------------------------

      private:
        moris::moris_index mCMElementIndex; /* Needed to access connectivity (verts) */

        xtk::Child_Mesh      *mChildMeshPtr;      /* Needed to access connectivity (verts) */
        xtk::Background_Mesh *mBackgroundMeshPtr; /* Needed to access coordinates */

        //------------------------------------------------------------------------------

    }; // class Cell_XTK_CM

    //------------------------------------------------------------------------------

}    // namespace xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_ */
