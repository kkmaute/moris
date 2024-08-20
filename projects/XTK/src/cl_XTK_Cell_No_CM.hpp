/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_No_CM.hpp
 *
 */
#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_NO_CM_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_CELL_NO_CM_HPP_

#include "cl_MTK_Cell.hpp"
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"         //MRS/CNT/src
#include "cl_MTK_Vertex.hpp"     //MTK/src
#include "cl_MTK_Enums.hpp"      //MTK/src
#include "cl_MTK_Cell_Info.hpp"

using namespace moris;

//------------------------------------------------------------------------------
namespace moris::xtk
{
    //------------------------------------------------------------------------------
    /**
     * \brief This is the XTK cell implementation when there is a child mesh to use.
     */

    class Cell_XTK_No_CM : public moris::mtk::Cell
    {
        //------------------------------------------------------------------------------

      private:
        // mtk::Cell_Info const * mCellInfo = nullptr;
        Vector< mtk::Vertex* > mCellVertices;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------

        /**
         * trivial constructor
         */
        Cell_XTK_No_CM(){};

        //------------------------------------------------------------------------------

        Cell_XTK_No_CM(
                moris::moris_id                   aElementId,
                moris::moris_index                aElementIndex,
                moris::moris_index                aElementOwner,
                std::shared_ptr< mtk::Cell_Info > aCellInfo,
                Vector< mtk::Vertex* >            aVertices );

        //------------------------------------------------------------------------------

        /**
         * Destructor. Must be virtual.
         */
        ~Cell_XTK_No_CM(){};

        //------------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this cell
         */
        uint
        get_number_of_vertices() const
        {
            return mCellVertices.size();
        }

        //------------------------------------------------------------------------------

        /**
         * fills a Vector with pointers to connected vertices
         */
        // FIXME: SDF's Facet_Vertex causes this to not be able to return a reference.
        Vector< mtk::Vertex* >
        get_vertex_pointers() const
        {
            return mCellVertices;
        }

        //------------------------------------------------------------------------------

        void
        set_vertex_pointers( Vector< mtk::Vertex* >& aVertexPointers )
        {
            mCellVertices = aVertexPointers;
        }

        //------------------------------------------------------------------------------

        void
        replace_vertex_pointer( mtk::Vertex* aVertex, moris_index aIndex )
        {
            mCellVertices( aIndex ) = aVertex;
        }

        //------------------------------------------------------------------------------

        void
        remove_vertex_pointer( moris_index aIndex )
        {
            mCellVertices.erase( aIndex );
        }

        //------------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        Matrix< DDRMat >
        get_vertex_coords() const;

        //------------------------------------------------------------------------------

    };    // class Cell_XTK_No_CM

    //------------------------------------------------------------------------------

}    // namespace moris::xtk

//------------------------------------------------------------------------------

#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_ */
