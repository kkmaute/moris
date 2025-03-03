/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_ISC.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_CELL_ISC_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CELL_ISC_HPP_

#include "cl_MTK_Cell.hpp"

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::mtk
{
    class Cell_ISC : public mtk::Cell
    {
        //------------------------------------------------------------------------------

      private:
        Vector< mtk::Vertex* > mCellVertices;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        /*
         * trivial constructor
         */
        Cell_ISC(){};

        //------------------------------------------------------------------------------

        /*
         * constructor
         */
        Cell_ISC( moris::moris_id                 aCellId,
                moris::moris_index                aCellIndex,
                moris::moris_id                   aCellOwner,
                std::shared_ptr< mtk::Cell_Info > aCellInfo,
                const Vector< mtk::Vertex* >&     aVertices );

        //------------------------------------------------------------------------------

        /**
         * Destructor
         */
        ~Cell_ISC() override{};

        //------------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this cell
         */
        uint
        get_number_of_vertices() const override
        {
            return mCellVertices.size();
        }

        //------------------------------------------------------------------------------

        /**
         * fills a Vector with pointers to connected vertices
         */
        // FIXME: SDF's Facet_Vertex causes this to not be able to return a reference.
        Vector< mtk::Vertex* >
        get_vertex_pointers() const override
        {
            return mCellVertices;
        }

        //------------------------------------------------------------------------------

        // TODO MESHCLEANUP
        void
        remove_vertex_pointer( moris_index aIndex ) override
        {
            std::cout << "In MTK Cell ISC" << '\n';
        }

        //------------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        Matrix< DDRMat >
        get_vertex_coords() const override;

        //------------------------------------------------------------------------------

    };    // class Cell_ISC

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_CELL_ISC_HPP_ */
