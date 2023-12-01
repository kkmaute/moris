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

#include "typedefs.hpp"
#include "cl_Cell.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        class Cell_ISC : public mtk::Cell
        {
            //------------------------------------------------------------------------------

          private:
            moris::Cell< mtk::Vertex* > mCellVertices;

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
                    moris::Cell< mtk::Vertex* >       aVertices );

            //------------------------------------------------------------------------------

            /**
             * Destructor
             */
            ~Cell_ISC(){};

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
             * fills a moris::cell with pointers to connected vertices
             */
            // FIXME: SDF's Facet_Vertex causes this to not be able to return a reference.
            moris::Cell< mtk::Vertex* >
            get_vertex_pointers() const
            {
                return mCellVertices;
            }

            //------------------------------------------------------------------------------

            // TODO MESHCLEANUP
            void
            remove_vertex_pointer( moris_index aIndex )
            {
                std::cout << "In MTK Cell ISC" << std::endl;
            }

            //------------------------------------------------------------------------------

            /**
             * returns a Mat of dimension
             * < number of vertices * number of dimensions >
             */
            Matrix< DDRMat >
            get_vertex_coords() const;

            //------------------------------------------------------------------------------

        };    // class Cell_ISC

        //------------------------------------------------------------------------------

    } /* end namespace mtk */
} /* end namespace moris */

#endif /* PROJECTS_MTK_SRC_CL_MTK_CELL_ISC_HPP_ */
