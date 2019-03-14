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
                MORIS_ERROR(0, "get_owner not implemented in XTK Cell");
                return 0;
            }

//------------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            //FIXME: SDF's Triangle_Vertex causes this to not be able to return a reference.
            moris::Cell< Vertex* >
            get_vertex_pointers() const
            {
                MORIS_ERROR(0, "get_vertex_pointers not implemented in XTK Cell");
                return moris::Cell< Vertex* >(0);
            }

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
                return Geometry_Type::TET;
            }

//------------------------------------------------------------------------------

            /**
             * returns the order of the element
             */
            virtual Interpolation_Order
            get_interpolation_order() const
            {
                return Interpolation_Order::LINEAR;
            }

//------------------------------------------------------------------------------

        private:
            moris::moris_id    mElementId;
            moris::moris_index mElementIndex;
            moris::moris_index mCMElementIndex;       /* Needed to access connectivity (verts) */
            xtk::Child_Mesh*      mChildMeshPtr;      /* Needed to access connectivity (verts) */
            xtk::Background_Mesh* mBackgroundMeshPtr; /* Needed to access coordinates */

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------



#endif /* PROJECTS_XTK_SRC_XTK_CL_MTK_CELL_XTK_IMPL_HPP_ */
