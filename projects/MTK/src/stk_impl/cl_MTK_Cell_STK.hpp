/*
 * cl_MTK_Cell_STK.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Cell.hpp" //MRS/CON/src
#include "cl_MTK_Vertex_STK.hpp" //MTK/src
#include "../cl_MTK_Cell.hpp" //MTK/src
#include "../cl_MTK_Enums.hpp" //MTK/src
#include "fn_print.hpp"


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

        class Cell_STK: public moris::mtk::Cell
        {
            enum CellType        mCellType;
            moris_id             mCellId;
            moris_index          mCellInd;
            moris::Cell<Vertex*> mCellVertices;
            Mesh_Temp*  mSTKMeshData;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Cell_STK(){};

//------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Cell_STK( enum CellType               aCellType,
                           moris_id               aCellId,
                           moris_index            aCellInd,
                           moris::Cell<Vertex*> const & aCellVertices,
                           Mesh_Temp* aStkImplementation):
                         mCellId(aCellId),
                         mCellInd(aCellInd),
                         mCellVertices(aCellVertices),
                         mSTKMeshData(aStkImplementation)
            {

            };
            //------------------------------------------------------------------------------

            /**
             * Destructor. Must be virtual.
             */
            ~Cell_STK(){};

            //------------------------------------------------------------------------------

            /**
             * returns the domain wide id of the cell
             *
             * @return luint ID
             */
            moris_id
            get_id() const
            {
                return mCellId;
            };

            //------------------------------------------------------------------------------

            /**
             * returns the domain wide id of the cell
             *
             * @return luint ID
             */
            moris_index
            get_index() const
            {
                return mCellInd;
            }

//------------------------------------------------------------------------------

            /**
             * tells how many vertices are connected to this cell
             */
            uint
            get_number_of_vertices() const
            {
              return mCellVertices.size();
            };

//------------------------------------------------------------------------------

            /**
             * returns the proc id of the owner of this cell
             * ( this information is needed for STK )
             */
            moris_id
            get_owner() const
            {
                return mSTKMeshData->get_entity_owner( mCellInd, EntityRank::ELEMENT);
            }

//------------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            moris::Cell< Vertex* >
            get_vertex_pointers() const
            {
                return mCellVertices;
            }


//------------------------------------------------------------------------------
            /**
             * returns a Matrix with IDs of connected vertices
             */
            Matrix< IdMat >
            get_vertex_ids() const
            {
                size_t tNumVertices = this->get_number_of_vertices();
                Matrix< IdMat > tVertexIds(1, tNumVertices);
                for(size_t i = 0; i<tNumVertices; i++)
                {
                    tVertexIds(i) = mCellVertices(i)->get_id();
                }
                return tVertexIds;
            }

            /**
             * returns a Mat with indices of connected vertices
             */
            Matrix< IndexMat >
            get_vertex_inds() const
            {
                size_t tNumVertices = this->get_number_of_vertices();
                Matrix< IndexMat > tVertexInds(1, tNumVertices);
                for(size_t i = 0; i<tNumVertices; i++)
                {
                    tVertexInds(i) = mCellVertices(i)->get_index();
                }
                return tVertexInds;
            }

//------------------------------------------------------------------------------

            /**
             * returns a Matrix of dimension
             * < number of vertices * number of dimensions >
             */
            Matrix< DDRMat >
            get_vertex_coords() const
            {
                size_t tNumVertices = this->get_number_of_vertices();
                Matrix< DDRMat > tVertexCoords(tNumVertices, mSTKMeshData->get_spatial_dim());
                for(size_t i = 0; i<tNumVertices; i++)
                {
                    tVertexCoords.get_row(i) = mCellVertices(i)->get_coords().get_row(0);
                }
                return tVertexCoords;
            }

//------------------------------------------------------------------------------

            /**
             * returns an enum that defines the geometry type of the element
             */
            Geometry_Type
            get_geometry_type() const
            {
                MORIS_ASSERT(false,"Not implemented");
                return Geometry_Type::UNDEFINED;
            }

//------------------------------------------------------------------------------

            /**
             * T-Matrix is calculated for vertices if this flag is set
             */
            void
            set_t_matrix_flag()
            {
                MORIS_ASSERT(false,"Not implemented");
            }

//------------------------------------------------------------------------------

            /**
             * T-Matrix is not calculated for vertices if this flag is not set
             */
            void
            unset_t_matrix_flag()
            {
                MORIS_ASSERT(false,"Not implemented");
            }

//------------------------------------------------------------------------------

            /**
             * query if the T-Matrix is to be calculated
             */
            bool
            get_t_matrix_flag() const
            {
                MORIS_ASSERT(false,"Not implemented");
                return false;
            }

//------------------------------------------------------------------------------

            /**
             * returns the order of the element
             */
            Interpolation_Order
            get_interpolation_order() const
            {
                MORIS_ASSERT(false,"Not implemented");
                return Interpolation_Order::UNDEFINED;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
//------------------------------------------------------------------------------



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_CELL_STK_HPP_ */
