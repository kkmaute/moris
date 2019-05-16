/*
 * cl_MTK_Cell_Proxy.hpp
 *
 *  Created on: Apr 26, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_TEST_CL_MTK_CELL_PROXY_HPP_
#define PROJECTS_MTK_TEST_CL_MTK_CELL_PROXY_HPP_

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        class Cell_Proxy : public Cell
        {
//------------------------------------------------------------------------------
        public:
            moris::moris_id          mId;
            moris::moris_id          mIndex;
            moris::moris_id          mOwner;
            moris::uint              mNumVerts;
            moris::Cell< Vertex* >   mVertices;
            moris::uint              mSpatialDim   = 3;
            enum Geometry_Type       mGeometryType = Geometry_Type::UNDEFINED;
            enum Interpolation_Order mInterpOrder  = Interpolation_Order::UNDEFINED;

//------------------------------------------------------------------------------

            Cell_Proxy (){};

//------------------------------------------------------------------------------

            ~Cell_Proxy (){};

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of the cell
             *
             * @return moris_id ID
             */
            moris_id
            get_id() const
            {
                return mId;
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
                return mIndex;
            }

//------------------------------------------------------------------------------

            /**
             * tells how many vertices are connected to this cell
             */
            uint
            get_number_of_vertices() const
            {
                return mVertices.size();
            }

//------------------------------------------------------------------------------

            /**
             * returns the proc id of the owner of this cell
             * ( this information is needed for STK )
             */
            moris_id
            get_owner() const
            {
                return mOwner;
            }

//------------------------------------------------------------------------------

            /**
             * fills a moris::cell with pointers to connected vertices
             */
            moris::Cell< Vertex* >
            get_vertex_pointers() const
            {
                return mVertices;
            }

//------------------------------------------------------------------------------

            /**
             * returns a Mat with IDs of connected vertices
             */
            Matrix< IdMat >
            get_vertex_ids() const
            {
                Matrix<IdMat> tVertexIds(mVertices.size(),1);

                for(moris::uint  i = 0 ; i < this->get_number_of_vertices(); i++)
                {
                    tVertexIds(i) = mVertices(i)->get_id();
                }

                return tVertexIds;
            }

//------------------------------------------------------------------------------

            /**
             * returns a Mat with indices of connected vertices
             */
             Matrix< IndexMat >
            get_vertex_inds() const
            {
                 Matrix<IndexMat> tVertexInd(mVertices.size(),1);

                 for(moris::uint  i = 0 ; i < this->get_number_of_vertices(); i++)
                 {
                     tVertexInd(i) = mVertices(i)->get_index();
                 }

                 return tVertexInd;
            }

//------------------------------------------------------------------------------

            /**
             * returns a Mat of dimension
             * < number of vertices * number of dimensions >
             */
            Matrix< DDRMat >
            get_vertex_coords() const
            {
                size_t tNumVertices = this->get_number_of_vertices();
                Matrix< DDRMat > tVertexCoords(tNumVertices, mSpatialDim);
                for(size_t i = 0; i<tNumVertices; i++)
                {
                    tVertexCoords.set_row(i,mVertices(i)->get_coords());
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
                return mGeometryType;
            }

//------------------------------------------------------------------------------

            /**
             * returns the order of the element
             */
            Interpolation_Order
            get_interpolation_order() const
            {
                return mInterpOrder;
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace mtk */
}


#endif /* PROJECTS_MTK_TEST_CL_MTK_CELL_PROXY_HPP_ */
