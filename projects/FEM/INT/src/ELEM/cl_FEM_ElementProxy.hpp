/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_ElementProxy.hpp
 *
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_ELEMENTPROXY_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_ELEMENTPROXY_HPP_

#include "cl_Matrix.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace fem
    {

        class ElementProxy : public mtk::Cell
        {
            //------------------------------------------------------------------------------

          private:
          
            moris::Cell< mtk::Vertex* >   mNodeList;
            mtk::Geometry_Type       mGeometryType       = mtk::Geometry_Type::UNDEFINED;
            mtk::Interpolation_Order mInterpolationOrder = mtk::Interpolation_Order::UNDEFINED;

            //------------------------------------------------------------------------------

          public:

            //------------------------------------------------------------------------------

            ElementProxy( moris::Cell< mtk::Vertex* > aNodeList,
                    mtk::Geometry_Type           aGeometryType,
                    mtk::Interpolation_Order     aInterpolationOrder )
            {
                mNodeList           = aNodeList;
                mGeometryType       = aGeometryType;
                mInterpolationOrder = aInterpolationOrder;
            }

            //------------------------------------------------------------------------------

            ~ElementProxy(){};

            //------------------------------------------------------------------------------

            moris_id
            get_id() const
            {
                MORIS_ERROR( false, " ElementProxy::get_id - not implemented. " );
                moris_id tId = 0;
                return tId;
            };

            //------------------------------------------------------------------------------

            moris_index
            get_index() const
            {
                MORIS_ERROR( false, " ElementProxy::get_index - not implemented. " );
                moris_index mIndex = 0;
                return mIndex;
            };

            //------------------------------------------------------------------------------

            uint
            get_number_of_vertices() const
            {
                return mNodeList.size();
            };

            //------------------------------------------------------------------------------

            moris_id
            get_owner() const
            {
                MORIS_ERROR( false, " ElementProxy::get_owner - not implemented " );
                moris_id tId = 0;
                return tId;
            };

            //------------------------------------------------------------------------------

            moris::Cell< mtk::Vertex* >
            get_vertex_pointers() const
            {
                return mNodeList;
            };

            //------------------------------------------------------------------------------

            // TODO MESHCLEANUP
            void
            remove_vertex_pointer( moris_index aIndex )
            {
                std::cout << "In MTK Cell Proxy" << std::endl;
            }


            //------------------------------------------------------------------------------

            Matrix< IdMat >
            get_vertex_ids() const
            {
                MORIS_ERROR( false, " ElementProxy::get_vertex_ids - not available for this element. " );
                return Matrix< IdMat >( 0, 0 );
            };

            //------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_vertex_inds() const
            {
                MORIS_ERROR( false, " ElementProxy::get_vertex_ids - not available for this element. " );
                return Matrix< IndexMat >( 0, 0 );
            };

            //------------------------------------------------------------------------------

            Matrix< DDRMat >
            get_vertex_coords() const
            {
                //            MORIS_ERROR( false, "get_vertex_coords() not available for this element " );
                //            return Matrix< DDRMat > (0,0);
                // FIXME: only for 2d
                uint tNumberOfVertices = this->get_number_of_vertices();

                Matrix< DDRMat > tVerticesCoords( tNumberOfVertices, 2 );
                for ( uint i = 0; i < tNumberOfVertices; i++ )
                {
                    Matrix< DDRMat > tVertexCoords = mNodeList( i )->get_coords();
                    for ( uint j = 0; j < 2; j++ )
                    {
                        tVerticesCoords( i, j ) = tVertexCoords( j );
                    }
                }
                return tVerticesCoords;
            };

            //------------------------------------------------------------------------------

            mtk::Geometry_Type
            get_geometry_type() const
            {
                return mGeometryType;
            };

            //------------------------------------------------------------------------------

            mtk::Interpolation_Order
            get_interpolation_order() const
            {
                // MORIS_ERROR( false, " ElementProxy::get_interpolation_order - not available for this element. " );
                return mInterpolationOrder;
            };

            //------------------------------------------------------------------------------

            mtk::Integration_Order
            get_integration_order() const
            {
                MORIS_ERROR( false, " ElementProxy::get_integration_order - not available for this element. " );
                return mtk::Integration_Order::UNDEFINED;
            };

            //------------------------------------------------------------------------------

        }; // class ElementProxy

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_ELEMENTPROXY_HPP_ */
