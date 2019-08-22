/*
 * cl_MTK_Vertex_Proxy.hpp
 *
 *  Created on: Apr 26, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_TEST_CL_MTK_VERTEX_PROXY_HPP_
#define PROJECTS_MTK_TEST_CL_MTK_VERTEX_PROXY_HPP_

#include "catch.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertex_Proxy: public Vertex
        {
        public:

            moris_id               mVertexId;
            moris_index            mVertexInd;
            moris_id               mOwner = 0;
            Matrix< DDRMat >       mVertexCoord;
            Vertex_Interpolation*  mVertexInterpolation = nullptr;



//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Vertex_Proxy(){};

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Vertex_Proxy(){};

//------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            Matrix< DDRMat >
            get_coords() const
            {
                return mVertexCoord;
            }

//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            moris_id
            get_id() const
            {
                return mVertexId;
            }


//------------------------------------------------------------------------------

            /**
             * returns the domain wide id of this vertex
             */
            moris_index
            get_index() const
            {
                return mVertexInd;
            }

//------------------------------------------------------------------------------

            /**
             * returns the id of the proc that owns this vertex
             */
            moris_id
            get_owner() const
            {
                return mOwner;
            }

//------------------------------------------------------------------------------
            Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex )
            {
                //MORIS_ERROR(0," Vertex interpolation not implemented");
                return mVertexInterpolation;
            }

//------------------------------------------------------------------------------
            const Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) const
            {
                //MORIS_ERROR(0," Vertex interpolation not implemented - const");
                return mVertexInterpolation;
            }


//------------------------------------------------------------------------------
        };

    } /* namespace mtk */
}




#endif /* PROJECTS_MTK_TEST_CL_MTK_VERTEX_PROXY_HPP_ */
