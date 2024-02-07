/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_Proxy.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_

#include "cl_MTK_Vertex.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Vertex_Interpolation_STK.hpp"
#include "cl_MTK_Mesh_Core.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertex_Proxy: public Vertex
        {
        private:

            moris_id               mVertexId;
            moris_index            mVertexInd;
            Mesh*                  mSTKMeshData;
            Vertex_Interpolation*  mVertexInterpolation = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Vertex_Proxy(moris_id aVertexId ):
                           mVertexId(aVertexId)
        {};

//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex_Proxy(){};

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Vertex_Proxy(){};

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
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_ */

