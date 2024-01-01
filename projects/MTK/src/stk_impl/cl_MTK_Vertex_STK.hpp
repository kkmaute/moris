/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Vertex_STK.hpp
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
        class Vertex_STK: public Vertex
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
            Vertex_STK(moris_id aVertexId,
                       moris_index aVertexInd,
                       Mesh* aStkImplementation):
                           mVertexId(aVertexId),
                           mVertexInd(aVertexInd),
                           mSTKMeshData(aStkImplementation)
        {};

//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex_STK(){};

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Vertex_STK(){};

//------------------------------------------------------------------------------

            void
            set_vertex_interpolation(Vertex_Interpolation * aVertexInterpolation)
            {
                MORIS_ASSERT(aVertexInterpolation!=nullptr,"Null pointer passed as vertex interpolation");
                mVertexInterpolation = aVertexInterpolation;
            }
//------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            Matrix< DDRMat >
            get_coords() const
            {
                return mSTKMeshData->get_node_coordinate(mVertexInd);
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
                return mSTKMeshData->get_entity_owner( mVertexInd, EntityRank::NODE);
            }

//------------------------------------------------------------------------------

            Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex )
            {
                //MORIS_ERROR(0," Vertex interpolation not implemented");
                return mVertexInterpolation;
            }

            const Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) const
            {
                //MORIS_ERROR(0," Vertex interpolation not implemented - const");
                return mVertexInterpolation;
            }

//------------------------------------------------------------------------------
        };

        class Vertex_Core_STK: public Vertex
        {
        private:

            moris_id               mVertexId;
            moris_index            mVertexInd;
            Mesh*         mSTKMeshData;
            Vertex_Interpolation*  mVertexInterpolation = nullptr;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             *  constructor
             */
            Vertex_Core_STK(moris_id aVertexId,
                            moris_index aVertexInd,
                            Mesh* aSTKCore):
                           mVertexId(aVertexId),
                           mVertexInd(aVertexInd),
                           mSTKMeshData(aSTKCore)
        {};

//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex_Core_STK(){};

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Vertex_Core_STK(){};

//------------------------------------------------------------------------------

            void
            set_vertex_interpolation(Vertex_Interpolation * aVertexInterpolation)
            {
                MORIS_ASSERT(aVertexInterpolation!=nullptr,"Null pointer passed as vertex interpolation");
                mVertexInterpolation = aVertexInterpolation;
            }
//------------------------------------------------------------------------------

            /**
             * returns a moris::Matrix with node coordinates
             */
            Matrix< DDRMat >
            get_coords() const
            {
                return mSTKMeshData->get_node_coordinate(mVertexInd);
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
                return mSTKMeshData->get_entity_owner( mVertexInd, EntityRank::NODE);
            }

//------------------------------------------------------------------------------

            Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex )
            {
                //MORIS_ERROR(0," Vertex interpolation not implemented");
                return mVertexInterpolation;
            }

            const Vertex_Interpolation * get_interpolation( const uint aBSplineMeshIndex ) const
            {
                //MORIS_ERROR(0," Vertex interpolation not implemented - const");
                return mVertexInterpolation;
            }

            bool has_interpolation( const uint aBSplineMeshIndex )
            {
                return true;
            }

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------

#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_ */

