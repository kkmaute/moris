/*
 * cl_MTK_Vertex_STK.hpp
 *
 *  Created on: Sep 17, 2018
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_
#define PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_

#include "../cl_MTK_Vertex.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Cell.hpp"
#include "../cl_Mesh.hpp"
#include "fn_print.hpp"

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------
        class Vertex_STK: public Vertex
        {
        private:

            moris_id    mVertexId;
            moris_index mVertexInd;
            Mesh*   mSTKMeshData;


            // TODO: remove as functions are filled in
            // Dummy member vars
            Matrix< DDRMat > mTMatrix;
            moris::Cell< Vertex* > mDUMMYAdof;


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
        {

        };

//------------------------------------------------------------------------------

            /**
             * trivial constructor
             */
            Vertex_STK(){};

//------------------------------------------------------------------------------

            /**
             * Destructor, virtual
             */
            ~Vertex_STK(){};

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
        };
//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
    //------------------------------------------------------------------------------



#endif /* PROJECTS_MTK_SRC_STK_IMPL_CL_MTK_VERTEX_STK_HPP_ */
