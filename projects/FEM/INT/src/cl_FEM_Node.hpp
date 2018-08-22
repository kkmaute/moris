/*
 * cl_FEM_Node.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: messe
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_NODE_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_NODE_HPP_

#include "typedefs.hpp"           //MRS/COR/src
#include "cl_MTK_Vertex.hpp"      //MTK/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Node
        {
            const mtk::Vertex & mVertex;
            sint                mID;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * constructor
             */
            Node( const mtk::Vertex & aVertex ) : mVertex( aVertex )
            {
                    mID = aVertex->get_id();
            }

//------------------------------------------------------------------------------

            /**
             * destructor
             */
            ~Node(){};

//------------------------------------------------------------------------------

            /**
             * returns the id of this node
             */
            auto
            get_id() const -> decltype( mID )
            {
                return mID;
            }

//------------------------------------------------------------------------------

            /**
             * returns the owner of this node
             */
            auto
            get_owner() const -> decltype( mVertex.get_owner() )
            {
                return mVertex.get_owner();
            }

//------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this node
             */
            auto
            get_t_matrix() const -> decltype( mVertex.get_t_matrix() )
            {
                return mVertex.get_t_matrix();
            }

//------------------------------------------------------------------------------

            /**
             * returns the B-Spline IDs of this node
             */
            auto
            get_adof_ids() const  -> decltype( mVertex.get_adof_ids() )
            {
                return mVertex.get_adof_ids();
            }

//------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this node
             */
            auto
            get_adof_owners() const -> decltype( mVertex.get_adof_owners() )
            {
                return mVertex.get_adof_owners();
            }

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_NODE_HPP_ */
