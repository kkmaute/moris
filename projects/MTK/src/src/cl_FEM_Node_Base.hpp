/*
 * cl_FEM_Node_Base.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_SRC_CL_FEM_NODE_BASE_HPP_
#define PROJECTS_FEM_SRC_CL_FEM_NODE_BASE_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Vertex.hpp"      //MTK/src

namespace moris
{
    namespace fem
    {
        class Node_Base
        {
        private:
            moris_index mNodeIndex;
            moris_index mNodeId;

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * destructor
             */
            virtual
            ~Node_Base(){};

//------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this node
             */

            virtual const Matrix< DDRMat > *
            get_t_matrix(  const uint aOrder ) const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_t_matrix()" );
                return nullptr;
            }

//------------------------------------------------------------------------------

            virtual Matrix< IdMat >
            get_adof_ids( const uint aOrder ) const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_adof_ids()" );
                return  Matrix< IdMat >(0,0);
            }

//------------------------------------------------------------------------------

            virtual Matrix< IndexMat >
            get_adof_indices(  const uint aOrder ) const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_adof_indices()" );
                return  Matrix< IndexMat >(0,0);
            }

//------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this node
             */

            virtual Matrix< IdMat >
            get_adof_owners(  const uint aOrder ) const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_adof_owners()" );
                return  Matrix< IdMat >(0,0);
            }

//------------------------------------------------------------------------------

            /**
             * get the ID of this node
             *
             */

            virtual moris_id
            get_id() const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_id()" );
                return gNoID;
            }

//------------------------------------------------------------------------------

            /**
             * get the ID of this node
             *
             */

            virtual moris_index
            get_index() const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_index()" );
                return gNoIndex;
            }

//------------------------------------------------------------------------------

            /**
             * set the ID of this node
             *
             * @param[ in ] aID  id for this node
             */

            void
            set_id( const moris_id aId )
            {
                mNodeId = aId;
            }

//------------------------------------------------------------------------------

            /**
             * set the index of this node
             *
             * @param[ in ] aIndex  id for this node
             */

            void
            set_index( const moris_id aIndex )
            {
                mNodeIndex = aIndex;
            }

//------------------------------------------------------------------------------

            /**
             * set the index of this node
             *
             * @param[ in ] aIndex  id for this node
             */

            virtual bool
            id_owned(  )
            {
                MORIS_ERROR( false, "Enterd virtual function Node_Base::id_owned()" );
                return false;
            }

//------------------------------------------------------------------------------
            /**
             * get vertex coordinates (if relevant)
             * @param[ out ] aIndex  id for this node
             */
            virtual void get_vertex_coords( Matrix< DDRMat > & aVertexCoords )
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_vertex_coords()" );
            }

//------------------------------------------------------------------------------

        };
    } /* namespace fem */
} /* namespace moris */



#endif /* PROJECTS_FEM_SRC_CL_FEM_NODE_BASE_HPP_ */
