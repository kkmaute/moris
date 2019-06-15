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

namespace moris
{
    namespace fem
    {
        class Node_Base
        {

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
             * @param[ in ] aID  id for this node
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
             * @param[ in ] aID  id for this node
             */

            virtual moris_index
            get_index() const
            {
                MORIS_ERROR( false, "Entered virtual function Node_Base::get_index()" );
                return gNoIndex;
            }

//------------------------------------------------------------------------------

        };
    } /* namespace fem */
} /* namespace moris */



#endif /* PROJECTS_FEM_SRC_CL_FEM_NODE_BASE_HPP_ */
