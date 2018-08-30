/*
 * cl_FEM_Node_Base.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_FEM_MSI_SRC_CL_MSI_NODE_BASE_HPP_
#define PROJECTS_FEM_MSI_SRC_CL_MSI_NODE_BASE_HPP_

#include "typedefs.hpp"
#include "cl_Mat.hpp"

namespace moris
{
    namespace fem
    {
        class Node_Base
        {
        public:

            /**
             * destructor
             */
            virtual ~Node_Base(){};

            /**
             * returns the owner of this node
             */
//            virtual moris::uint get_owner() const
//            {
//                MORIS_ERROR( false, "get_index() not available for node object.");
//                return 0;
//            }

//------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this node
             */

            virtual const moris::Mat < moris::real > * get_t_matrix() const =0;


//------------------------------------------------------------------------------

            virtual moris::Mat < moris::sint > get_adof_ids() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this node
             */

            virtual moris::Mat < moris::uint > get_adof_owners() const = 0;

//------------------------------------------------------------------------------

            /**
             * set the ID of this node
             *
             * @param[ in ] aID  id for this node
             */
//            virtual void set_id( const luint & aID )
//            {
//                MORIS_ERROR( false, "get_owner() not available for node object.");
//            }

//------------------------------------------------------------------------------

            /**
             * get the ID of this node
             *
             * @param[ in ] aID  id for this node
             */

            virtual moris::luint get_id() const = 0;

//------------------------------------------------------------------------------

            /**
             * set the ID of this node
             *
             * @param[ in ] aID  id for this node
             */
//            virtual void set_index( const sint & aIndex )
//            {
//                MORIS_ERROR( false, "get_coords() not available for node object.");
//            }

//------------------------------------------------------------------------------

            /**
             * get the Index of this node
             *
             * @param[ in ] aID  id for this node
             */
//            virtual moris:: sint get_index() const
//            {
//                return 0;
//            }
        };
    } /* namespace fem */
} /* namespace moris */



#endif /* PROJECTS_FEM_MSI_SRC_CL_MSI_NODE_BASE_HPP_ */
