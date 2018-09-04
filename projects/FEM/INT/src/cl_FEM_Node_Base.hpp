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

//------------------------------------------------------------------------------

        class Node_Base
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * destructor
             */
            virtual ~Node_Base(){};

//------------------------------------------------------------------------------

            /**
             * returns the T-Matrix of this node
             */

            virtual const moris::Mat < moris::real > * get_t_matrix() const =0;

//------------------------------------------------------------------------------

            virtual moris::Mat < sint > get_adof_ids() const = 0;

//------------------------------------------------------------------------------

            virtual moris::Mat < sint > get_adof_indices() const = 0;

//------------------------------------------------------------------------------

            /**
             * returns the proc owners of the IDs of this node
             */

            virtual moris::Mat < moris::uint > get_adof_owners() const = 0;

//------------------------------------------------------------------------------

            /**
             * get the ID of this node
             *
             * @param[ in ] aID  id for this node
             */

            virtual sint get_id() const = 0;

//------------------------------------------------------------------------------


            /**
             * get the Ind of this node
             *
             * @param[ in ] aInd  ind for this node
             */

            virtual sint get_index() const = 0;
        };
    } /* namespace fem */
} /* namespace moris */



#endif /* PROJECTS_FEM_MSI_SRC_CL_MSI_NODE_BASE_HPP_ */
