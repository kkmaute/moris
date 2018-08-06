/*
 * cl_Nodal_Pdof_Host.hpp
 *
 *  Created on: Jul 15, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NODAL_PDOF_HOST_HPP_
#define SRC_FEM_CL_NODAL_PDOF_HOST_HPP_

#include "cl_Pdof_Host.hpp"

namespace moris
{
    namespace MSI
    {
    class Nodal_Pdof_Host //:public Pdof_Host
    {
    private:

    public:
        Nodal_Pdof_Host( )
        {
            //NodeId = aNodeId;
        };

        ~Nodal_Pdof_Host()
        {};
    };
    }
}



#endif /* SRC_FEM_CL_NODAL_PDOF_HOST_HPP_ */
