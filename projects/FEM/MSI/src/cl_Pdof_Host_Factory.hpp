/*
 * cl_Pdof_Host_Factory.hpp
 *
 *  Created on: Jul 15, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_PDOF_HOST_FACTORY_HPP_
#define SRC_FEM_CL_PDOF_HOST_FACTORY_HPP_

#include "cl_Nodal_Pdof_Host.hpp"

namespace moris
{
    namespace MSI
    {
    class Pdof_Host_Factory
    {
    private:
    public:
        Pdof_Host_Factory()
        {};

        ~Pdof_Host_Factory()
        {};

//        Pdof_Host * create_pdof_host( const moris::uint & aNodeId )
//        {
//            Pdof_Host * tPdofHost;
//
////            switch(0)
////            {
////            case (0):
////                tPdofHost = new Nodal_Pdof_Host( aNodeId );
////                break;
////            default:
////                MORIS_ASSERT( false, "No Pdof_Host type specified." );
////                break;
////            }
//            return tPdofHost;
//        };
    };
    }
}

#endif /* SRC_FEM_CL_PDOF_HOST_FACTORY_HPP_ */
