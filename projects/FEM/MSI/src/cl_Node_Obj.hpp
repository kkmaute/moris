/*
 * cl_Node_Obj.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_NODE_OBJ_HPP_
#define SRC_FEM_CL_NODE_OBJ_HPP_

#include "linalg.hpp"
#include "cl_Pdof_Host.hpp"

namespace moris
{
    namespace MSI
    {
    class Node_Obj
    {
    private:
        moris::uint                mNodeId;
        moris::Mat < moris::sint > mAdofs;
        moris::Mat < moris::real > mMatrix;

    public:
        //moris::Mat < moris::uint > mAdofs;

        Node_Obj()
        {
        };

        Node_Obj( const moris::uint              & aNodeId,
                  const moris::Mat< moris::sint> & aAdofs,
                  const moris::Mat< moris::real> & aMatrix ) : mNodeId( aNodeId ),
                                                               mAdofs( aAdofs ),
                                                               mMatrix ( aMatrix )
        {

        };

        ~Node_Obj()
        {};

        const moris::uint get_node_id()
        {
            return mNodeId;
        }

        const moris::Mat < moris::sint > & get_adofs()
        {
            return mAdofs;
        };

        moris::Mat < moris::real > * get_Tmatrix()
        {
            return & mMatrix;
        };

    };
    }
}



#endif /* SRC_FEM_CL_NODE_OBJ_HPP_ */
