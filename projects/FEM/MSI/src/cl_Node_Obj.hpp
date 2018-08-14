/*
 * cl_Node_Obj.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: schmidt
 */

#ifndef SRC_FEM_CL_NODE_OBJ_HPP_
#define SRC_FEM_CL_NODE_OBJ_HPP_

#include "linalg.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_Pdof_Host.hpp"

namespace moris
{
    namespace MSI
    {
    // FIXME change to FEM node as soon as FEM Node exists
    class Node_Obj : public mtk::Vertex
    {
    private:
        moris::luint                mNodeId;
        moris::Mat < moris::sint > mAdofs;
        moris::Mat < moris::real > mMatrix;
        moris::Mat < moris::uint > mAdofOwningProcessor;


    public:
        //moris::Mat < moris::uint > mAdofs;

        Node_Obj()
        {
        };

        Node_Obj( const moris::luint             & aNodeId,
                  const moris::Mat< moris::sint> & aAdofs,
                  const moris::Mat< moris::real> & aMatrix,
                  const moris::Mat< moris::uint> & aOwningProcessorList ) : mNodeId( aNodeId ),
                                                                            mAdofs( aAdofs ),
                                                                            mMatrix ( aMatrix ),
                                                                            mAdofOwningProcessor( aOwningProcessorList )
        {

        };

        ~Node_Obj()
        {};

        moris::luint get_id() const
        {
            return mNodeId;
        }

        moris::Mat < moris::sint > get_adof_ids() const
        {
            return mAdofs;
        };

        const moris::Mat < moris::real > * get_t_matrix() const
        {
            return & mMatrix;
        };

        moris::Mat < moris::uint > get_adof_owners() const
        {
            return mAdofOwningProcessor;
        };

        //------------------------------------------------------------------------------

        moris::Mat< moris::real >
        get_coords() const
        {
            MORIS_ERROR( false, "get_coords() not available for node object.");
            return moris::Mat< moris::real >(0,0);
        }

        //------------------------------------------------------------------------------

        moris::Cell< moris::mtk::Vertex* >
        get_adof_pointers()
        {
            MORIS_ERROR( false, "get_adof_pointers() not available for node object.");
            return moris::Cell< moris::mtk::Vertex* >(0);
        }

        moris::uint
        get_owner() const
        {
            MORIS_ERROR( false, "get_owner() not available for node object.");
            return MORIS_UINT_MAX;
        }

    };
    }
}



#endif /* SRC_FEM_CL_NODE_OBJ_HPP_ */
