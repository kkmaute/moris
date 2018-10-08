/*
 * cl_MSI_Node_Proxy.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: schmidt
 */
#ifndef SRC_MSI_CL_NODE_OBJ_HPP_
#define SRC_MSI_CL_NODE_OBJ_HPP_

#include "cl_FEM_Node_Base.hpp"

namespace moris
{
    namespace MSI
    {
    class Node_Proxy : public moris::fem::Node_Base
    {
    private:
        moris::sint                mNodeId;
        moris::sint                mNodeOwner;
        Matrix< IdMat > mAdofIds;
        Matrix< IndexMat >  mAdofInd;
        Matrix< DDRMat >  mMatrix;
        Matrix< IdMat > mAdofOwningProcessor;


    public:
        Node_Proxy(){};

        Node_Proxy( const moris::luint             & aNodeId,
                  const Matrix< IdMat > & aAdofIds,
                  const Matrix< IndexMat > & aAdofInd,
                  const Matrix< DDRMat > & aMatrix,
                  const Matrix< IdMat > & aOwningProcessorList ) : mNodeId( aNodeId ),
                                                                            mAdofIds( aAdofIds ),
                                                                            mAdofInd( aAdofInd ),
                                                                            mMatrix ( aMatrix ),
                                                                            mAdofOwningProcessor( aOwningProcessorList )
        {};

        ~Node_Proxy() {};

        moris::sint get_id() const { return mNodeId; }
        moris::sint get_index() const { return mNodeId; }

        moris::sint get_owner() const { return mNodeOwner; }

        Matrix< IdMat > get_adof_ids() const { return mAdofIds; };

        Matrix< IndexMat > get_adof_indices() const { return mAdofInd; };

        const Matrix< DDRMat > * get_t_matrix() const { return & mMatrix; };

        Matrix< IndexMat > get_adof_owners() const { return mAdofOwningProcessor; };
    };
    }
}

#endif /* SRC_MSI_CL_NODE_OBJ_HPP_ */
