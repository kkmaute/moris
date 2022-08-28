/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Node_Proxy.hpp
 *
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

        Node_Proxy( const moris::luint             & aNodeId ) : mNodeId( aNodeId )
        {};

        ~Node_Proxy() {};

        moris::sint get_id() const { return mNodeId; }
        moris::sint get_index() const { return mNodeId; }

        moris::sint get_owner() const { return mNodeOwner; }

        Matrix< IdMat > get_adof_ids( const uint aOrder ) const { return mAdofIds; };

        Matrix< IndexMat > get_adof_indices( const uint aOrder ) const { return mAdofInd; };

        const Matrix< DDRMat > * get_t_matrix( const uint aOrder ) const { return & mMatrix; };

        Matrix< IndexMat > get_adof_owners( const uint aOrder ) const { return mAdofOwningProcessor; };
    };
    }
}

#endif /* SRC_MSI_CL_NODE_OBJ_HPP_ */

