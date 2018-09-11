/*
 * cl_MSI_Node_Obj.hpp
 *
 *  Created on: Jul 24, 2018
 *      Author: schmidt
 */
#ifndef SRC_MSI_CL_NODE_OBJ_HPP_
#define SRC_MSI_CL_NODE_OBJ_HPP_

#include "linalg.hpp"
#include "cl_FEM_Node_Base.hpp"

namespace moris
{
    namespace MSI
    {
    class Node_Obj : public moris::fem::Node_Base
    {
    private:
        moris::sint                mNodeId;
        moris::Mat < moris::sint > mAdofIds;
        moris::Mat < moris::sint > mAdofInd;
        moris::Mat < moris::real > mMatrix;
        moris::Mat < moris::uint > mAdofOwningProcessor;


    public:
        Node_Obj(){};

        Node_Obj( const moris::luint             & aNodeId,
                  const moris::Mat< moris::sint> & aAdofIds,
                  const moris::Mat< moris::sint> & aAdofInd,
                  const moris::Mat< moris::real> & aMatrix,
                  const moris::Mat< moris::uint> & aOwningProcessorList ) : mNodeId( aNodeId ),
                                                                            mAdofIds( aAdofIds ),
                                                                            mAdofInd( aAdofInd ),
                                                                            mMatrix ( aMatrix ),
                                                                            mAdofOwningProcessor( aOwningProcessorList )
        {};

        ~Node_Obj() {};

        moris::sint get_id() const { return mNodeId; }
        moris::sint get_index() const { return mNodeId; }

        moris::Mat < moris::sint > get_adof_ids() const { return mAdofIds; };

        moris::Mat < moris::sint > get_adof_indices() const { return mAdofInd; };

        const moris::Mat < moris::real > * get_t_matrix() const { return & mMatrix; };

        moris::Mat < moris::uint > get_adof_owners() const { return mAdofOwningProcessor; };
    };
    }
}

#endif /* SRC_MSI_CL_NODE_OBJ_HPP_ */
