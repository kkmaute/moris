/*
 * cl_Mesh_Node_Set_Input.hpp
 *
 *  Created on: Oct 2, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_NODE_SET_INPUT_HPP_
#define SRC_MESH_CL_MESH_NODE_SET_INPUT_HPP_

#include "containers/cl_XTK_Cell.hpp"
#include "assert/fn_xtk_assert.hpp"

namespace mesh
{



template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Node_Set_Input
{

public:

    Node_Set_Input(Integer const & aNumNodes,
                   Integer const & aNumNodeSets,
                   Integer const & aMaxStringLength,
                   Integer const & aModelDimension = 3,
                   Integer const & aNumIntegerScalarFields = 0,
                   Integer const & aNumIntegerVectorFields = 0,
                   Integer const & aNumRealScalarFields = 0,
                   Integer const & aNumRealVectorFields = 0):
                       mHasFields(false)
    {
        mNodeIds.reserve(aNumNodes);
        mNodeSetNames.reserve(aNumNodeSets*aMaxStringLength);

        if(aNumIntegerScalarFields !=0 || aNumIntegerVectorFields !=0 || aNumRealScalarFields !=0 || aNumRealVectorFields != 0)
        {
            mHasFields = true;
        }

        // Allocate fields memory
        mRealFieldNames.reserve((aNumRealScalarFields+aNumRealVectorFields)*aMaxStringLength);
        mRealFieldData.reserve(aModelDimension*aNumRealVectorFields + aNumRealScalarFields );
        mIntegerFieldNames.reserve((aNumIntegerScalarFields+aNumIntegerVectorFields)*aMaxStringLength);
        mIntegerFieldData.reserve(aModelDimension*aNumIntegerVectorFields + aNumIntegerScalarFields );
    }


    void
    add_node_set_names(xtk::Cell<std::string> const & aNodeSetNames)
    {
        mNodeSetNames.append(aNodeSetNames);
    }

    void
    add_node_set_name(std::string const & aNodeSetName)
    {
        mNodeSetNames.push_back(aNodeSetName);
    }

    void
    add_node_ids(xtk::Cell<Integer> const & aNodeIds)
    {
        mNodeIds.append(aNodeIds);
    }

    void
    add_node_id(Integer const & aNodeId)
    {
        mNodeIds.push_back(aNodeId);
    }

    void
    add_node_ids(moris::Matrix< Integer_Matrix > const & aNodeIds)
    {
        for(Integer i = 0; i<aNodeIds.n_cols(); i++)
        {
            mNodeIds.push_back(aNodeIds(0,i));
        }
    }

    xtk::Cell<std::string> const &
    get_node_set_names() const
    {
        return mNodeSetNames;
    }

    Integer const &
    get_node_id(Integer const & aIndex) const
    {
        return mNodeIds(aIndex);
    }

    xtk::Cell<Integer> const &
    get_node_ids() const
    {
        return mNodeIds;
    }

    bool has_node_parts()
    {
        if(mNodeSetNames.size() == 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    // Field information

    void
    add_integer_field_name_and_data(std::string const & aIntegerFieldName,
                                    xtk::Cell<moris::Matrix< Integer_Matrix >> const & aIntegerFieldData)
    {

        XTK_ASSERT(aIntegerFieldData.size() == get_num_nodes_in_node_set(),"Number of nodes in node set does not match size of data provided");
        mIntegerFieldNames.push_back(aIntegerFieldName);
        mIntegerFieldData.push_back(aIntegerFieldData);

    }


    Integer
    get_num_integer_fields() const
    {
        return mIntegerFieldNames.size();
    }

    std::string const &
    get_integer_field_name(Integer const & aFieldIndex) const
    {
        XTK_ASSERT(get_num_integer_fields()>aFieldIndex,"Requested field index out of bounds");

        return mIntegerFieldNames(aFieldIndex);
    }

    moris::Matrix< Integer_Matrix > const &
    get_integer_field_data(Integer const & aFieldIndex,
                           Integer const & aNodeIndex) const
    {
        XTK_ASSERT(get_num_integer_fields()>aFieldIndex,"Requested field index out of bounds");
        XTK_ASSERT(get_num_nodes_in_node_set()>aNodeIndex,"Requested node index out of bounds");

        return (*mIntegerFieldData(aFieldIndex)(aNodeIndex));
    }


    void
    add_real_field_name_and_data(std::string const & aIntegerFieldName,
                                 xtk::Cell<moris::Matrix< Real_Matrix >> const & aRealFieldData)
    {

        XTK_ASSERT(aRealFieldData.size() == get_num_nodes_in_node_set(),"Number of nodes in node set does not match size of data provided");
        mRealFieldNames.push_back(aIntegerFieldName);
        mRealFieldData.push_back(aRealFieldData);

    }


    Integer
    get_num_real_fields() const
    {
        return mRealFieldNames.size();
    }

    std::string const &
    get_real_field_name(Integer const & aFieldIndex) const
    {
        XTK_ASSERT(get_num_integer_fields()>aFieldIndex,"Requested field index out of bounds");

        return mRealFieldNames(aFieldIndex);
    }

    moris::Matrix< Real_Matrix > const &
    get_real_field_data(Integer const & aFieldIndex,
                        Integer const & aNodeIndex) const
    {
        XTK_ASSERT(get_num_real_fields()>aFieldIndex,"Requested field index out of bounds");
        XTK_ASSERT(get_num_nodes_in_node_set()>aNodeIndex,"Requested node index out of bounds");

        return mRealFieldData(aFieldIndex)(aNodeIndex)  ;
    }


    Integer
    get_num_nodes_in_node_set() const
    {
        return mNodeIds.size();
    }

    void
    print() const
    {
        std::cout<<"Parts: ";
        for(Integer iNames = 0; iNames<mNodeSetNames.size(); iNames++)
        {
            std::cout<< mNodeSetNames(iNames) << "  ";
        }
        std::cout<<std::endl;

        std::cout<<"Node Ids: "<<std::endl;

        for(Integer iNodes = 0; iNodes<mNodeIds.size(); iNodes++)
        {
            std::cout<<mNodeIds(iNodes)<<std::endl;
        }
    }

private:

    // Node Ids and Sets to place these nodes into
    xtk::Cell<Integer> mNodeIds;
    xtk::Cell<std::string> mNodeSetNames;

    // Fields (indexed same as the node ids above)
    bool mHasFields;

    // Real type fields
    xtk::Cell<std::string> mRealFieldNames;
    // Outer cell corresponds to field
    // Inner cell corresponds to node indexed the same as the ids)
    // Matrix corresponds to the field data
    xtk::Cell<xtk::Cell<moris::Matrix< Real_Matrix >>> mRealFieldData;

    // Integer type fields
    xtk::Cell<std::string> mIntegerFieldNames;

    // Outer cell corresponds to field
    // Inner cell corresponds to node indexed the same as the ids)
    // Matrix corresponds to the field data
    xtk::Cell<xtk::Cell<moris::Matrix< Integer_Matrix >>> mIntegerFieldData;



};

}



#endif /* SRC_MESH_CL_MESH_NODE_SET_INPUT_HPP_ */
