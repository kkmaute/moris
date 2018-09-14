/*
 * cl_Side_Set_Input.hpp
 *
 *  Created on: Aug 30, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_SIDE_SET_INPUT_HPP_
#define SRC_MESH_CL_MESH_SIDE_SET_INPUT_HPP_

#include "assert/fn_xtk_assert.hpp"

#include "containers/cl_XTK_Cell.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include <memory>

using namespace xtk;

namespace mesh
{
template<typename Integer, typename Integer_Matrix>
class Side_Set_Input
{
public:
    Side_Set_Input(Integer const & aNumSides,
                   Integer const & aMaxStringLength)
    {
        mElementIdAndSideOrdinal.reserve(aNumSides); // 2 for the pair
        mSideSetName.reserve(aMaxStringLength);
        mFaceNodes = moris::Matrix< Integer_Matrix >(1,1);
    }



    void add_element_id_and_side_ordinal(moris::Matrix< Integer_Matrix > const & aElementIds,
                                         moris::Matrix< Integer_Matrix > const & aSideOrdinals)
    {
        Integer tNumElementIds = aElementIds.n_cols();
        Integer tNumSideOrdinals = aSideOrdinals.n_cols();
        XTK_ASSERT(tNumElementIds==tNumSideOrdinals, "The number of element IDs needs to be equal to the number of side ordinals ");

        for(Integer iAdd = 0; iAdd<tNumElementIds; iAdd++)
        {
            mElementIdAndSideOrdinal.push_back(std::pair<Integer,Integer>(aElementIds(0,iAdd),aSideOrdinals(0,iAdd)));
        }
    }




    void
    add_element_id_and_side_ordinal_with_face_nodes(moris::Matrix< Integer_Matrix > const & aElementIds,
                                                    moris::Matrix< Integer_Matrix > const & aSideOrdinals,
                                                    moris::Matrix< Integer_Matrix > const & aFaceNodes)
     {
         Integer tNumElementIds   = aElementIds.n_cols();
         Integer tNumSideOrdinals = aSideOrdinals.n_cols();
         XTK_ASSERT(tNumElementIds == tNumSideOrdinals, "The number of element IDs needs to be equal to the number of side ordinals ");

         // Resize mFaceNodes
         Integer tOldSize = mElementIdAndSideOrdinal.size();
         Integer tNumNodesPerFace = aFaceNodes.n_cols();
         mFaceNodes->resize(tOldSize+tNumElementIds,tNumNodesPerFace);


         for(Integer iAdd = 0; iAdd<tNumElementIds; iAdd++)
         {
             mElementIdAndSideOrdinal.push_back(std::pair<Integer,Integer>(aElementIds(0,iAdd),aSideOrdinals(0,iAdd)));
             replace_row( iAdd, aFaceNodes, tOldSize+iAdd, *mFaceNodes);
         }

     }


    void add_element_id_and_side_ordinal(Integer const & aElementId,
                                         Integer const & aSideOrdinal,
                                         Integer const & aRowIndex)
    {
        mElementIdAndSideOrdinal.push_back(std::pair<Integer,Integer>(aElementId,aSideOrdinal));
    }

    void add_element_id_and_side_ordinal_with_face_nodes(Integer const & aElementId,
                                                         Integer const & aSideOrdinal,
                                                         Integer const & aRowIndex,
                                                         moris::Matrix< Integer_Matrix > const & aFaceNodes)
    {
        mElementIdAndSideOrdinal.push_back(std::pair<Integer,Integer>(aElementId,aSideOrdinal));

        // Resize mFaceNodes
        // FIXME: DYNAMIC ALLOCATION ALLOCATE A LOT AND SHRINK TO FIT
         Integer tOldSize = mElementIdAndSideOrdinal.size();
         Integer tNumNodesPerFace = aFaceNodes.n_cols();
         mFaceNodes->resize(tOldSize+1,tNumNodesPerFace);

         replace_row( aRowIndex, aFaceNodes, tOldSize, *mFaceNodes);

    }


    void set_side_set_name(std::string const & aSideSetName)
    {
        XTK_ASSERT(mSideSetName.length()==0,"Side set name already set, This would overwrite the previous name");
        mSideSetName = aSideSetName;
    }


    Integer get_num_of_sides() const
    {
        return mElementIdAndSideOrdinal.size();
    }

    Integer get_num_side_sets() const
    {
        return 1;
    }

    std::string const & get_side_set_name() const
    {
        return mSideSetName;
    }


    Integer const & get_element_id(Integer const & aPairIndex) const
    {
        XTK_ASSERT(aPairIndex< mElementIdAndSideOrdinal.size(),"Pair index is out of bounds");
        return mElementIdAndSideOrdinal(aPairIndex).first;
    }

    Integer const & get_side_ordinal(Integer aPairIndex) const
    {
        XTK_ASSERT(aPairIndex< mElementIdAndSideOrdinal.size(),"Pair index is out of bounds");
        return mElementIdAndSideOrdinal(aPairIndex).second;
    }

    moris::Matrix< Integer_Matrix > const & get_side_nodes(Integer const & aPairIndex) const
        {
            return *mFaceNodes(aPairIndex);
        }

    bool has_sides() const
    {
        if(mElementIdAndSideOrdinal.size()!=0)
        {
            return true;
        }

        else
            return false;
    }

    void print()
    {
        std::cout<<" Part Name: "<< mSideSetName <<std::endl;

        Integer tNumPairs = mElementIdAndSideOrdinal.size();
        for(Integer iPair = 0; iPair<tNumPairs; iPair++)
        {
            std::cout<<"Element Id: "<< mElementIdAndSideOrdinal(iPair).first<< " Ordinal: " << mElementIdAndSideOrdinal(iPair).second<< " | Nodes: ";

            if(mFaceNodes.size() != 0 )
            {
            for(Integer iNode = 0; iNode<mFaceNodes(iPair)->n_cols(); iNode++)
            {
                std::cout<<(*mFaceNodes(iPair))(0,iNode)<<" ";
            }
            std::cout<<std::endl;
            }
        }
    }

private:
    std::string mSideSetName;
    Cell<std::pair<Integer,Integer>> mElementIdAndSideOrdinal;
    moris::Matrix< Integer_Matrix > mFaceNodes;
};

}


#endif /* SRC_MESH_CL_MESH_SIDE_SET_INPUT_HPP_ */
