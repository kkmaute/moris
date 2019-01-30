/*
 * cl_XTK_Interface_Element.hpp
 *
 *  Created on: Jan 16, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_INTERFACE_ELEMENT_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_INTERFACE_ELEMENT_HPP_

#include "cl_Matrix.hpp"
#include "assert.hpp"

namespace xtk
{

class Interface_Element
{
public:
    Interface_Element():
        mElementIndexPair(0,0),
        mElementSideOrdinals(0,0)
    {

    }

    /*
     * Provide this interface element with the element index pair and corresponding side ordinals which these two elements
     * share
     */
    void
    set_element_pair_and_side_ordinal(moris::Matrix<moris::IndexMat> aElementIndexPair,
                                      moris::Matrix<moris::IndexMat> aElementSideOrdinals)
    {
        // Verify lengths
        MORIS_ASSERT(aElementIndexPair.numel() == 2,   "Provided aElementIndexPair vector needs to have length 2");
        MORIS_ASSERT(aElementSideOrdinals.numel() == 2,"Provided aElementSideOrdinals vector needs to have length 2");

        // Verify that we haven't already set these
        MORIS_ASSERT(mElementIndexPair.numel() == 0,"This interface element has already been given an element index pair");
        MORIS_ASSERT(mElementSideOrdinals.numel() == 0,"This interface element has already been given an element side ordinal pair");

        // Copy the vectors
        mElementIndexPair = aElementIndexPair.copy();
        mElementSideOrdinals = aElementSideOrdinals.copy();
    }


    moris::Matrix<moris::IndexMat> const &
    get_element_indices_pair() const
    {
        MORIS_ASSERT(mElementIndexPair.numel() == 2,"This interface element has not been given an element index pair");
        return mElementIndexPair;
    }

    moris::Matrix<moris::IndexMat> const &
    get_element_pair_side_ordinals() const
    {
        MORIS_ASSERT(mElementSideOrdinals.numel() == 2,"This interface element has not been given an element index pair");
        return mElementSideOrdinals;
    }


private:
    moris::Matrix<moris::IndexMat> mElementIndexPair;
    moris::Matrix<moris::IndexMat> mElementSideOrdinals;
};

}


#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_INTERFACE_ELEMENT_HPP_ */
