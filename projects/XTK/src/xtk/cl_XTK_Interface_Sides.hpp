/*
 * cl_XTK_Interface_Sides.hpp
 *
 *  Created on: Sep 20, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_INTERFACE_SIDES_HPP_
#define SRC_XTK_CL_XTK_INTERFACE_SIDES_HPP_

#include "linalg/cl_XTK_Matrix_Base.hpp"
namespace xtk
{

template <typename Integer>
class Interface_Sides
{
public:
    Interface_Sides(Integer  const & aMaxNumInterfaceSides)
    {
        mInterfaceSideIndices.reserve(aMaxNumInterfaceSides);
    }

    /*
     * This add the child mesh local index for a side which lives on the interface
     */
    template<typename Integer_Matrix>
    void
    add_interface_side_with_side_index(moris::Matrix< Integer_Matrix > const & aSideIndices)
    {
          Integer tNumSidesToAdd = aSideIndices.n_cols();

          for(Integer i = 0; i < tNumSidesToAdd; i++)
          {
              mInterfaceSideIndices.push_back(aSideIndices(0,i));
          }
    }


    Cell<Integer> const &
    get_interface_sides() const
    {
        return mInterfaceSideIndices;
    }

    Integer
    get_num_interface_sides() const
    {
        return mInterfaceSideIndices.size();
    }

private:
    Cell<Integer> mInterfaceSideIndices;
};

}
#endif /* SRC_XTK_CL_XTK_INTERFACE_SIDES_HPP_ */
