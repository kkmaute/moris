/*
 * cl_Pairing.hpp
 *
 *  Created on: Jun 26, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOOLS_FN_PAIRING_HPP_
#define SRC_TOOLS_FN_PAIRING_HPP_

#include "linalg/cl_XTK_Matrix_Base.hpp"

namespace xtk
{

template<typename Integer>
static Integer
cantor_pairing(Integer  const aPt1,
               Integer  const aPt2)
{
    Integer tID = ((aPt1*aPt1+3*aPt1+2*aPt1*aPt2+aPt2+aPt2*aPt2)/2);
    return tID;
}
}


#endif /* SRC_TOOLS_FN_PAIRING_HPP_ */
