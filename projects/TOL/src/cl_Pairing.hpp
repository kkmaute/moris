/*
 * cl_Pairing.hpp
 *
 *  Created on: Mar 29, 2017
 *      Author: doble
 */

#ifndef SRC_TOOLS_CL_PAIRING_HPP_
#define SRC_TOOLS_CL_PAIRING_HPP_

#include "cl_Mat.hpp" // LNA/src

namespace moris
{
    class Pairing
    {
    public:
        static moris::uint
        cantor_pairing(moris::Mat<moris::uint>  & aIDPair)
        {
            moris::uint tID = ((aIDPair(0,0)*aIDPair(0,0)+3*aIDPair(0,0)+2*aIDPair(0,0)*aIDPair(0,1)+aIDPair(0,1)+aIDPair(0,1)*aIDPair(0,1))/2);
            return tID;
        }

        static moris::uint
        cantor_pairing(moris::uint   aPt1,
                       moris::uint   aPt2)
        {
            moris::uint tID = ((aPt1*aPt1+3*aPt1+2*aPt1*aPt2+aPt2+aPt2*aPt2)/2);
            return tID;
        }
    };
}


#endif /* SRC_TOOLS_CL_PAIRING_HPP_ */
