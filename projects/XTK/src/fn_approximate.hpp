/*
 * fn_approximate.hpp
 *
 *  Created on: Jan 11, 2018
 *      Author: doble
 */

#ifndef SRC_TOOLS_FN_APPROXIMATE_HPP_
#define SRC_TOOLS_FN_APPROXIMATE_HPP_

#include<limits>

namespace xtk
{
inline
bool approximate(size_t const & aNum1,
                 size_t const & aNum2)
{
    bool tEqual = false;

    if(aNum1 == aNum2)
    {
        tEqual = true;
    }

    return tEqual;
}

inline
bool approximate(sint const & aNum1,
                 sint const & aNum2)
{
    bool tEqual = false;

    if(aNum1 == aNum2)
    {
        tEqual = true;
    }

    return tEqual;
}

inline
bool approximate(real const & aNum1,
                 real const & aNum2,
                 real aErrorFactor = 1000.0)
{
    if(std::abs(aNum1-aNum2)<=std::numeric_limits<real>::epsilon()*aErrorFactor)
    {
        return true;
    }

    return false;
}

}


#endif /* SRC_TOOLS_FN_APPROXIMATE_HPP_ */
