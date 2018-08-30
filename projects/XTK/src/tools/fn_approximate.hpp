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
bool approximate(size_t const & aNum1,
                 size_t const & aNum2)
{
    bool tEqual = false;
    std::cout<<"HERE"<<std::endl;

    if(aNum1 == aNum2)
    {
        tEqual = true;
    }

    return tEqual;
}


bool approximate(sint const & aNum1,
                 sint const & aNum2)
{
    bool tEqual = false;
    std::cout<<"HERE"<<std::endl;

    if(aNum1 == aNum2)
    {
        tEqual = true;
    }

    return tEqual;
}

bool approximate(real const & aNum1,
                 real const & aNum2,
                 real aErrorFactor = 1000.0)
{
    std::cout<<"HERE"<<std::endl;
    if(std::abs(aNum1-aNum2)<=std::numeric_limits<real>::epsilon()*aErrorFactor)
    {
        return true;
    }

    return false;
}

}


#endif /* SRC_TOOLS_FN_APPROXIMATE_HPP_ */
