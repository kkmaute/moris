/*
 * cl_GE_Discrete.hpp
 *
 *  Created on: Mar 21, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_CL_GE_DISCRETE_HPP_
#define PROJECTS_GEN_SRC_CL_GE_DISCRETE_HPP_

// GE includes
#include "cl_GE_Geometry.hpp"


namespace moris
{
namespace ge
{
    class Discrete : public Geometry
    {
    private:

    protected:

    public:

        Discrete()
    {
        std::cout<<"Discrete constructor"<<std::endl;
    };
        ~Discrete(){};

    };
} /* namespace ge */
} /* namespace moris */



#endif /* PROJECTS_GEN_SRC_CL_GE_DISCRETE_HPP_ */
