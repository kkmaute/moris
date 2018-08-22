/*
 * cl_Geometry_Engine_Enums.hpp
 *
 *  Created on: Oct 4, 2017
 *      Author: ktdoble
 */

#ifndef SRC_GEOMENG_CL_MGE_ENUMS_HPP_
#define SRC_GEOMENG_CL_MGE_ENUMS_HPP_



namespace xtk
{
    enum class Root_Finding_Algorithm
    {
        // Not iterative
        Single_Linear_Interpolation,

        // Iterative Methods
        Bisection_Method,
        Newton_Method
    };

}


#endif /* SRC_GEOMENG_CL_MGE_ENUMS_HPP_ */
