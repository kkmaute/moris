/*
 * cl_Geometry_Engine_Enums.hpp
 *
 *  Created on: Oct 4, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_ENUMS_HPP_
#define PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_ENUMS_HPP_



namespace ge
{
    enum class Root_Finding_Algorithm
    {
        // Not iterative
        Single_Linear_Interpolation,

        // Iterative Methods
        Bisection_Method,
        Newton_Method
    };
    enum class Field_Type
    {
        End_Enum
    };

}


#endif /* PROJECTS_GEN_SRC_NEW_GEOMENG_CL_GEN_ENUMS_HPP_ */
