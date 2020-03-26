/*
 * cl_GEN_Enums.hpp
 *
 *  Created on: Dec 19, 2019
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_SRC_ADDITIONAL_CL_GEN_ENUMS_HPP_
#define PROJECTS_GEN_SRC_ADDITIONAL_CL_GEN_ENUMS_HPP_

namespace moris
{
    namespace ge
    {
    //------------------------------------------------------------------------------
//        enum class GEN_PDV
//        {
//            XCOORD,
//            YCOORD,
//            ZCOORD,
//
//            RADIUS,
//            DENSITY0,
//            DENSITY1,
//            DENSITY2,
//            TEMPERATURE0,
//            TEMPERATURE1,
//            TEMPERATURE2,
//
//            END_ENUM
//        };
    //------------------------------------------------------------------------------
        enum class Root_Finding_Algorithm
        {
            // Not iterative
            Single_Linear_Interpolation,

            // Iterative Methods
            Bisection_Method,
            Newton_Method,

            END_ENUM
        };

    }   // end ge namespace
}       // end moris namespace
#endif /* PROJECTS_GEN_SRC_ADDITIONAL_CL_GEN_ENUMS_HPP_ */
