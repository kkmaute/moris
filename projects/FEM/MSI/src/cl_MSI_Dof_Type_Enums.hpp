/*
 * cl_MSI_Dof_Type_Enums.hpp
 *
 *  Created on: Jul 23, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_
#define SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_

namespace moris
{
namespace MSI
{
    enum class Dof_Type
    {
        UX,     //< X-Displacement
        UY,     //< Y-Displacement
        UZ,     //< Z-Displacement
        TEMP,   //< Temperature degree of freedom
        L2,     //< Least Squares type
        MAPPING_DOF,
        LS1,    //< Level set
        LS2,    //< Level set
        NLSX,   //< X-Level set normal
        NLSY,   //< Y-Level set normal
        NLSZ,   //< Z-Level set normal
        VX,     //< X-Velocity
        VY,     //< Y-Velocity
        VZ,     //< Z-Velocity
        UNDEFINED, //< Undefined
        END_ENUM//
    };

    enum class Dv_Type
    {
        LS1,        //< Level set 1
        LS2,        //< Level set 2
        DENSITY,    //< Density
        UNDEFINED, //< Undefined
        END_ENUM//
    };
}
}

#endif /* SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_ */
