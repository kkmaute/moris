/*
 * cl_VIS_Enums.hpp
 *
 *  Created on: Dec 10, 2019
 *      Author: schmidt
 */
#ifndef SRC_CL_VIS_ENUMS_HPP_
#define SRC_CL_VIS_ENUMS_HPP_

namespace moris
{
namespace vis
{
    enum class Output_Type
    {
        UNDEFINED, //< Undefined
        STRAIN_ENERGY,
        VOLUME,
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
        DOF, // Dof
        L2_ERROR_ANALYTIC, // L2 error for dof
        H1_ERROR_ANALYTIC, // H1 error for dof
        H1_SEMI_ERROR, // H1-semi error for dof
        END_ENUM//
    };

    enum class VIS_Mesh_Type
    {
        UNDEFINED,
        STANDARD,
        OVERLAPPING_INTERFACE,
        FULL_DISCONTINOUS,
        END_ENUM//
    };

    enum class Field_Type
    {
        UNDEFINED,
        NODAL,
        NODAL_IP,
        ELEMENTAL,
        GLOBAL,
        END_ENUM//
    };

}
}

#endif /* SRC_CL_VIS_ENUMS_HPP_ */
