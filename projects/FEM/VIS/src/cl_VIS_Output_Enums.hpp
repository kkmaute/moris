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
        UNDEFINED, //< Undefined
        END_ENUM//
    };

    enum class VIS_Mesh_Type
    {
        STANDARD,
        OVERLAPPING_INTERFACE,
        FULL_DISCONTINOUS,
        END_ENUM//
    };

    enum class Field_Type
    {
        NODAL,
        NODAL_IP,
        ELEMENTAL,
        GLOBAL,
        END_ENUM//
    };

}
}

#endif /* SRC_CL_VIS_ENUMS_HPP_ */
