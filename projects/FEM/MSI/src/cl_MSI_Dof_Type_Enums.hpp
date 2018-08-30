/*
 * cl_Dof_Type_Enums.hpp
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
        END_ENUM//
    };
}
}

#endif /* SRC_FEM_CL_DOF_TYPE_ENUMS_HPP_ */
