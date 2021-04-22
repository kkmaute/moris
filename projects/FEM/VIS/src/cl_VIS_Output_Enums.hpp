/*
 * cl_VIS_Enums.hpp
 *
 *  Created on: Dec 10, 2019
 *      Author: schmidt
 */
#ifndef SRC_CL_VIS_ENUMS_HPP_
#define SRC_CL_VIS_ENUMS_HPP_

#include "cl_Map.hpp"

namespace moris
{
namespace vis
{
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

    inline
    moris::map< std::string, enum vis::Field_Type > get_vis_field_type_map()
    {
        moris::map< std::string, enum vis::Field_Type > tVisFieldTypeMap;

        tVisFieldTypeMap["UNDEFINED"] = vis::Field_Type::UNDEFINED;
        tVisFieldTypeMap["NODAL"]     = vis::Field_Type::NODAL;
        tVisFieldTypeMap["NODAL_IP"]  = vis::Field_Type::NODAL_IP;
        tVisFieldTypeMap["ELEMENTAL"] = vis::Field_Type::ELEMENTAL;
        tVisFieldTypeMap["GLOBAL"]    = vis::Field_Type::GLOBAL;
        tVisFieldTypeMap["END_ENUM"]  = vis::Field_Type::END_ENUM;

        return tVisFieldTypeMap;
    }
}
}

#endif /* SRC_CL_VIS_ENUMS_HPP_ */
