/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Output_Enums.hpp
 *
 */

#ifndef SRC_CL_VIS_ENUMS_HPP_
#define SRC_CL_VIS_ENUMS_HPP_

#include "cl_Map.hpp"

namespace moris::vis
{
    enum class VIS_Mesh_Type
    {
        UNDEFINED,
        STANDARD,
        STANDARD_WITH_OVERLAP,
        FULL_DISCONTINUOUS,
        FULL_DISCONTINUOUS_WITH_OVERLAP,
        END_ENUM    //
    };

    enum class Field_Type
    {
        UNDEFINED,
        NODAL,
        NODAL_IP,
        ELEMENTAL_INT,
        ELEMENTAL_AVG,
        FACETED_INT,
        FACETED_AVG,
        GLOBAL,
        END_ENUM    //
    };

    inline moris::map< std::string, enum vis::Field_Type >
    get_vis_field_type_map()
    {
        moris::map< std::string, enum vis::Field_Type > tVisFieldTypeMap;

        tVisFieldTypeMap[ "UNDEFINED" ]     = vis::Field_Type::UNDEFINED;
        tVisFieldTypeMap[ "NODAL" ]         = vis::Field_Type::NODAL;
        tVisFieldTypeMap[ "NODAL_IP" ]      = vis::Field_Type::NODAL_IP;
        tVisFieldTypeMap[ "ELEMENTAL_INT" ] = vis::Field_Type::ELEMENTAL_INT;
        tVisFieldTypeMap[ "ELEMENTAL_AVG" ] = vis::Field_Type::ELEMENTAL_AVG;
        tVisFieldTypeMap[ "FACETED_INT" ]   = vis::Field_Type::FACETED_INT;
        tVisFieldTypeMap[ "FACETED_AVG" ]   = vis::Field_Type::FACETED_AVG;
        tVisFieldTypeMap[ "GLOBAL" ]        = vis::Field_Type::GLOBAL;
        tVisFieldTypeMap[ "UNDEFINED" ]     = vis::Field_Type::END_ENUM;

        return tVisFieldTypeMap;
    }
}    // namespace moris::vis

#endif /* SRC_CL_VIS_ENUMS_HPP_ */
