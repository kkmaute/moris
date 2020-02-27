/*
 * cl_GEN_Dv_Enums.hpp
 *
 *  Created on: Jan 15, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_GEN_GEN_CORE_SRC_CL_GEN_DV_ENUMS_HPP_
#define PROJECTS_GEN_GEN_CORE_SRC_CL_GEN_DV_ENUMS_HPP_

enum class GEN_DV
{
    UNDEFINED,
    XCOORD,
    YCOORD,
    ZCOORD,
    DENSITY0,
    DENSITY1,
    DENSITY2,
    LS1,
    LS2,
    LS3,
    RADIUS,
    TEMPERATURE0,
    TEMPERATURE1,
    TEMPERATURE2,
    MODULUS0,
    MODULUS1,
    MODULUS2,
    END_ENUM
};

moris::map< std::string, enum GEN_DV > get_dv_type_map()
{
    moris::map< std::string, enum GEN_DV > tDvTypeMap;

    tDvTypeMap["XCOORD"]       = GEN_DV::XCOORD;
    tDvTypeMap["YCOORD"]       = GEN_DV::YCOORD;
    tDvTypeMap["ZCOORD"]       = GEN_DV::ZCOORD;
    tDvTypeMap["DENSITY0"]     = GEN_DV::DENSITY0;
    tDvTypeMap["DENSITY1"]     = GEN_DV::DENSITY1;
    tDvTypeMap["DENSITY2"]     = GEN_DV::DENSITY2;
    tDvTypeMap["LS1"]          = GEN_DV::LS1;
    tDvTypeMap["LS2"]          = GEN_DV::LS2;
    tDvTypeMap["LS3"]          = GEN_DV::LS3;
    tDvTypeMap["RADIUS"]       = GEN_DV::RADIUS;
    tDvTypeMap["TEMPERATURE0"] = GEN_DV::TEMPERATURE0;
    tDvTypeMap["TEMPERATURE1"] = GEN_DV::TEMPERATURE1;
    tDvTypeMap["TEMPERATURE2"] = GEN_DV::TEMPERATURE2;
    tDvTypeMap["MODULUS0"]     = GEN_DV::MODULUS0;
    tDvTypeMap["MODULUS1"]     = GEN_DV::MODULUS1;
    tDvTypeMap["MODULUS2"]     = GEN_DV::MODULUS2;
    return tDvTypeMap;
}

//------------------------------------------------------------------------------
enum class Mat_Type
{
    UNDEFINED,
    PHASE0,
    PHASE1,
    PHASE2,
    END_MAT_TYPE
};

#endif /* PROJECTS_GEN_GEN_CORE_SRC_CL_GEN_DV_ENUMS_HPP_ */
