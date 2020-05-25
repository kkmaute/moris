#ifndef MORIS_CL_PDV_Type_ENUMS_HPP_
#define MORIS_CL_PDV_Type_ENUMS_HPP_

#include "cl_Map.hpp"

namespace moris
{
    enum class PDV_Type
    {
        X_COORDINATE,
        Y_COORDINATE,
        Z_COORDINATE,
        DENSITY,
        TEMPERATURE,
        ELASTIC_MODULUS,
        LS1,
        LS2,
        UNDEFINED
    };

    moris::map< std::string, PDV_Type > get_pdv_type_map()
    {
        moris::map< std::string, PDV_Type > tPdvTypeMap;

        tPdvTypeMap["X_COORDINATE"]     = PDV_Type::X_COORDINATE;
        tPdvTypeMap["Y_COORDINATE"]     = PDV_Type::Y_COORDINATE;
        tPdvTypeMap["Z_COORDINATE"]     = PDV_Type::Z_COORDINATE;
        tPdvTypeMap["DENSITY"]          = PDV_Type::DENSITY;
        tPdvTypeMap["TEMPERATURE"]      = PDV_Type::TEMPERATURE;
        tPdvTypeMap["ELASTIC_MODULUS"]  = PDV_Type::ELASTIC_MODULUS;
        tPdvTypeMap["LS1"]              = PDV_Type::LS1;
        tPdvTypeMap["LS2"]              = PDV_Type::LS2;
        tPdvTypeMap[""]                 = PDV_Type::UNDEFINED;
        return tPdvTypeMap;
    }

    enum class Phase_Type
    {
        PHASE0,
        PHASE1,
        PHASE2,
        UNDEFINED,
        END_MAT_TYPE
    };
}


#endif /* MORIS_CL_PDV_ENUMS_HPP_ */
