#ifndef MORIS_CL_PDV_ENUMS_HPP_
#define MORIS_CL_PDV_ENUMS_HPP_

enum class PDV
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

moris::map< std::string, PDV > get_dv_type_map()
{
    moris::map< std::string, PDV > tPdvTypeMap;

    tPdvTypeMap["X_COORDINATE"]     = PDV::X_COORDINATE;
    tPdvTypeMap["Y_COORDINATE"]     = PDV::Y_COORDINATE;
    tPdvTypeMap["Z_COORDINATE"]     = PDV::Z_COORDINATE;
    tPdvTypeMap["DENSITY"]          = PDV::DENSITY;
    tPdvTypeMap["TEMPERATURE"]      = PDV::TEMPERATURE;
    tPdvTypeMap["ELASTIC_MODULUS"]  = PDV::ELASTIC_MODULUS;
    tPdvTypeMap["LS1"]              = PDV::LS1;
    tPdvTypeMap["LS2"]              = PDV::LS2;
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

#endif /* MORIS_CL_PDV_ENUMS_HPP_ */
