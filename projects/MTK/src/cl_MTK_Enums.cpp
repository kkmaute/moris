#include "cl_MTK_Enums.hpp"

       moris::map< std::string, moris::mtk::Master_Slave > 
       moris::mtk::get_master_type_map()
        {
            moris::map< std::string, moris::mtk::Master_Slave  > tMasterTypeMap;

            tMasterTypeMap["MASTER"]    = Master_Slave::MASTER;
            tMasterTypeMap["SLAVE"]     = Master_Slave::SLAVE;
            tMasterTypeMap["UNDEFINED"] = Master_Slave::UNDEFINED;
            return tMasterTypeMap;
        }

        moris::map< std::string, moris::mtk::Primary_Void > 
        moris::mtk::get_primary_type_map()
        {
            moris::map< std::string, moris::mtk::Primary_Void > tPrimaryTypeMap;

            tPrimaryTypeMap["PRIMARY"]   = moris::mtk::Primary_Void::PRIMARY;
            tPrimaryTypeMap["VOID"]      = moris::mtk::Primary_Void::VOID;
            tPrimaryTypeMap["INTERP"]    = moris::mtk::Primary_Void::INTERP;
            tPrimaryTypeMap["UNDEFINED"] = moris::mtk::Primary_Void::UNDEFINED;
            return tPrimaryTypeMap;
        }

        moris::map< std::string, moris::mtk::Field_Type > 
        moris::mtk::get_field_type_map()
        {
            moris::map< std::string, moris::mtk::Field_Type > tFieldTypeMap;

            tFieldTypeMap["FIELD_1"]   = moris::mtk::Field_Type::FIELD_1;
            tFieldTypeMap["FIELD_2"]   = moris::mtk::Field_Type::FIELD_2;
            tFieldTypeMap["FIELD_3"]   = moris::mtk::Field_Type::FIELD_3;
            tFieldTypeMap["FIELD_4"]   = moris::mtk::Field_Type::FIELD_4;
            tFieldTypeMap["FIELD_5"]   = moris::mtk::Field_Type::FIELD_5;
            tFieldTypeMap[""]          = moris::mtk::Field_Type::UNDEFINED;
            return tFieldTypeMap;
        }

        moris::map< std::string, moris::mtk::Field_Entity_Type >
        moris::mtk::get_field_entity_type_map()
        {
            moris::map< std::string, moris::mtk::Field_Entity_Type > tFieldEntityTypeMap;

            tFieldEntityTypeMap["NODAL"]       = moris::mtk::Field_Entity_Type::NODAL;
            tFieldEntityTypeMap["ELEMENTAL"]   = moris::mtk::Field_Entity_Type::ELEMENTAL;
            tFieldEntityTypeMap[""]            = moris::mtk::Field_Entity_Type::UNDEFINED;
            return tFieldEntityTypeMap;
        }

