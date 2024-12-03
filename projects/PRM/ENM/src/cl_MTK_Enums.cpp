/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Enums.cpp
 *
 */

#include "cl_MTK_Enums.hpp"

namespace moris::mtk
{
    std::string get_enum_str( enum CellTopology aCellTopology )
    {
        switch ( aCellTopology )
        {
            case CellTopology::TRI3:
                return "TRI3";
            case CellTopology::QUAD4:
                return "QUAD4";
            case CellTopology::QUAD8:
                return "QUAD8";
            case CellTopology::TET4:
                return "TET4";
            case CellTopology::TET10:
                return "TET10";
            case CellTopology::HEX8:
                return "HEX8";
            case CellTopology::HEX20:
                return "HEX20";
            case CellTopology::HEX27:
                return "HEX27";
            case CellTopology::HEX64:
                return "HEX64";
            case CellTopology::PRISM6:
                return "PRISM6";
            case CellTopology::UNDEFINED:
                return "UNDEFINED";
            default:
                return "invalid cell topology enum provided";
        }
    }

    //------------------------------------------------------------------------------

    uint get_order_from_topology( CellTopology aCellTopology )
    {
        switch ( aCellTopology )
        {
            case CellTopology::TRI3:
            case CellTopology::QUAD4:
            case CellTopology::TET4:
            case CellTopology::PRISM6:
            case CellTopology::HEX8:
                return 1;
            case CellTopology::QUAD8:
            case CellTopology::TET10:
            case CellTopology::HEX20:
            case CellTopology::HEX27:
                return 2;
            case CellTopology::HEX64:
                return 3;
            case CellTopology::UNDEFINED:
            default:
                return 0;
        }
    }

    std::string get_enum_str( EntityRank aCellTopology )
    {
        switch ( aCellTopology )
        {
            case EntityRank::NODE:
                return "NODE";
            case EntityRank::EDGE:
                return "EDGE";
            case EntityRank::FACE:
                return "FACE";
            case EntityRank::ELEMENT:
                return "ELEMENT";
            case EntityRank::BSPLINE:
                return "BSPLINE";
            case EntityRank::BSPLINE_2:
                return "BSPLINE_2";
            case EntityRank::BSPLINE_3:
                return "BSPLINE_3";
            case EntityRank::UNDEFINED:
            default:
                return "UNDEFINED";
        }
    }

    //------------------------------------------------------------------------------

    EntityRank get_entity_rank_from_str( const std::string& aString )
    {
        if ( aString.compare( "NODE" ) == 0 || aString.compare( "node" ) == 0 ) { return EntityRank::NODE; }
        else if ( aString.compare( "EDGE" ) == 0 || aString.compare( "edge" ) == 0 )
        {
            return EntityRank::EDGE;
        }
        else if ( aString.compare( "FACE" ) == 0 || aString.compare( "face" ) == 0 )
        {
            return EntityRank::FACE;
        }
        else if ( aString.compare( "ELEMENT" ) == 0 || aString.compare( "element" ) == 0 )
        {
            return EntityRank::ELEMENT;
        }
        else if ( aString.compare( "BSPLINE" ) == 0 || aString.compare( "bspline" ) == 0 )
        {
            return EntityRank::BSPLINE;
        }
        else if ( aString.compare( "BSPLINE_2" ) == 0 || aString.compare( "bspline_2" ) == 0 )
        {
            return EntityRank::BSPLINE_2;
        }
        else if ( aString.compare( "BSPLINE_3" ) == 0 || aString.compare( "bspline_3" ) == 0 )
        {
            return EntityRank::BSPLINE_3;
        }
        else if ( aString.compare( "UNDEFINED" ) == 0 || aString.compare( "invalid" ) == 0 )
        {
            return EntityRank::INVALID;
        }
        else
        {
            MORIS_ERROR( false, "Invliad entity rank string" );
            return EntityRank::INVALID;
        };
    }

    //------------------------------------------------------------------------------

    EntityRank get_entity_rank_from_index( moris_index aEntityRankIndex )
    {
        if ( aEntityRankIndex == 0 ) { return EntityRank::NODE; }
        else if ( aEntityRankIndex == 1 )
        {
            return EntityRank::EDGE;
        }
        else if ( aEntityRankIndex == 2 )
        {
            return EntityRank::FACE;
        }
        else if ( aEntityRankIndex == 3 )
        {
            return EntityRank::ELEMENT;
        }
        else
        {
            MORIS_ERROR( false, "Invalid index entity rank string" );
            return EntityRank::INVALID;
        };
    }

    //------------------------------------------------------------------------------

    map< std::string, Leader_Follower >
    get_leader_type_map()
    {
        map< std::string, Leader_Follower > tLeaderTypeMap;

        tLeaderTypeMap[ "LEADER" ]    = Leader_Follower::LEADER;
        tLeaderTypeMap[ "FOLLOWER" ]  = Leader_Follower::FOLLOWER;
        tLeaderTypeMap[ "UNDEFINED" ] = Leader_Follower::UNDEFINED;
        return tLeaderTypeMap;
    }

    //------------------------------------------------------------------------------

    Vector< Leader_Follower > get_leader_follower( SetType aSetType )
    {
        if ( aSetType == SetType::DOUBLE_SIDED_SIDESET )
        {
            return { Leader_Follower::LEADER, Leader_Follower::FOLLOWER };
        }
        else
        {
            return { Leader_Follower::LEADER };
        }
    }

    //------------------------------------------------------------------------------

    map< std::string, Primary_Void >
    get_primary_type_map()
    {
        map< std::string, Primary_Void > tPrimaryTypeMap;

        tPrimaryTypeMap[ "PRIMARY" ]   = Primary_Void::PRIMARY;
        tPrimaryTypeMap[ "VOID" ]      = Primary_Void::VOID;
        tPrimaryTypeMap[ "INTERP" ]    = Primary_Void::INTERP;
        tPrimaryTypeMap[ "UNDEFINED" ] = Primary_Void::UNDEFINED;
        return tPrimaryTypeMap;
    }

    //------------------------------------------------------------------------------

    map< std::string, Field_Type >
    get_field_type_map()
    {
        map< std::string, Field_Type > tFieldTypeMap;

        tFieldTypeMap[ "FIELD_0" ] = Field_Type::FIELD_0;
        tFieldTypeMap[ "FIELD_1" ] = Field_Type::FIELD_1;
        tFieldTypeMap[ "FIELD_2" ] = Field_Type::FIELD_2;
        tFieldTypeMap[ "FIELD_3" ] = Field_Type::FIELD_3;
        tFieldTypeMap[ "FIELD_4" ] = Field_Type::FIELD_4;
        tFieldTypeMap[ "FIELD_5" ] = Field_Type::FIELD_5;
        tFieldTypeMap[ "FIELD_6" ] = Field_Type::FIELD_6;
        tFieldTypeMap[ "FIELD_7" ] = Field_Type::FIELD_7;
        tFieldTypeMap[ "FIELD_8" ] = Field_Type::FIELD_8;
        tFieldTypeMap[ "FIELD_9" ] = Field_Type::FIELD_9;
        tFieldTypeMap[ "" ]        = Field_Type::UNDEFINED;
        return tFieldTypeMap;
    }

    //------------------------------------------------------------------------------

    map< std::string, Field_Entity_Type >
    get_field_entity_type_map()
    {
        map< std::string, Field_Entity_Type > tFieldEntityTypeMap;

        tFieldEntityTypeMap[ "NODAL" ]     = Field_Entity_Type::NODAL;
        tFieldEntityTypeMap[ "ELEMENTAL" ] = Field_Entity_Type::ELEMENTAL;
        tFieldEntityTypeMap[ "" ]          = Field_Entity_Type::UNDEFINED;
        return tFieldEntityTypeMap;
    }

    std::string get_leader_follower_string( Leader_Follower aLeaderFollower )
    {
        switch ( aLeaderFollower )
        {
            case Leader_Follower::LEADER:
                return "leader";
            case Leader_Follower::FOLLOWER:
                return "follower";
            case Leader_Follower::UNDEFINED:
                return "undefined";
            default:
                return "invalid leader follower enum provided";
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk
