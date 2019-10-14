/*
 * cl_Mesh_Enums.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: doble
 */

#ifndef MORIS_MESH_CL_MESH_ENUMS_HPP_
#define MORIS_MESH_CL_MESH_ENUMS_HPP_

#include "assert.hpp"

//namespace moris
//{

enum class MeshType
{

    MTK,           //< Wrapper around STK mesh database
    STK,           //< Wrapper around STK mesh database
    HMR,           //  Wrapper around HMR mesh database
    XTK,           //  Wrapper around XTK mesh database
    END_ENUM
};

enum class CellTopology
{
    LINE2,
    TRI3,
    QUAD4,
    QUAD9,
    TET4,
    TET10,
    HEX8,
    HEX27,
    HEX64,
    PRISM6,
    INVALID,
    END_ENUM
};

inline
const std::string get_enum_str(enum CellTopology aCellTopoEnum)
{
    switch (aCellTopoEnum)
    {

        case CellTopology::TRI3: return "TRI3";
        case CellTopology::QUAD4: return "QUAD4";
        case CellTopology::TET4:  return "TET4";
        case CellTopology::TET10: return "TET10";
        case CellTopology::HEX8: return "HEX8";
        case CellTopology::HEX27: return "HEX27";
        case CellTopology::HEX64: return "HEX64";
        case CellTopology::PRISM6: return "PRISM6";
        case CellTopology::INVALID: return "INVALID";
        case CellTopology::END_ENUM: return "END_ENUM";
       default: return "invalid cell topology enum provided";
    }
}



namespace moris
{
    enum class EntityRank
    {
        NODE      ,   // Indicates the entity has rank NODE
        EDGE      ,   // Indicates the entity has rank EDGE
        FACE      ,   // Indicates the entity has rank FACE
        ELEMENT   ,// Indicates the entity has rank ELEMENT
        BSPLINE_1 , // Indicates the entity has rank BSPLINE
        BSPLINE_2 , // Indicates the entity has rank BSPLINE
        BSPLINE_3 , // Indicates the entity has rank BSPLINE
        INVALID   , // Indicates the entity is invalid
        END_ENUM    //
    };

inline
const std::string get_enum_str(enum EntityRank aCellTopoEnum)
{
    switch (aCellTopoEnum)
    {

        case EntityRank::NODE     : return "NODE"; break;
        case EntityRank::EDGE     : return "EDGE"; break;
        case EntityRank::FACE     : return "FACE"; break;
        case EntityRank::ELEMENT  : return "ELEMENT"; break;
        case EntityRank::BSPLINE_1: return "BSPLINE_1"; break;
        case EntityRank::BSPLINE_2: return "BSPLINE_2"; break;
        case EntityRank::BSPLINE_3: return "BSPLINE_3"; break;
        case EntityRank::INVALID  : return "INVALID"; break;
        default: return "Unrecognized Enum provided to get_enum_str";
    }
}
}
#endif /* MORIS_MESH_CL_MESH_ENUMS_HPP_ */
