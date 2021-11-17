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
    VIS,           //  Wrapper around XTK mesh database
    END_ENUM
};

enum class CellTopology
{
    LINE2,
    TRI3,
    TRI6,
    TRI10,
    QUAD4,
    QUAD9,
    QUAD16,
    TET4,
    TET10,
    HEX8,
    HEX27,
    HEX64,
    PRISM6,
    INVALID,
    END_ENUM
};

enum class CellShape
{
    // the following cell shapes increase specificity towards the top.
    // ie. Rectangular will also be parallel and straight. but a parallel cell shape will not be rectangular/
    RECTANGULAR, // rectangular and aligned with global csys quad or hex
    PARALLEL,    // Parallelogram
    STRAIGHT,    // all straight edges or planar faces faces
    GENERAL,
    EMPTY,
    INVALID,
    SIMPLEX,   // LINE, TRI, TET
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
        BSPLINE , // Indicates the entity has rank BSPLINE
        BSPLINE_2 , // Indicates the entity has rank BSPLINE
        BSPLINE_3 , // Indicates the entity has rank BSPLINE
        INVALID   , // Indicates the entity is invalid
        END_ENUM    //
    };

    enum class SetType
    {
        BULK                 ,   // Bulk Set
        SIDESET              ,   // SideSet
        DOUBLE_SIDED_SIDESET ,   // DoubleSided SideSet
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
        case EntityRank::BSPLINE  : return "BSPLINE"; break;
        case EntityRank::BSPLINE_2: return "BSPLINE_2"; break;
        case EntityRank::BSPLINE_3: return "BSPLINE_3"; break;
        case EntityRank::INVALID  : return "INVALID"; break;
        default: return "Unrecognized Enum provided to get_enum_str";
    }
}
inline
enum EntityRank get_entity_rank_from_str(std::string const & aString)
{
        if      (aString.compare("NODE"     ) == 0 || aString.compare("node"     ) == 0 ){ return EntityRank::NODE     ;}
        else if (aString.compare("EDGE"     ) == 0 || aString.compare("edge"     ) == 0 ){ return EntityRank::EDGE     ;}
        else if (aString.compare("FACE"     ) == 0 || aString.compare("face"     ) == 0 ){ return EntityRank::FACE     ;}
        else if (aString.compare("ELEMENT"  ) == 0 || aString.compare("element"  ) == 0 ){ return EntityRank::ELEMENT  ;}
        else if (aString.compare("BSPLINE"  ) == 0 || aString.compare("bspline"  ) == 0 ){ return EntityRank::BSPLINE  ;}
        else if (aString.compare("BSPLINE_2") == 0 || aString.compare("bspline_2") == 0 ){ return EntityRank::BSPLINE_2;}
        else if (aString.compare("BSPLINE_3") == 0 || aString.compare("bspline_3") == 0 ){ return EntityRank::BSPLINE_3;}
        else if (aString.compare("INVALID"  ) == 0 || aString.compare("invalid"  ) == 0 ){ return EntityRank::INVALID  ;}
        else{ MORIS_ERROR(0,"Invliad entity rank string"); return EntityRank::INVALID; };
}
inline
enum EntityRank get_entity_rank_from_index(moris_index aEntityRankIndex)
{
        if      (aEntityRankIndex == 0) { return EntityRank::NODE     ;}
        else if (aEntityRankIndex == 1) { return EntityRank::EDGE     ;}
        else if (aEntityRankIndex == 2) { return EntityRank::FACE     ;}
        else if (aEntityRankIndex == 3) { return EntityRank::ELEMENT  ;}
        else{ MORIS_ERROR(0,"Invalid index entity rank string"); return EntityRank::INVALID; };
}
}

#endif /* MORIS_MESH_CL_MESH_ENUMS_HPP_ */
