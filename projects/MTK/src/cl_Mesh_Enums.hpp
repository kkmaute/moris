/*
 * cl_tools_Enums.hpp
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
    END_ENUM
};

enum class CellTopology
{
    TET4,
    HEX8,
    END_ENUM
};


namespace moris
{
    enum class EntityRank
    {
        NODE,   // Indicates the entity has rank NODE
        EDGE,   // Indicates the entity has rank EDGE
        FACE,   // Indicates the entity has rank FACE
        ELEMENT,// Indicates the entity has rank ELEMENT
        BSPLINE_1, // Indicates the entity has rank BSPLINE
        BSPLINE_2, // Indicates the entity has rank BSPLINE
        BSPLINE_3, // Indicates the entity has rank BSPLINE
        INVALID, // Indicates the entity is invalid
        END_ENUM//
    };


    uint
    entity_rank_to_order( const EntityRank aEntityRank )
    {
        switch( aEntityRank )
        {
            case( EntityRank::BSPLINE_1 ) :
            {
                return 1;
                break;
            }
            case( EntityRank::BSPLINE_2 ) :
            {
                return 2;
                break;
            }
            case( EntityRank::BSPLINE_3 ) :
            {
                return 3;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "entity_rank_to_order() : invalid input" );
                return 0;
                break;
            }
        }
    }

}
//enum DerivativeOrder
//{
//    ZEROTH_ORDER,           // zeroth order derivative
//    FIRST_ORDER,        // first order derivative
//    SECOND_ORDER        // second order derivative
//};
//}   // namespace moris
#endif /* MORIS_MESH_CL_MESH_ENUMS_HPP_ */
