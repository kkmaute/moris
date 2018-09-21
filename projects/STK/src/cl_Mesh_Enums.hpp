/*
 * cl_tools_Enums.hpp
 *
 *  Created on: Oct 20, 2016
 *      Author: doble
 */

#ifndef MORIS_MESH_CL_MESH_ENUMS_HPP_
#define MORIS_MESH_CL_MESH_ENUMS_HPP_
//namespace moris
//{

enum class MeshType
{

    MTK,           //< Wrapper around STK mesh database
    HMR,           //  Wrapper around HMR mesh database
    FTK            //< Wrapper around FutureToolKit mesh database

};

enum class EntityRank
{
    NODE,   // Indicates the entity has rank NODE
    EDGE,   // Indicates the entity has rank EDGE
    FACE,   // Indicates the entity has rank FACE
    ELEMENT,// Indicates the entity has rank ELEMENT
    INVALID,// Indicates the entity is invalid
    END_ENUM//
};

//enum DerivativeOrder
//{
//    ZEROTH_ORDER,           // zeroth order derivative
//    FIRST_ORDER,        // first order derivative
//    SECOND_ORDER        // second order derivative
//};
//}   // namespace moris
#endif /* MORIS_MESH_CL_MESH_ENUMS_HPP_ */
