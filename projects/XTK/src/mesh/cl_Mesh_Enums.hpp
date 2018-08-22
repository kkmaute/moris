/*
 * cl_Mesh_Enums.hpp
 *
 *  Created on: Jun 19, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_ENUMS_HPP_
#define SRC_MESH_CL_MESH_ENUMS_HPP_

enum class EntityRank
{
    NODE, EDGE, FACE, ELEMENT, END_ENUM, INVALID_RANK
};

enum class EntityTopology
{
    TET_4,
    TET_10,
    HEXA_8,
    INVALID
};


#endif /* SRC_MESH_CL_MESH_ENUMS_HPP_ */
