/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Mesh_Enums.hpp
 *
 */

#ifndef MORIS_MESH_CL_MESH_ENUMS_HPP_
#define MORIS_MESH_CL_MESH_ENUMS_HPP_

#include "assert.hpp"

// namespace moris
//{

enum class MeshType
{
    MTK,    //< Wrapper around STK mesh database
    STK,    //< Wrapper around STK mesh database
    HMR,    //  Wrapper around HMR mesh database
    XTK,    //  Wrapper around XTK mesh database
    VIS,    //  Wrapper around XTK mesh database
    END_ENUM
};

enum class CellTopology
{
    LINE2,
    TRI3,
    TRI6,
    TRI10,
    QUAD4,
    QUAD8,
    QUAD9,
    QUAD16,
    TET4,
    TET10,
    HEX8,
    HEX20,
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
    RECTANGULAR,    // rectangular and aligned with global csys quad or hex
    PARALLEL,       // Parallelogram
    STRAIGHT,       // all straight edges or planar faces faces
    GENERAL,
    EMPTY,
    INVALID,
    SIMPLEX,    // LINE, TRI, TET
    END_ENUM
};

/**
 * Gets a string representing the given CellTopology enum
 *
 * @param aCellTopology CellTopology enum
 * @return String with enum name
 */
std::string get_enum_str( CellTopology ACellTopology );

/**
 * Gets just the order from a cell topology
 *
 * @param aCellTopology CellTopology enum
 * @return Polynomial order
 */
uint get_order_from_topology( CellTopology aCellTopology );

namespace moris
{
    enum class EntityRank
    {
        NODE,         // Indicates the entity has rank NODE
        EDGE,         // Indicates the entity has rank EDGE
        FACE,         // Indicates the entity has rank FACE
        ELEMENT,      // Indicates the entity has rank ELEMENT
        BSPLINE,      // Indicates the entity has rank BSPLINE
        BSPLINE_2,    // Indicates the entity has rank BSPLINE
        BSPLINE_3,    // Indicates the entity has rank BSPLINE
        INVALID,      // Indicates the entity is invalid
        END_ENUM      //
    };

    enum class SetType
    {
        BULK,                    // Bulk Set
        SIDESET,                 // SideSet
        DOUBLE_SIDED_SIDESET,    // DoubleSided SideSet
        END_ENUM                 //
    };

    enum class ClusterType
    {
        CELL_CLUSTER,            // Cell Cluster
        SIDE_CLUSTER,            // Side Cluster
        DOUBLE_SIDED_CLUSTER,    // DoubleSided Cluster
        END_ENUM
    };

    /**
     * Gets a string representing the given EntityRank enum
     *
     * @param aEntityRank EntityRank enum
     * @return String with enum name
     */
    std::string get_enum_str( EntityRank aEntityRank );

    /**
     * Gets an EntityRank enum from a given string
     *
     * @param aString String with enum name
     * @return EntityRank enum
     */
    EntityRank get_entity_rank_from_str( std::string const & aString );

    /**
     * Gets the entity rank from a given index
     * @param aEntityRankIndex
     * @return
     */
    EntityRank get_entity_rank_from_index( moris_index aEntityRankIndex );
}

#endif /* MORIS_MESH_CL_MESH_ENUMS_HPP_ */

