/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Enums.hpp
 *
 */

#ifndef SRC_MESH_CL_MTK_ENUMS_HPP_
#define SRC_MESH_CL_MTK_ENUMS_HPP_

#include "assert.hpp"
#include "moris_typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Vector.hpp"

namespace moris::mtk
{

    enum class MeshType
    {
        MTK,    //< Wrapper around STK mesh database
        STK,    //< Wrapper around STK mesh database
        HMR,    //  Wrapper around HMR mesh database
        XTK,    //  Wrapper around XTK mesh database
        VIS,    //  Wrapper around XTK mesh database
        UNDEFINED
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
        UNDEFINED
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
        UNDEFINED
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
        UNDEFINED     //
    };

    enum class SetType
    {
        BULK,                    // Bulk Set
        SIDESET,                 // SideSet
        DOUBLE_SIDED_SIDESET,    // DoubleSided SideSet
        UNDEFINED                //
    };

    enum class ClusterType
    {
        CELL,            // Cell Cluster
        SIDE,            // Side Cluster
        DOUBLE_SIDED,    // DoubleSided Cluster
        UNDEFINED
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
     * @return EntityRank enum
     */
    EntityRank get_entity_rank_from_index( moris_index aEntityRankIndex );

    //------------------------------------------------------------------------------

    enum class Interpolation_Type
    {
        CONSTANT, // constant interpolation
        LAGRANGE, // the most common finite element types
        BEZIER,   // Bezier type elements
        UNDEFINED
    };

    //------------------------------------------------------------------------------

    enum class Integration_Type
    {
        CONSTANT,
        GAUSS, // Gauss ( Quad and Hex ), Dunavant ( Tri ), Hammer ( Tet )
        UNDEFINED
    };

    //------------------------------------------------------------------------------

    enum class Integration_Order
    {
        POINT,
        BAR_1,
        BAR_2,
        BAR_3,
        BAR_4,
        BAR_5,
        BAR_6,
        QUAD_1x1,
        QUAD_2x2,
        QUAD_3x3,
        QUAD_4x4,
        QUAD_5x5,
        TRI_1,
        TRI_3,
        TRI_4,
        TRI_6,
        TRI_7,
        TRI_12,
        TRI_13,
        TRI_16,
        TRI_19,
        TRI_25,
        HEX_1x1x1,
        HEX_2x2x2,
        HEX_3x3x3,
        HEX_4x4x4,
        HEX_5x5x5,
        TET_1,
        TET_4,
        TET_5,
        TET_10,
        TET_11,
        TET_15,
        TET_20,
        TET_35,
        TET_56,
        UNDEFINED
    };

    //------------------------------------------------------------------------------

    enum class Geometry_Type
    {
            POINT, // point
            LINE,  // 1D line or curve
            QUAD,  // rectangle
            TRI,   // triangle
            HEX,   // quadrangle
            TET,   // tetrahedron
            PENTA, // pentahedron
            UNDEFINED
    };

    //------------------------------------------------------------------------------

    enum class Interpolation_Order
    {
            CONSTANT,
            LINEAR,
            QUADRATIC,
            SERENDIPITY,
            CUBIC,
            UNDEFINED
    };

    //------------------------------------------------------------------------------

    enum class Leader_Follower
    {
            LEADER,
            FOLLOWER,
            UNDEFINED
    };

    /**
     * Gets the map that can be used to go from a parameter list value to the internal enum Leader_Follower
     *
     * @return map from std::string to mtk::Leader_Follower
     */
    moris::map< std::string, mtk::Leader_Follower > get_leader_type_map();

    /**
     * Gets a cell of leader/follower enums to loop over depending on the type of a given set
     *
     * @param aSetType Set enum
     * @return Leader/follower cell, either { Leader } or { Leader, Follower }
     */
    Vector< Leader_Follower > get_leader_follower( SetType aSetType );

    //------------------------------------------------------------------------------

    enum class Primary_Void
    {
            PRIMARY,
            VOID,
            INTERP,
            UNDEFINED
    };

    /**
     * Gets the map that can be used to go from a parameter list value to the internal enum Primary_Void
     *
     * @return map from std::string to mtk::Primary_Void enum
     */
    moris::map< std::string, mtk::Primary_Void > get_primary_type_map();

    //------------------------------------------------------------------------------

    enum class Field_Type
    {
            FIELD_1,
            FIELD_2,
            FIELD_3,
            FIELD_4,
            FIELD_5,
            UNDEFINED
    };

    /**
     * Gets the map that can be used to go from a parameter list value to the internal enum Field_Type
     *
     * @return map from std::string to mtk::Field_Type
     */
    moris::map< std::string, Field_Type > get_field_type_map();

    //------------------------------------------------------------------------------

    enum class Field_Entity_Type
    {
            NODAL,
            ELEMENTAL,
            UNDEFINED
    };

    enum class Field_Implementation
    {
            MTK,
            FEM,
            GEN,
            UNDEFINED
    };

    /**
     * Gets the map that can be used to go from a parameter list value to the internal enum Field_Entity_Type
     *
     * @return map from std::string to mtk::Field_Entity_Type
     */
    moris::map< std::string, Field_Entity_Type > get_field_entity_type_map();

    //------------------------------------------------------------------------------
}

#endif /* SRC_MESH_CL_MTK_ENUMS_HPP_ */

