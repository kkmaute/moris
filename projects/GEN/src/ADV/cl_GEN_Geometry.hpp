/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Geometry.hpp
 *
 */

#pragma once

#include <memory>
#include "cl_Matrix.hpp"
#include "cl_GEN_Node_Manager.hpp"

#include "cl_Cell.hpp" // TODO remove

// Forward declarations
namespace moris
{
    template< typename T > class Cell;
    namespace mtk
    {
        class Field;
        enum class Geometry_Type;
    }
}

namespace moris::ge
{
    // Forward declare intersection node classes
    class Intersection_Node;
    class Parent_Node;

    // Geometric location, for determining where a node is relative to a specific geometry
    enum class Geometric_Region : signed char
    {
        NEGATIVE = -1,
        INTERFACE = 0,
        POSITIVE = 1
    };

    class Geometry
    {
      public:

        /**
         * Constructor
         */
        Geometry();

        /**
         * Default destructor
         */
        virtual ~Geometry() = default;

        /**
         * Sets a new node manager (from the geometry engine, if it was created after this geometry).
         * Default implementation does nothing.
         *
         * @param aNodeManager Geometry engine node manager
         */
        virtual void set_node_manager( Node_Manager& aNodeManager );

        /**
         * Gets if this geometry depends on ADVs.
         *
         * @return ADV dependence
         */
        virtual bool depends_on_advs() = 0;

        /**
         * Gets the geometric region of a node, based on this geometry.
         *
         * @param aNodeIndex Node index
         * @param aNodeCoordinates Node coordinates
         * @return Geometric region enum
         */
        virtual Geometric_Region get_geometric_region(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aNodeCoordinates ) = 0;

        /**
         * Creates an intersection node based on the given information. The intersection node may or may not represent an intersection;
         * that is, its position may lie outside of the edge definition based on the given nodal coordinates. This information can be
         * requested from the created intersection node.
         *
         * @param aNodeIndex Node index of the new intersection node
         * @param aBaseNodes Base nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         * @return New intersection node
         */
        virtual Intersection_Node* create_intersection_node(
                uint                     aNodeIndex,
                const Cell< Node* >&     aBaseNodes,
                const Parent_Node&       aFirstParentNode,
                const Parent_Node&       aSecondParentNode,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder ) = 0;

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return Cell of MTK fields for remeshing
         */
        virtual Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() = 0;
    };
}
