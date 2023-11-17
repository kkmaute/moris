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

// Forward declarations
namespace moris
{
    template< typename T > class Cell;
    namespace mtk
    {
        class Field;
    }
}

namespace moris::ge
{
    // Forward declare intersection node
    class Intersection_Node;

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
         * @param aEdgeFirstNodeIndex First node index on the intersection edge
         * @param aEdgeSecondNodeIndex Second node index on the intersection edge
         * @param aEdgeFirstIntersectionNode First intersection node on the intersection edge, if it is also an intersection
         * @param aEdgeSecondIntersectionNode Second intersection node on the intersection edge, if it is also an intersection
         * @param aEdgeFirstNodeLocalCoordinates Local coordinates of the first node inside the background element
         * @param aEdgeSecondNodeLocalCoordinates Local coordinates of the second node inside the background element
         * @param aEdgeFirstNodeGlobalCoordinates Global coordinates of the first node
         * @param aEdgeSecondNodeGlobalCoordinates Global coordinates of the second node
         * @param aBackgroundElementNodeIndices Node indices of the background element
         * @param aBackgroundElementNodeCoordinates Node coordinates of the background element
         * @return Created intersection node
         */
        virtual std::shared_ptr< Intersection_Node > create_intersection_node(
                uint                                 aEdgeFirstNodeIndex,
                uint                                 aEdgeSecondNodeIndex,
                std::shared_ptr< Intersection_Node > aEdgeFirstIntersectionNode,
                std::shared_ptr< Intersection_Node > aEdgeSecondIntersectionNode,
                const Matrix< DDRMat >&              aEdgeFirstNodeLocalCoordinates,
                const Matrix< DDRMat >&              aEdgeSecondNodeLocalCoordinates,
                const Matrix< DDRMat >&              aEdgeFirstNodeGlobalCoordinates,
                const Matrix< DDRMat >&              aEdgeSecondNodeGlobalCoordinates,
                const Matrix< DDUMat >&              aBackgroundElementNodeIndices,
                const Cell< Matrix< DDRMat > >&      aBackgroundElementNodeCoordinates ) = 0;

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return Cell of MTK fields for remeshing
         */
        virtual Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() = 0;
    };
}
