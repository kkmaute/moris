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
#include "cl_GEN_Design.hpp"
#include "cl_MTK_Mesh_Pair.hpp"

#include "cl_Cell.hpp"    // TODO remove

// Forward declarations
namespace moris
{
    template< typename T >
    class Cell;
    namespace mtk
    {
        class Field;
        enum class Geometry_Type;
    }    // namespace mtk
}    // namespace moris

namespace moris::ge
{
    // Forward declare intersection node classes
    class Intersection_Node;
    class Parent_Node;

    // Geometric location, for determining where a node is relative to a specific geometry
    enum class Geometric_Region : signed char
    {
        NEGATIVE  = -1,
        INTERFACE = 0,
        POSITIVE  = 1
    };

    class Geometry : public Design
    {
      private:
        real mIntersectionTolerance;

      public:
        /**
         * Constructor
         */
        explicit Geometry( Design_Parameters aParameters, real aIntersectionTolerance );

        /**
         * Default destructor
         */
        ~Geometry() override = default;

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
        virtual bool depends_on_advs() const = 0;

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
         * @param aBackgroundNodes Background nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         * @return New intersection node
         */
        virtual Intersection_Node* create_intersection_node(
                uint                     aNodeIndex,
                const Cell< Node* >&     aBackgroundNodes,
                const Parent_Node&       aFirstParentNode,
                const Parent_Node&       aSecondParentNode,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder ) = 0;

        /**
         * Computes the local coordinate along a parent edge of an intersection node created using this geometry.
         *
         * @param aBackgroundNodes Background nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @return Parent edge local coordinate, between -1 and 1
         */
        virtual real compute_intersection_local_coordinate(
                const Cell< Node* >& aBackgroundNodes,
                const Parent_Node&   aFirstParentNode,
                const Parent_Node&   aSecondParentNode ) = 0;

        /**
         * Gets an MTK field, if this geometry uses one that needs to be remapped to a new mesh
         *
         * @return Cell of MTK fields for remeshing
         */
        virtual Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() = 0;

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        virtual void import_advs( sol::Dist_Vector* aOwnedADVs ) = 0;

        /**
         * Resets all nodal information, including child nodes. This should be called when a new XTK mesh is being
         * created.
         *
         * @param aInterpolationMesh Interpolation mesh containing new nodal data
         */
        virtual void reset_nodal_data( mtk::Interpolation_Mesh* aInterpolationMesh ) = 0;

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMeshPair The mesh pair where the discretization information can be obtained
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         */
        virtual void discretize(
                mtk::Mesh_Pair          aMeshPair,
                sol::Dist_Vector*       aOwnedADVs,
                const Matrix< DDSMat >& aSharedADVIds,
                uint                    aADVOffsetID ) = 0;

        /**
         * If intended for this field, maps the field to B-spline coefficients or stores the nodal field values in a stored field object.
         *
         * @param aMTKField Input MTK field to map based on
         * @param aOwnedADVs Pointer to the owned distributed ADVs
         * @param aSharedADVIds All owned and shared ADV IDs for this B-spline field
         * @param aADVOffsetID Offset in the owned ADV IDs for pulling ADV IDs
         */
        virtual void discretize(
                std::shared_ptr< mtk::Field > aMTKField,
                mtk::Mesh_Pair                aMeshPair,
                sol::Dist_Vector*             aOwnedADVs,
                const Matrix< DDSMat >&       aSharedADVIds,
                uint                          aADVOffsetID ) = 0;

        /**
         * Used to print geometry information to exodus files and print debug information.
         *
         * @param aNodeIndex decides the point at which the field value is printed. If the node is a derived node, the value is interpolated from the parents.
         * @param aCoordinates The field location to get the value from.
         * @return the value of the property field at the requested location
         */
        virtual void get_design_info(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Cell< real >&           aOutputDesignInfo ) = 0;

        /**
         * Gets the intersection tolerance for this geometry
         *
         * @return tolerance for creating intersection nodes
         */
        virtual real get_intersection_tolerance()
        {
            return mIntersectionTolerance;
        }
    };
}    // namespace moris::ge
