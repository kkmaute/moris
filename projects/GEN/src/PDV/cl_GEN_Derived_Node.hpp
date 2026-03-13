/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Derived_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Interpolation_Function.hpp"

namespace moris::mtk
{
    enum class Geometry_Type;
    enum class Interpolation_Order;
}    // namespace moris::mtk

namespace moris::gen
{
    // Forward declare basis node
    class Background_Node;
    class Basis_Node;
    class Geometry;
    class Node_Manager;

    class Derived_Node : public Node
    {
      private:
        // Vector< Basis_Node > mBackgroundNodes; // brendan delete
        const mtk::Cell&    mBackgroundElement;
        const Node_Manager& mNodeManager;
        Matrix< DDRMat >    mGlobalCoordinates;
        Matrix< DDRMat >    mParametricCoordinates;
        static inline bool  gOverrideLinearInterpolation = false;

      public:
        /**
         * Derived node constructor, using other nodal information with background nodes.
         *
         * @param aIndex Node index
         * @param aBackgroundNodes Background nodes
         * @param aParametricCoordinates Parametric coordinates inside the background element
         * @param aGeometryType Geometry type of the background element
         * @param aInterpolationOrder Interpolation order of the background element
         */
        Derived_Node(
                uint                    aIndex,
                const mtk::Cell&        aBackgroundElement,
                const Node_Manager&     aNodeManager,
                const Matrix< DDRMat >& aParametricCoordinates );

        /**
         * Destructor
         */
        ~Derived_Node();

        /**
         * Gets the global coordinates of this node
         *
         * @return Node coordinates
         */
        const Matrix< DDRMat >& get_global_coordinates() const override;

        /**
         * Get the value of a coordinate of this node
         *
         * @param aCoordinateIndex index of the coordinate, obtained from casting the related PDV coordinate type
         * @return Coordinate value
         */
        real get_coordinate_value( uint aCoordinateIndex ) const override;

        /**
         * Gets the parametric coordinates of this derived node relative to its locators
         *
         * @return Parametric coordinates
         */
        const Matrix< DDRMat >& get_parametric_coordinates() const;

        /**
         * Gets the basis nodes of this derived node
         *
         * @return Basis nodes
         */
        Vector< const Background_Node* > get_background_nodes() const;

        /**
         * Gets the locator nodes of this derived node.
         * Locator nodes are the most derived basis nodes that can determine the location of this node.
         * For derived nodes, these are the background nodes, and for intersection nodes these are its parents.
         *
         * @return Locator nodes
         */
        Vector< Basis_Node > get_locator_nodes() const override;

        /**
         * Gets the nodes that provided the field values for this basis node to be created. By default, these is the locator nodes.
         * For linear intersection nodes, these are the parent nodes
         */
        Vector< Basis_Node > get_field_basis_nodes() const override;

        /**
         * Gets if this derived node can be determined that it is on a specific interface without any field evaluation.
         *
         * @param aGeometry Potential interface geometry
         * @return If this node is on the requested interface
         */
        virtual bool is_on_interface( const Geometry& aGeometry ) const;

        /**
         * Gets the number of PDVs on this derived node.
         *
         * @return Number of PDVs
         */
        virtual uint get_num_pdvs() const;

        /**
         * Sets the starting index to be able to use the coordinates of this node as PDVs
         *
         * @param aPDVStartingID The global index of the first PDV on the host
         */
        virtual void set_starting_pdv_id( moris_id aPDVStartingID );

        /**
         * Get the starting global index for the coordinate PDVs
         *
         * @return The global index of the first PDV on the host
         */
        virtual moris_id get_starting_pdv_id() const;

        /**
         * Set the node ID for this node.
         *
         * @param aNodeID Node ID
         */
        virtual void set_id( moris_id aNodeID );

        /**
         * Set the owning processor for this node.
         *
         * @param aNodeOwner Owning processor
         */
        virtual void set_owner( moris_index aNodeOwner );

        /**
         * Get the ID for this node.
         *
         * @return Node ID
         */
        virtual moris_id get_id() const;

        /**
         * Get the owning processor for this node.
         *
         * @return Owning processor
         */
        virtual moris_index get_owner() const;

        /**
         * Gets the background element's nodal coordinates
         */
        Matrix< DDRMat > get_background_element_nodal_coordinates() const;

        /**
         * Gets the geometry type of the background element for this derived node.
         */
        mtk::Geometry_Type get_background_element_geometry_type() const;

        /**
         * Gets the interpolation order of the background element for this derived node.
         */
        mtk::Interpolation_Order get_background_element_interpolation_order() const;

        /**
         * Gets the interpolation order of the locator nodes for this derived node. By default, this is the same as the background element,
         *  but for some derived nodes such as multilinear intersection nodes, this is linear even if the background element is higher order.
         */
        virtual mtk::Interpolation_Order get_locator_interpolation_order() const;

        /**
         * Sets whether to override the locator interpolation order to linear for all derived nodes. This is used for multilinear intersection nodes, which need to have linear interpolation between their locator nodes
         * even if the background element is higher order, in order to be able to use multilinear interpolation for field values.
         *
         * @param aOverride If true, the locator interpolation order will be overridden to linear for all derived nodes
         */
        static void override_linear_interpolation( bool aOverride = true );
    };
}    // namespace moris::gen
