/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node.hpp
 *
 */

#pragma once

#include "cl_GEN_Derived_Node.hpp"

namespace moris::ge
{
    // Forward declare necessary classes
    class Geometry;
    class Basis_Node;
    class Parent_Node;

    class Intersection_Node : public Derived_Node
    {
      private:
        Cell< Basis_Node > mParentNodes;
        std::shared_ptr< Geometry > mInterfaceGeometry;
        moris_id mPDVStartingID;
        bool     mPDVStartingIDSet = false;

        moris_id    mNodeID    = -1;
        moris_index mNodeOwner = -1;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBaseNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry
         */
        Intersection_Node(
                uint                        aNodeIndex,
                const Cell< Node* >&        aBaseNodes,
                const Parent_Node&          aFirstParentNode,
                const Parent_Node&          aSecondParentNode,
                real                        aLocalCoordinate,
                mtk::Geometry_Type          aBackgroundGeometryType,
                mtk::Interpolation_Order    aBackgroundInterpolationOrder,
                std::shared_ptr< Geometry > aInterfaceGeometry );

        /**
         * Gets if this node's position depends on ADVs. This means either the interface geometry or the parent nodes depend on ADVs.
         *
         * @return ADV dependence
         */
        bool depends_on_advs() const override;

        /**
         * Gets the locator nodes of this derived node.
         * For intersection nodes, these are its parents.
         *
         * @return Locator nodes
         */
        const Cell< Basis_Node >& get_locator_nodes() const override;

        /**
         * Gets if this intersection node can be determined that it is on a specific interface without any field evaluation.
         *
         * @param aGeometry Potential interface geometry
         * @return If this node is on the requested interface
         */
        bool is_on_interface( Geometry* aGeometry ) const override;

        /**
         * Returns if the parent edge is intersected (if the local coordinate of the intersection lies between
         * -1 and 1)
         *
         * @return If the edge is intersected
         */
        bool parent_edge_is_intersected();

        /**
         * Returns if the first parent used to create this node is on the geoemtry interface already.
         *
         * @return If the first parent is on the interface
         */
        bool is_first_parent_on_interface();

        /**
         * Returns if the second parent used to create this node is on the geometry interface already.
         *
         * @return If the second parent is on the interface
         */
        bool is_second_parent_on_interface();

        /**
         * Gets the local coordinate of this intersection node inside of the parent edge.
         *
         * @return Local coordinate
         */
        real get_local_coordinate() const;

        /**
         * Get the value of a coordinate of this node
         *
         * @param aCoordinateIndex index of the coordinate, obtained from casting the related PDV coordinate type
         * @return Coordinate value
         */
        real get_coordinate_value( uint aCoordinateIndex );

        /**
         * Gets the number of PDVs on this intersection node.
         *
         * @return Number of PDVs
         */
        uint get_num_pdvs();

        /**
         * Sets the starting index to be able to use the intersection coordinates of this node as PDVs
         *
         * @param aPDVStartingID The global index of the first PDV on the host
         */
        void set_starting_pdv_id( moris_id aPDVStartingID );

        /**
         * Get the starting global index for the intersection coordinate PDVs
         *
         * @return The global index of the first PDV on the host
         */
        moris_id get_starting_pdv_id();

        /**
         * Set the node ID for this node.
         *
         * @param aNodeID Node ID
         */
        void set_id( moris_id aNodeID );

        /**
         * Set the owning processor for this node.
         *
         * @param aNodeOwner Owning processor
         */
        void set_owner( moris_index aNodeOwner );

        /**
         * Get the ID for this node.
         *
         * @return Node ID
         */
        moris_id get_id();

        /**
         * Get the owning processor for this node.
         *
         * @return Owning processor
         */
        moris_index get_owner();

      protected:

        /**
         * Gets the first parent node of this intersection node.
         *
         * @return First parent node
         */
        const Basis_Node& get_first_parent_node() const;

        /**
         * Gets the second parent node of this intersection node.
         *
         * @return Second parent node
         */
        const Basis_Node& get_second_parent_node() const;

        /**
         * Function for appending to the depending ADV IDs member variable, eliminating duplicate code
         *
         * @param aCombinedIDs Combined IDs
         * @param aIDsToAdd IDs to add
         */
        static void join_adv_ids(
                Matrix< DDSMat >&       aCombinedIDs,
                const Matrix< DDSMat >& aIDsToAdd );
    };
}