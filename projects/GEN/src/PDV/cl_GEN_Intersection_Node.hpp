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
#include "cl_GEN_Basis_Node.hpp"

namespace moris::ge
{
    // Forward declare parent node base class
    class Geometry;
    class Parent_Node;

    class Intersection_Node : public Derived_Node
    {
      protected:
        Matrix< DDRMat > mParentVector;
        Matrix< DDSMat > mCoordinateDeterminingADVIDs;

      private:
        Cell< Basis_Node > mParentNodes;
        bool mIsIntersected;
        std::shared_ptr< Geometry > mInterfaceGeometry;
        moris_id mPDVStartingID;
        bool     mPDVStartingIDSet = false;

        moris_id    mNodeID    = -1;
        moris_index mNodeOwner = -1;

      public:
        /**
         * Constructor
         *
         * @param aLocalCoordinate Local coordinate inside of parent edge
         * @param aFirstParentNode First parent node if it is also an intersection node, otherwise nullptr
         * @param aSecondParentNode Second parent node if it is also an intersection node, otherwise nullptr
         * @param aFirstParentNodeIndex Index of the first parent of this node
         * @param aSecondParentNodeIndex Index of the second parent of this node
         * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent of this node
         * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent of this node
         * @param aAncestorNodeIndices Node indices of the ancestors of this intersection node
         * @param aAncestorNodeCoordinates Coordinates of the ancestors of this intersection node
         * @param aAncestorBasisFunction Basis function of the ancestor topology
         */
        Intersection_Node(
                uint                        aNodeIndex,
                const Cell< Node* >&        aBaseNodes,
                const Parent_Node&          aFirstParentNode,
                const Parent_Node&          aSecondParentNode,
                real                        aLocalCoordinate,
                mtk::Geometry_Type          aBaseGeometryType,
                std::shared_ptr< Geometry > aInterfaceGeometry );

        /**
         * Gets if this node's position depends on ADVs. This means either the interface geometry or the parent nodes depend on ADVs.
         *
         * @return ADV dependence
         */
        bool depends_on_advs() override;

        /**
         * Gets the locator nodes of this derived node.
         * For intersection nodes, these are its parents.
         *
         * @return Locator nodes
         */
        const Cell< Basis_Node >& get_locator_nodes() override;

        /**
         * Gets if this intersection node can be determined that it is on a specific interface without any field evaluation.
         *
         * @param aGeometry Potential interface geometry
         * @return If this node is on the requested interface
         */
        bool is_on_interface( Geometry* aGeometry ) override;

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
        real get_local_coordinate();

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

        // FIXME where this is called in the geometry engine, it just uses it to grab an intersection node. We can just provide this information directly instead.
        moris_index
        get_first_parent_node_index()
        {
            return this->get_first_parent_node().get_index();
        }

        moris_index
        get_second_parent_node_index()
        {
            return this->get_second_parent_node().get_index();
        }

      protected:
        /**
         * Computes basic member data for all intersection node derived classes.
         * Must be called by lowest level child class constructors.
         */
        void initialize();

        /**
         * Gets the first parent node of this intersection node.
         *
         * @return First parent node
         */
        Basis_Node& get_first_parent_node();

        /**
         * Gets the second parent node of this intersection node.
         *
         * @return Second parent node
         */
        Basis_Node& get_second_parent_node();

        /**
         * Function for appending to the depending ADV IDs member variable, eliminating duplicate code
         *
         * @param aIDsToAdd IDs to add
         */
        void join_adv_ids( const Matrix< DDSMat >& aIDsToAdd );

        /**
         * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
         * values on each of its ancestors.
         *
         * @param aAncestorIndex Ancestor index
         * @return Local coordinate sensitivity
         */
        virtual real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) = 0;

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its first parent.
         *
         * @return Local coordinate sensitivity
         */
        virtual Matrix< DDRMat > get_dxi_dcoordinate_first_parent() = 0;

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its second parent.
         *
         * @return Local coordinate sensitivity
         */
        virtual Matrix< DDRMat > get_dxi_dcoordinate_second_parent() = 0;

      private:

        /**
         * Determines if the parent nodes are intersected.
         * Used by initialize() to set mIsIntersected. Implementation provided by child class.
         *
         * @return if the parent nodes are intersected
         * @return false if there is no intersection detected
         */
        virtual bool determine_is_intersected() = 0;
    };
}
