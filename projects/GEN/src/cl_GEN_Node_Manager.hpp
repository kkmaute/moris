/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Phase_Table.hpp
 *
 */

#pragma once

#include "cl_Cell.hpp"
#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Derived_Node.hpp"

// Forward declare mtk mesh
namespace moris::mtk
{
    class Mesh;
}

namespace moris::ge
{
    class Node_Manager
    {
      protected:
        bool mMeshGiven = false;

      private:
        Cell< Background_Node > mBackgroundNodes;
        Cell< Derived_Node* > mDerivedNodes;

      public:
        /**
         * Constructor with a given mesh
         *
         * @param aMesh Mesh pointer for creating background nodes
         */
        explicit Node_Manager( mtk::Mesh* aMesh );

        /**
         * Destructor
         */
        ~Node_Manager();

        /**
         * Resets the stored background node information with a given MTK mesh.
         *
         * @param aMesh Input mesh with background nodes
         */
        void reset_background_nodes( mtk::Mesh* aMesh );

        /**
         * Gets the number of background nodes stored by the node manager.
         *
         * @return Number of background nodes
         */
        uint get_number_of_background_nodes();

        /**
         * Gets the total number of nodes stored by the node manager.
         *
         * @return Total number of nodes
         */
        uint get_total_number_of_nodes();

        /**
         * Gets a node stored in this manager.
         *
         * @param aNodeIndex Node index
         * @return Node pointer
         */
        const Node& get_node( uint aNodeIndex );

        /**
         * Gets if the given node index refers to a background node.
         *
         * @return If node is a background node.
         */
        bool is_background_node( uint aNodeIndex ) const;

        /**
         * Gets a background node stored in this manager.
         *
         * @param aBackgroundNodeIndex background node index
         * @return Node pointer
         */
        Background_Node& get_background_node( uint aBackgroundNodeIndex );

        /**
         * Creates a new derived node and adds it to this manager with the given parameters.
         *
         * @param aBackgroundNodes Background nodes
         * @param aParametricCoordinates Parametric coordinates inside the background element
         * @param aGeometryType Geometry type of the background element
         * @param aInterpolationOrder Interpolation order of the background element
         */
        void create_derived_node(
                const Cell< Node* >&     aBackgroundNodes,
                const Matrix< DDRMat >&  aParametricCoordinates,
                mtk::Geometry_Type       aGeometryType,
                mtk::Interpolation_Order aInterpolationOrder );

        /**
         * Adds a derived node to this manager.
         *
         * @param aDerivedNode New derived node
         */
        void add_derived_node( Derived_Node* aDerivedNode );

        /**
         * Gets a derived node from this manager (const version).
         *
         * @param aDerivedNodeIndex Node index (this must be the index of a derived node, or an error will be thrown)
         * @return Derived node
         */
        const Derived_Node& get_derived_node( uint aDerivedNodeIndex ) const;

        /**
         * Update the queued intersection node with its node ID and node owner.
         *
         * @param aNodeIndex Node index
         * @param aNodeID Node ID
         * @param aNodeOwner Node owner
         */
        void update_derived_node(
                uint        aNodeIndex,
                moris_id    aNodeID,
                moris_index aNodeOwner );

        /**
         * Gets the coordinate value at a specified index of a stored derived node.
         *
         * @param aNodeIndex Node index
         * @param aCoordinateIndex Coordinate index
         * @return Coordinate value
         */
        real get_node_coordinate_value(
                uint aNodeIndex,
                uint aCoordinateIndex );

        /**
         * Gets if a stored derived node depends on ADVs
         *
         * @param aNodeIndex Node index
         * @return Dependency on ADVs
         */
        bool node_depends_on_advs( uint aNodeIndex );

        /**
         * Gets the nubmer of PDVs contained by a stored derived node.
         *
         * @param aNodeIndex Node index
         * @return Number of PDVs on the node
         */
        uint get_number_of_derived_node_pdvs( uint aNodeIndex );

        /**
         * Sets the starting PDV ID of a stored derived node.
         *
         * @param aNodeIndex Node index
         * @param aStartingPDVID Starting PDV ID
         */
        void set_derived_node_starting_pdv_id(
                uint     aNodeIndex,
                moris_id aStartingPDVID );

        /**
         * Gets the starting PDV ID of a stored derived node.
         *
         * @param aNodeIndex Node index
         * @return Starting PDV ID
         */
        moris_id get_derived_node_starting_pdv_id( uint aNodeIndex );

        /**
         * Gets the ID of a stored derived node, or -1 if it has not been set.
         *
         * @param aNodeIndex Node index
         * @return Node ID (not PDV ID)
         */
        moris_id get_derived_node_id( uint aNodeIndex );

        /**
         * Gets the owning processor index of a stored derived node.
         *
         * @param aNodeIndex Node index
         * @return Owning processor
         */
        moris_index get_derived_node_owner( uint aNodeIndex );

        /**
         * Appends the sensitivities of the given node's global coordinates with respect to ADVs.
         *
         * @param aNodeIndex Derived node index
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         */
        void append_dcoordinate_dadv_from_derived_node(
                uint                    aNodeIndex,
                Matrix< DDRMat >&       aCoordinateSensitivities,
                const Matrix< DDRMat >& aSensitivityFactor );

        /**
         * Gets the ADV IDs that determine the coordinates of a derived node
         *
         * @param aNodeIndex Node index
         * @return ADV ID vector
         */
        Matrix< DDSMat > get_coordinate_determining_adv_ids_from_derived_node( uint aNodeIndex );

        /**
         * Gets a trivial node manager. This is useful for allowing geometry construction before
         * the geometry engine has been created.
         *
         * @return Single trivial node manager
         */
        static Node_Manager& get_trivial_instance();

      private:

        /**
         * Gets a derived node from this manager. Helper function for getting references in this class.
         *
         * @param aDerivedNodeIndex Node index (this must be the index of a derived node, or an error will be thrown)
         * @return Derived node
         */
        Derived_Node& get_derived_node( uint aDerivedNodeIndex );

        /**
         * Deletes all stored nodes, as cleanup operation.
         */
        void delete_all_nodes();
    };
}
