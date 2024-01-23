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
      private:
        Cell< Background_Node > mBackgroundNodes;
        Cell< Derived_Node* > mDerivedNodes;
        bool mMeshGiven = false;

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
        void reset_base_nodes( mtk::Mesh* aMesh );

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
         * Gets a derived node from this manager
         *
         * @param aDerivedNodeIndex Node index (this must be the index of a derived node, or an error will be thrown)
         * @return Derived node
         */
        const Derived_Node& get_derived_node( uint aDerivedNodeIndex ) const;

        /**
         * Gets a trivial node manager. This is useful for allowing geometry construction before
         * the geometry engine has been created.
         *
         * @return Single trivial node manager
         */
        static Node_Manager& get_trivial_instance();

      private:
        /**
         * Deletes all stored nodes, as cleanup operation.
         */
        void delete_all_nodes();
    };
}
