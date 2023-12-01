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
#include "cl_GEN_Node.hpp"
#include "cl_GEN_Intersection_Node.hpp"

// Forward declare mtk mesh
namespace moris::mtk
{
    class Mesh;
}

namespace moris::ge
{
    // Forward declare derived node
    class Base_Node;
    class Derived_Node;
    class Intersection_Node;

    class Node_Manager
    {
      private:
        Cell< Base_Node* > mBaseNodes;
        Cell< Derived_Node* > mDerivedNodes;
        Cell< Intersection_Node* > mIntersectionNodes;

      public:
        /**
         * Constructor with a given mesh
         *
         * @param aMesh Mesh pointer for creating base nodes
         */
        explicit Node_Manager( mtk::Mesh* aMesh );

        /**
         * Destructor
         */
        ~Node_Manager();

        /**
         * Resets the stored base node information with a given MTK mesh.
         *
         * @param aMesh Input mesh with base nodes
         */
        void reset_base_nodes( mtk::Mesh* aMesh );

        /**
         * Gets the number of base nodes currently stored by this node manager.
         *
         * @return Number of base nodes
         */
        uint get_number_of_base_nodes();

        /**
         * Gets the total number of nodes stored by the node manager.
         *
         * @return Total number of nodes
         */
        uint get_total_number_of_nodes();

        /**
         * Gets a node stored in this manager.
         *
         * @param aIndex Node index
         * @return Node pointer
         */
        Base_Node* get_base_node( uint aIndex );

        /**
         * Adds a derived node to this manager.
         *
         * @param aDerivedNode New derived node
         */
        void add_derived_node( Derived_Node* aDerivedNode );

        /**
         * Gets a derived node from this manager
         *
         * @param aIndex Node index (this must be the index of a derived node, or an error will be thrown)
         * @return Derived node
         */
        Derived_Node* get_derived_node( uint aIndex );

        /**
         * Adds an intersection node to this manager.
         *
         * @param aIntersectionNode New derived node
         */
        void add_intersection_node( Intersection_Node* aIntersectionNode );

        /**
         * Gets all the intersection nodes stored in the node manager that depend on ADVs.
         * Useful for the PDV host manager.
         *
         * @return Vector of intersection nodes
         */
        const Cell< Intersection_Node* >& get_intersection_nodes();

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
