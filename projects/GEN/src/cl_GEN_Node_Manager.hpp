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
#include "cl_GEN_Base_Node.hpp"
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
        Cell< Base_Node* > mBaseNodes;
        Cell< Derived_Node* > mDerivedNodes;
        bool mMeshGiven = false;

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
        Node* get_node( uint aNodeIndex );

        /**
         * Gets if the given node index refers to a base node.
         *
         * @return If node is a base node.
         */
        bool is_base_node( uint aNodeIndex );

        /**
         * Gets a base node stored in this manager.
         *
         * @param aBaseNodeIndex Base node index
         * @return Node pointer
         */
        Base_Node* get_base_node( uint aBaseNodeIndex );

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
        Derived_Node* get_derived_node( uint aDerivedNodeIndex );

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
