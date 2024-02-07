/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Background_Mesh.hpp
 *
 */

#ifndef SRC_XTK_CL_XTK_BACKGROUND_MESH_HPP_
#define SRC_XTK_CL_XTK_BACKGROUND_MESH_HPP_

#include <unordered_map>
#include <utility>

// Mesh Includes:
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Enums.hpp"

// Assertion Includes:
#include "fn_assert.hpp"

// other
#include "cl_TOL_Memory_Map.hpp"
#include "containers.hpp"

// Linear Algebra Includes
#include "cl_Matrix.hpp"
#include "fn_print.hpp"
#include "fn_isvector.hpp"

// XTK Includes:
#include "cl_XTK_External_Mesh_Data.hpp"
#include "cl_XTK_Downward_Inheritance.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Cell_CM.hpp"
#include "cl_MTK_Vertex_XTK_Impl.hpp"
#include "cl_MTK_Vertex_Interpolation_XTK_Impl.hpp"

// general geometry engine class
#include "cl_GEN_Geometry_Engine.hpp"

namespace xtk
{
    class Background_Mesh
    {
      public:
        Background_Mesh(){};

        Background_Mesh( moris::mtk::Interpolation_Mesh* aMeshData );

        Background_Mesh(
                moris::mtk::Interpolation_Mesh* aMeshData,
                moris::gen::Geometry_Engine*     aGeometryEngine );

        ~Background_Mesh();
        /*!
         * Get number of entities in the background mesh and
         * the number of entities XTK has created
         *
         * NOTE: this function includes all nodes (background nodes and interface nodes)
         * but only includes background elements
         *
         */
        moris::size_t
        get_num_entities( mtk::EntityRank aEntityRank ) const;

        /*!
         * Get number of entities in the background mesh
         */
        moris::size_t
        get_num_entities_background( mtk::EntityRank aEntityRank ) const;

        /*!
         * get the owning proc of a vertex
         */
        moris::moris_index
        get_vertex_owner( moris::moris_index aVertexIndex ) const;

        /*!
         * Get the vertex ownership of both background and XTK created nodes
         */
        moris::Matrix< moris::IndexMat >
        get_vertices_owner( moris::Matrix< moris::IndexMat > const & aVertexIndices ) const;

        /*!
         * Get an mtk vertex
         */
        Vector< moris::mtk::Vertex const * >
        get_mtk_vertices( Matrix< IndexMat > const & aVertexIndices );

        /*!
         * Get an mtk vertex
         */
        moris::mtk::Vertex&
        get_mtk_vertex( moris::moris_index aVertexIndex );

        /*!
         * Get an mtk vertex as the XTK child vertex.
         * allows internal use of more vertex functionality
         */
        moris::mtk::Vertex_XTK&
        get_mtk_vertex_xtk( moris::moris_index aVertexIndex );

        /*!
         * Get an mtk cell
         */
        moris::mtk::Cell&
        get_mtk_cell( moris::moris_index aCellIndex );

        /*!
         * get a const mtk cell
         */
        moris::mtk::Cell const &
        get_mtk_cell( moris::moris_index aCellIndex ) const;

        /*!
         *  Get an offset which is the first id allocated (assumed ids are grouped)
         */
        moris::moris_id
        allocate_entity_ids( moris::size_t aNumReqs,
                mtk::EntityRank            aChildEntityRank );

        /*!
         * Get the first available index
         */
        moris::moris_index
        get_first_available_index( mtk::EntityRank aEntityRank ) const;

        /*!
         * Increment the first available index by aNewFirstAvailableIndex
         */
        void
        update_first_available_index(
                moris::size_t   aNewFirstAvailableIndex,
                mtk::EntityRank aEntityRank );

        /*!
         * Create a batch of new nodes
         */
        void
        batch_create_new_nodes(
                Vector< moris_index > const &                    aNewNodeIds,
                Vector< moris_index > const &                    aNewNodeIndices,
                Vector< moris_index > const &                    aNewNodeOwningProc,
                Vector< moris::Matrix< moris::DDRMat > > const & aNewNodeCoordinates );

        /*!
         * Batch create a copy node (used in unzipping)
         */
        void
        batch_create_new_nodes_as_copy_of_other_nodes(
                moris::Matrix< moris::IndexMat > const & aExistingNodeIndices,
                moris::Matrix< moris::IndexMat > const & aNewNodeIds,
                moris::Matrix< moris::IndexMat > const & aNewNodeIndices );

        void
        allocate_external_node_to_child_mesh_associations();

        /*!
         * Create association of external nodes to their child mesh index,
         * The node index vector does not necessarily need to be only external nodes
         * but only the ones which are external will be associated to a child mesh
         */
        void
        associate_external_nodes_to_child_mesh(
                moris::moris_index                       aChildMeshIndex,
                moris::Matrix< moris::IndexMat > const & aNodeIndices );

        /*!
         * Get the child mesh indices that a node belongs to.
         * Only implemented for nodes created during the decomposition process, called the
         * external nodes
         */
        moris::Matrix< moris::IndexMat >
        get_node_child_mesh_assocation( moris::moris_index aNodeIndex ) const;

        /*!
         * Get the global entity id from local index
         */
        moris::size_t
        get_glb_entity_id_from_entity_loc_index(
                moris::size_t       aEntityIndex,
                mtk::EntityRank     aEntityRank,
                moris_index const & aMeshIndex = 0 ) const;
        /*!
         * From a vector of entity ids and ranks, return the global ids of these entities
         */
        moris::Matrix< moris::IdMat >
        get_glb_entity_id_from_entity_loc_index_range(
                moris::Matrix< moris::IndexMat > const & tEntityIndices,
                mtk::EntityRank                          aEntityRank ) const;
        /*!
         * Convert local entity indices to global entity ids
         */
        void
        convert_loc_entity_ind_to_glb_entity_ids(
                mtk::EntityRank                   aEntityRank,
                moris::Matrix< moris::IndexMat >& aEntityIndices ) const;

        /*!
         * Returns whether a node was in the original background mesh
         */
        bool
        is_background_node( moris::moris_index aNodeIndex );

        /*!
         * Returns whether a cell was in the original background mesh
         */
        bool
        is_background_cell( moris::moris_index aCellIndex ) const;

        /*!
         * Return all node coordinates ordered by local indices
         */
        moris::Matrix< moris::DDRMat >
        get_all_node_coordinates_loc_inds() const;

        /*!
         * Return a coordinate matrix for the specified node indices
         */
        moris::Matrix< moris::DDRMat >
        get_selected_node_coordinates_loc_inds( moris::Matrix< moris::IndexMat > const & aNodeIndices ) const;

        /*!
         * Returns the local to global map (only implemented in for nodes)
         */
        moris::Matrix< moris::IdMat >
        get_local_to_global_map( mtk::EntityRank aEntityRank ) const;

        /*!
         * Return a vector of all non-intersected elements'
         * element to node connectivity
         */
        moris::Matrix< moris::IdMat >
        get_full_non_intersected_node_to_element_glob_ids() const;

        /*!
         * Package and return all intersected element to node connectivity
         * sorted by phase
         */
        Vector< moris::Matrix< moris::IdMat > >
        get_non_intersected_element_to_node_by_phase( moris::uint aNumPhases );

        /*!
         * Return all ids of non-intersected elements
         */
        moris::Matrix< moris::IdMat >
        get_all_non_intersected_elements() const;

        /*!
         * Return all ids of non-intersected elements
         */
        Vector< moris::Matrix< moris::IdMat > >
        get_all_non_intersected_elements_by_phase( uint aNumPhases ) const;

        /*!
         * Return all non-intersected element proc local indices
         */
        moris::Matrix< moris::IndexMat >
        get_all_non_intersected_elements_loc_inds() const;

        // -------------------------------------------------------------------

        // Functions for setting up downard inheritance, where downward
        // downward inheritance is the relationship between a parent element
        // and its children elements
        // -------------------------------------------------------------------
        void
        register_new_downward_inheritance( Vector< std::pair< moris::moris_index, moris::moris_index > > const & aNewElementToChildMeshPairs );

        /*
         * used after clean up of child mesh and deletion of child meshes to recompute the downward inheritance
         */
        void
        setup_downward_inheritance( Cut_Mesh& aCutMesh );

        /*!
         * returns whether a given entity has any children entities
         */
        bool
        entity_has_children( moris::size_t aEntityIndex,
                mtk::EntityRank            aEntityRank ) const;

        /*!
         * returns the child mesh index of entity with children
         */
        moris::moris_index const & child_mesh_index(
                moris::size_t   aEntityIndex,
                mtk::EntityRank aEntityRank );

        // -------------------------------------------------------------------
        // Functions related to setting and accessing interface node information
        // -------------------------------------------------------------------

        void
        initialize_interface_node_flags(
                moris::size_t const & aNumNodes,
                moris::size_t const & aNumGeometry );

        /*!
         * Allocates additional space in interface node flags
         */
        void
        allocate_space_in_interface_node_flags(
                moris::size_t aNumNodes,
                moris::size_t aNumGeometry );

        /*!
         * Marks a node as an interface node for a given geometry index
         */
        void
        mark_node_as_interface_node(
                moris::moris_index aNodeIndex,
                moris::size_t      aGeomIndex );

        /*!
         * Marks a node as an interface node for a given geometry index
         */
        void
        mark_nodes_as_interface_node_loc_inds(
                moris::Matrix< moris::IndexMat > aNodeIndices,
                moris::size_t                    aGeomIndex );

        /*!
         * Returns whether a node is an interface node for a given geometry index
         */
        bool
        is_interface_node(
                moris::moris_index aNodeIndex,
                moris::size_t      aGeomIndex ) const;

        /*!
         * get the interface nodes with respect to a given geometry index
         */
        moris::Matrix< moris::IndexMat >
        get_interface_nodes_loc_inds( moris::moris_index aGeometryIndex ) const;

        /*!
         * get the interface nodes with respect to a given geometry index
         */
        Vector< moris::Matrix< moris::IdMat > >
        get_interface_nodes_loc_inds() const;

        /*!
         * get the interface nodes with respect to a given geometry index
         */
        moris::Matrix< moris::IdMat >
        get_interface_nodes_glob_ids( moris::moris_index aGeometryIndex ) const;

        //
        moris::Matrix< moris::IndexMat >
        restrict_vertex_list_to_owned_by_this_proc_loc_inds( moris::Matrix< moris::IndexMat > const & aNodeIndexList ) const;

        /*!
         * get the interface nodes with respect to a given geometry index
         */
        Vector< moris::Matrix< moris::IdMat > >
        get_interface_nodes_glob_ids();

        void
        print_interface_node_flags();

        void
        print_vertex_map();

        //--------------------------------------------------------------------------------
        // Memory Functions
        //--------------------------------------------------------------------------------

        /*!
         * @brief get the memory usage of background mesh
         */
        moris::Memory_Map
        get_memory_usage();

        // -------------------------------------------------------------------
        // Functions related to setting and accessing element phase indices
        // -------------------------------------------------------------------

        /*!
         * Allocate space for element phase indices
         */
        void
        initialize_element_phase_indices( moris::size_t const & aNumElements );

        /*!
         * Set the phase index value of element with element index. This is relative to each geometry.
         */
        void
        set_element_phase_index(
                moris::size_t const & aElementIndex,
                moris::size_t const & aElementPhaseIndex );

        /*!
         * Get the phase index value of element with element index
         */
        moris::moris_index const &
        get_element_phase_index( moris::size_t const & aElementIndex ) const;

        /*!
         * Multiple version of above
         */
        moris::Matrix< moris::IndexMat >
        get_element_phase_inds( moris::Matrix< moris::DDSTMat > const & aElementInds );

        moris_index
        get_loc_entity_ind_from_entity_glb_id(
                moris_id        aEntityId,
                mtk::EntityRank aEntityRank ) const;

        std::unordered_map< moris_id, moris_index >
        get_vertex_glb_id_to_loc_vertex_ind_map() const;

        /*!
         * add child element to mtk cells
         */
        void
        add_child_element_to_mtk_cells(
                moris::moris_index aElementIndex,
                moris::moris_index aElementId,
                moris::moris_index aElementOwner,
                moris::moris_index aCMElementIndex,
                Child_Mesh*        aChildMeshPtr );

        void
        add_new_cell_to_mesh( moris::mtk::Cell* aCell );

        void
        add_cells_to_global_to_local_map(
                Matrix< IndexMat > const & aCellIndices,
                Matrix< IdMat > const &    aCellIds );

        /*!
         * return a pointer to a cell
         */
        const moris::mtk::Cell*
        get_child_element_mtk_cell_ptr( moris::moris_index aElementIndex ) const;

        /*!
         * return a reference to an mtk cell
         */
        moris::mtk::Cell&
        get_child_element_mtk_cell( moris::moris_index aElementIndex )
        {
            auto tIter = mChildMtkCellMap.find( aElementIndex );

            MORIS_ASSERT( mChildMtkCellMap.find( aElementIndex ) != mChildMtkCellMap.end(), "Element index is not in the map" );

            moris::moris_index tIndex = tIter->second;
            return *mChildMtkCells( tIndex );
        }

        /*!
         * return a const reference to an mtk cell
         */
        moris::mtk::Cell const &
        get_child_element_mtk_cell( moris::moris_index aElementIndex ) const;

        // -------------------------------------------------------------------
        // Access underlying mesh data functions

        /*!
         * return the background mesh data
         */
        moris::mtk::Interpolation_Mesh&
        get_mesh_data();

        /*!
         * return const reference to background mesh data
         */
        moris::mtk::Interpolation_Mesh const &
        get_mesh_data() const;

        // -------------------------------------------------------------------

        mtk::CellTopology
        get_parent_cell_topology() const;

        Matrix< IdMat > const &
        get_communication_table() const;

        void
        add_proc_to_comm_table( moris_index aProcRank );

        void
        remove_cells_from_mesh( Vector< moris_index > const & aCellsToRemove,
                Vector< moris_index >&                        aOldIndexToNewCellIndex );
        /*!
         * Sets up the entity local to global maps
         */
        void
        setup_local_to_global_maps();

      private:
        // Background mesh data
        moris::mtk::Interpolation_Mesh* mMeshData;

        // External Entity information
        // The background mesh remains constant, and every new entity created is stored within XTK
        // child meshes
        Mesh_External_Entity_Data mExternalMeshData;

        // Downward inheritance pairs (links elements in XTK mesh to indices in Child Meshes)
        Downward_Inheritance< moris::moris_index, moris::moris_index > mElementDownwardInheritance;

        // Local to Global Id Entity Matrix
        Vector< Vector< moris::moris_index > > mEntityLocaltoGlobalMap;

        // communication map
        moris::Matrix< IdMat > mCommunicationMap;

        // Elements constructed by the decomposition process mtk Cells
        Mini_Map< moris_id, moris_index > mChildMtkCellMap; /* To go from cell index to location in child cell ptrs*/
        Vector< moris::mtk::Cell* >  mChildMtkCells;

        // Vertex constructed by the decomposition process
        std::unordered_map< moris_id, moris_index > mVertexGlbToLocalMap;
        Vector< moris::mtk::Vertex_XTK >       mXtkMtkVertices;

        // Associate external node indices to the child meshes they belong to
        // Row - External node index
        // Col - Child Mesh Index
        moris::Matrix< moris::IndexMat > mNodeIndexToChildMeshIndex;

        // Element Phase Index ordered by processor local indices
        moris::Matrix< moris::IndexMat > mElementPhaseIndex;

        // Nodal Phase Index
        // Note the exact phase value is located in the geometry index.
        // Columns - Geometry Index
        // Rows - Node Index
        // If Val = 0; This means the node is not an interface node for a given geometry
        // If Val = 1; This means the node is an interface node for a given geometry
        moris::Matrix< moris::IndexMat > mInterfaceNodeFlag;

        moris::Matrix< moris::DDRMat >
        get_all_node_coordinates_loc_inds_background_mesh() const;

        /*!
         * Allocate space in the downward inheritance map, one for each element in background mesh
         */
        void intialize_downward_inheritance();

        /*!
         * initialize xtk mtk vertices
         */
        void
        initialize_background_mesh_vertices();
    };
}    // namespace xtk

#endif /* SRC_XTK_CL_XTK_BACKGROUND_MESH_HPP_ */
