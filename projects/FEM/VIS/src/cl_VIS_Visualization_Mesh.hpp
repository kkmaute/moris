/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Visualization_Mesh.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"

#include "cl_VIS_Vertex_Visualization.hpp"
#include "cl_VIS_Cell_Visualization.hpp"
#include "cl_VIS_Cell_Cluster_Visualization.hpp"

// Logging package
#include "cl_Logger.hpp"

namespace moris::vis
{
    class Visualization_Mesh : public mtk::Mesh
    {
        friend class VIS_Factory;

        /// @brief which cells are to be used for output
        const bool mOnlyPrimary;

        /// @brief global list of all vertices on the mesh
        moris::Cell< mtk::Vertex* > mVertices;

        /// @brief global list of all cells on the mesh
        moris::Cell< mtk::Cell* > mCells;

        /// @brief lists of sets for various types
        moris::Cell< moris::mtk::Set* > mListOfBlocks;
        moris::Cell< moris::mtk::Set* > mListOfSideSets;
        moris::Cell< moris::mtk::Set* > mListOfDoubleSideSets;

        /// @brief sequential list of all sets
        moris::Cell< moris::mtk::Set* > mListOfAllSets;
        Cell< moris_index >             mSetLocalToTypeIndex;

        /// @brief Names of sets in mesh
        moris::Cell< std::string > mAllSetNames;
        moris::Cell< std::string > mBlockSetNames;
        moris::Cell< std::string > mSideSetNames;
        moris::Cell< std::string > mDoubleSideSetNames;

        /// @brief lists of clusters in sets for each type
        moris::Cell< moris::Cell< const mtk::Cluster* > > mClustersOnBlockSets;
        moris::Cell< moris::Cell< const mtk::Cluster* > > mClustersOnSideSets;
        moris::Cell< moris::Cell< const mtk::Cluster* > > mClustersOnDoubleSideSets;

        /// @brief single side clusters that correspond to the double sided side clusters
        moris::Cell< moris::Cell< const mtk::Cluster* > > mLeaderSideClusters;
        moris::Cell< moris::Cell< const mtk::Cluster* > > mFollowerSideClusters;

        map< std::string, moris_index > mSetNameToIndexMap;

        /// @brief flag indicating whether the mesh is ready for output
        bool mMeshIsFinalized = false;

      public:

        /**
         * Constructor
         *
         * @param aOnlyPrimary Whether or not to get only primary cells per set
         */
        explicit Visualization_Mesh( bool aOnlyPrimary );

        /**
         * Destructor, cleans up stored pointers
         */
        ~Visualization_Mesh() override;

        // ##############################################
        //  Cell Cluster Access
        // ##############################################

        /**
         * Gets set names for a given entity type
         *
         * @param aSetEntityRank Entity rank on the set
         * @return Retrieved set names
         */
        moris::Cell< std::string >
        get_set_names( mtk::EntityRank aSetEntityRank ) const override;

        /**
         * Get the number of nodes on the mesh.
         *
         * @return Number of nodes
         */
        uint
        get_num_nodes() const override;

        /**
         * Gets the number of elements on the mesh.
         *
         * @return Number of elements
         */
        uint
        get_num_elems() const override;

        /**
         * Returns all MTK cells in a given block set
         *
         * @param aSetName Set name for getting mtk cells
         * @return Found MTK cells
         */
        moris::Cell< mtk::Cell const * >
        get_set_cells( std::string aSetName ) const override;

        /**
         * Gets all vertices on this mesh
         *
         * @return Vertex pointers
         */
        moris::Cell< moris::mtk::Vertex const * >
        get_all_vertices() const override;

        /**
         * Gets the total number of blocks on this mesh
         *
         * @return Number of blocks
         */
        moris::uint
        get_num_blocks() const;

        /**
         * Gets the total number of sets on this mesh
         *
         * @return Number of sets
         */
        moris::uint
        get_num_sets() const override;

        /**
         * Gets an MTK set by a given index
         *
         * @param aSetIndex
         * @return MTK set at the given index
         */
        mtk::Set*
        get_set_by_index( moris_index aSetIndex ) const override;

        /**
         * Get the set pointer for set with a given name
         *
         * @param aSetLabel name of the set to be retrieved
         * @return moris::mtk::Set* set pointer
         */
        moris::mtk::Set*
        get_set_by_name( std::string aSetLabel ) const override;

        /**
         * Gets a set index by a given name
         *
         * @param aSetLabel Set name
         * @return Set index
         */
        moris_index
        get_set_index_by_name( const std::string& aSetLabel ) const;

        /**
         * Gets element indices in a block set.
         *
         * @param aSetIndex Block set index
         * @return Element indices in the set
         */
        Matrix< IndexMat >
        get_element_indices_in_block_set( uint aSetIndex ) override;

        /**
         * Gets this mesh's type, in this case VIS
         *
         * @return VIS mesh type
         */
        mtk::MeshType
        get_mesh_type() const override;

        /**
         * VIS mesh does not have a communication table, calling this function will raise an error.
         */
        Matrix< IdMat >
        get_communication_table() const override;

        /**
         * Gets the number of spatial dimensions on this mesh
         *
         * @return Number of spatial dimensions
         */
        uint
        get_spatial_dim() const override;

        /**
         * VIS mesh does not support number of entities, calling this function will raise an error.
         */
        uint
        get_num_entities(
                mtk::EntityRank aEntityRank,
                moris_index     aDiscretizationIndex = 0 ) const override;

        /**
         * Gets the nodes connected to a given element
         *
         * @param aElementIndex Element index
         * @return Nodes that are a part of the element
         */
        Matrix< IndexMat >
        get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const override;

        /**
         * Will be moved soon, right now errors out
         */
        Matrix< IndexMat >
        get_entity_connected_to_entity_loc_inds(
                moris_index     aEntityIndex,
                mtk::EntityRank aInputEntityRank,
                mtk::EntityRank aOutputEntityRank,
                moris_index     aDiscretizationIndex = 0 ) const override;
        /**
         * VIS mesh does not support getting neighboring elements, calling this function will raise an error.
         */
        Matrix< IndexMat >
        get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const override;

        /**
         * VIS mesh does not support getting neighboring elements, calling this function will raise an error.
         */
        Matrix< IndexMat >
        get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const override;

        /**
         * Get the spatial coordinates of a node.
         *
         * @param aNodeIndex Node index
         * @return Node coordinates
         */
        Matrix< DDRMat >
        get_node_coordinate( moris_index aNodeIndex ) const override;

        /**
         * Get a global entity ID from an entity rank and local index.
         *
         * @param aEntityIndex Local entity index
         * @param aEntityRank Entity rank
         * @param aDiscretizationIndex Discretization mesh index
         * @return Global entity ID
         */
        moris_id
        get_glb_entity_id_from_entity_loc_index(
                moris_index     aEntityIndex,
                mtk::EntityRank aEntityRank,
                moris_index     aDiscretizationIndex = 0 ) const override;

        /**
         * Get a local entity ID from an entity rank and global ID
         *
         * @param aEntityId Global entity ID
         * @param aEntityRank Entity rank
         * @param aDiscretizationIndex Discretization mesh index
         * @return Local entity index
         */
        moris_index
        get_loc_entity_ind_from_entity_glb_id(
                moris_id        aEntityId,
                mtk::EntityRank aEntityRank,
                moris_index     aDiscretizationIndex = 0 ) const override;

        /**
         * Gets the owner of a node.
         *
         * @param aNodeIndex Node index
         * @return Node owner
         */
        uint
        get_node_owner( moris_index aNodeIndex ) const override;

        /**
         * Gets the owner of an element.
         *
         * @param aElementIndex Element index
         * @return Element owner
         */
        uint
        get_element_owner( moris_index aElementIndex ) const override;

        /**
         * Gets the block set topology for a named set
         *
         * @param aSetName Set name
         * @return CellTopology of the set
         */
        mtk::CellTopology
        get_blockset_topology( const std::string& aSetName ) override;

        /**
         * Gets the integration cell shape for a named set
         *
         * @param aSetName Set name
         * @return CellTopology of the set
         */
        mtk::CellShape
        get_IG_blockset_shape( const std::string& aSetName ) override;

        /**
         * Gets the interpolation cell shape for a named set
         *
         * @param aSetName Set name
         * @return CellTopology of the set
         */
        mtk::CellShape
        get_IP_blockset_shape( const std::string& aSetName ) override;

        /**
         * Gets the sideset element indices and ordinals for a given set name
         *
         * @param aSetName Set name
         * @param aElemIndices Element indices on the sideset
         * @param aSidesetOrdinals Sideset ordinals
         */
        void
        get_sideset_elems_loc_inds_and_ords(
                const std::string&  aSetName,
                Matrix< IndexMat >& aElemIndices,
                Matrix< IndexMat >& aSidesetOrdinals ) const override;

      private:

        /**
         * Constructs a list of all sets as well as a map going from set name to position in the overall set list
         */
        void
        collect_all_sets();

        /**
         * Finalizes this VIS mesh
         */
        void
        finalize();

    };

}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_VIS_VISUALIZATION_MESH_HPP_ */
