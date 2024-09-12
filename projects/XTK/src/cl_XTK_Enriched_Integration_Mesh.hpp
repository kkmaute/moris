/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_XTK_Enriched_Integration_Mesh.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTEGRATION_MESH_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTEGRATION_MESH_HPP_

#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_MTK_Vertex.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "cl_XTK_Field.hpp"
#include <unordered_map>
#include "cl_TOL_Memory_Map.hpp"
#include "cl_XTK_Subphase_Group.hpp"

using namespace moris;

namespace moris::xtk
{
    class Model;
    class Cell_Cluster;
    class Side_Cluster;
    class Interpolation_Cell_Unzipped;
    class Ghost_Stabilization;
    class Enrichment;
    class Cut_Integration_Mesh;
    class Cell_Cluster_Group;
    class Side_Cluster_Group;

    class Enriched_Integration_Mesh : public mtk::Integration_Mesh
    {
        //------------------------------------------------------------------------------

      protected:
        friend class Enrichment;
        friend class Ghost_Stabilization;

        Model                *mModel;
        Cut_Integration_Mesh *mCutIgMesh;

        // mesh indices
        Matrix< IndexMat >                             mBsplineMeshIndices;
        std::unordered_map< moris_index, moris_index > mMeshIndexToLocMeshIndex;
        bool                                           mLocalMeshIndexMapIsSet = false;

        // Cell Clusters
        Vector< std::shared_ptr< xtk::Cell_Cluster > > mCellClusters;

        // cluster groups, first array index corresponds to the discretization mesh index
        Vector< Vector< std::shared_ptr< xtk::Cell_Cluster_Group > > > mCellClusterGroups;
        Vector< Vector< std::shared_ptr< xtk::Side_Cluster_Group > > > mSideClusterGroups;
        Vector< Vector< std::shared_ptr< xtk::Side_Cluster_Group > > > mDblSideClusterGroups;

        // Vertex Set
        std::unordered_map< std::string, moris_index > mVertexSetLabelToOrd;
        Vector< std::string >                          mVertexSetNames;
        Vector< Vector< mtk::Vertex * > >              mVerticesInVertexSet;
        Vector< Matrix< IndexMat > >                   mVertexSetColors;

        // Block sets containing Cell Clusters
        std::unordered_map< std::string, moris_index > mBlockSetLabelToOrd;
        Vector< std::string >                          mBlockSetNames;
        Vector< mtk::CellTopology >                    mBlockSetTopology;
        Vector< Vector< xtk::Cell_Cluster const * > >  mPrimaryBlockSetClusters;
        Vector< Matrix< IndexMat > >                   mBlockSetColors;  /*Bulk phases*/
        Vector< Vector< moris_index > >                mColorsBlockSets; /*transpose of mBlockSetColors*/

        // side sets
        std::unordered_map< std::string, moris_index >           mSideSideSetLabelToOrd;
        Vector< std::string >                                    mSideSetLabels;
        Vector< Vector< std::shared_ptr< xtk::Side_Cluster > > > mSideSets;
        Vector< Matrix< IndexMat > >                             mSideSetColors;  /*Bulk phases of cells attached to side*/
        Vector< Vector< moris_index > >                          mColorsSideSets; /*transpose of mSideSetColors*/

        // double side sets
        std::unordered_map< std::string, moris_index >                  mDoubleSideSetLabelToOrd;
        Vector< std::string >                                           mDoubleSideSetLabels;
        Vector< Vector< std::shared_ptr< mtk::Double_Side_Cluster > > > mDoubleSideSets;                 // outer cell: index of the dbl-side set || inner cell: all dbl-side clusters that live on this dbl-side set
        Vector< Vector< moris_index > >                                 mDoubleSideSetsLeaderIndex;      // outer cell: index of the dbl-side set || inner cell: indices of the leader of side clusters that live on this dbl-side set
        Vector< Vector< moris_index > >                                 mDoubleSideSetsFollowerIndex;    // outer cell: index of the dbl-side set || inner cell: indices of all the follower of side clusters that live on this dbl-side set
        Vector< std::shared_ptr< mtk::Double_Side_Cluster > >           mDoubleSideClusters;
        Vector< std::shared_ptr< xtk::Side_Cluster > >                  mDoubleSideSingleSideClusters;    // lefts and rights of the double side sets
        Matrix< IndexMat >                                              mBulkPhaseToDblSideIndex;
        Vector< Matrix< IndexMat > >                                    mLeaderDoubleSideSetColor;
        Vector< Matrix< IndexMat > >                                    mFollowerDoubleSideSetColor;
        Vector< Vector< moris_index > >                                 mColorLeaderDoubleSideSet;      // transpose of mLeaderDoubleSideSetColor
        Vector< Vector< moris_index > >                                 mColorFollowerDoubleSideSet;    // transpose of mFollowerDoubleSideSetColor

        // Fields
        Vector< xtk::Field >                                               mFields;                        // list of global fields
        Vector< Vector< xtk::Field > >                                     mSideSetFields;                 // outer cell: set ordinal || inner cell: list of fields on that set
        Vector< std::unordered_map< std::string, moris_index > >           mGlobalSetFieldLabelToIndex;    // outer cell: index indicating set type
        Vector< Vector< std::unordered_map< std::string, moris_index > > > mSetWiseFieldLabelToIndex;      // outer cell: set type || inner cell: set ordinal

        // Sub phase index to Cell Cluster Index (these only include the standard cluster i.e. non-ghost clusters.)
        Matrix< IndexMat >              mSubphaseIndexToClusterIndex;      // input: enr IP cell (= cluster) index || output: subphase index
        Vector< Vector< moris_index > > mClusterIndexToSubphaseIndices;    // input: enr IP cell (= cluster) index || output: List of subphase indices in cluster

        // a connectivity pointer used by all transition cells
        mtk::Cell_Info *mCellInfo;

        //------------------------------------------------------------------------------

      public:
        //------------------------------------------------------------------------------
        Enriched_Integration_Mesh( Model *aXTKModel );

        Enriched_Integration_Mesh(
                Model                    *aXTKModel,
                const Matrix< IndexMat > &aBsplineMeshIndices );
        //------------------------------------------------------------------------------
        ~Enriched_Integration_Mesh() override;
        //------------------------------------------------------------------------------
        // MTK Mesh Core Functionality (see base class mtk::Mesh for documentation)
        //------------------------------------------------------------------------------
        mtk::MeshType                               get_mesh_type() const override;
        moris::uint                                 get_spatial_dim() const override;
        uint                                        get_num_entities( mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const override;
        Matrix< IndexMat >                          get_entity_connected_to_entity_loc_inds( moris_index aEntityIndex, mtk::EntityRank aInputEntityRank, mtk::EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const override;
        Matrix< IndexMat >                          get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const override;
        moris_id                                    get_glb_entity_id_from_entity_loc_index( moris_index aEntityIndex, mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const override;
        std::unordered_map< moris_id, moris_index > get_vertex_glb_id_to_loc_vertex_ind_map() const override;
        moris_index                                 get_loc_entity_ind_from_entity_glb_id( moris_id aEntityId, mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const override;
        Vector< mtk::Vertex const * >               get_all_vertices() const override;
        Matrix< IdMat >                             get_entity_connected_to_entity_glob_ids( moris_id aEntityId, mtk::EntityRank aInputEntityRank, mtk::EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const override;
        Matrix< DDRMat >                            get_node_coordinate( moris_index aNodeIndex ) const override;
        mtk::Vertex                                &get_mtk_vertex( moris_index aVertexIndex ) override;
        mtk::Vertex const                          &get_mtk_vertex( moris_index aVertexIndex ) const override;
        mtk::Cell                                  &get_writable_mtk_cell( moris_index aElementIndex ) override;
        mtk::Cell                                  &get_mtk_cell( moris_index aElementIndex ) override;
        mtk::Cell const                            &get_mtk_cell( moris_index aElementIndex ) const override;
        Matrix< IdMat >                             get_communication_table() const override;
        Vector< std::string >                       get_set_names( mtk::EntityRank aSetEntityRank ) const override;
        mtk::CellTopology                           get_blockset_topology( const std::string &aSetName ) override;
        mtk::CellShape                              get_IG_blockset_shape( const std::string &aSetName ) override;
        mtk::CellShape                              get_IP_blockset_shape( const std::string &aSetName ) override;
        Matrix< IndexMat >                          get_set_entity_loc_inds( mtk::EntityRank aSetEntityRank, const std::string &aSetName ) const override;
        Matrix< IndexMat >                          get_element_indices_in_block_set( uint aSetIndex ) override;
        void                                        get_sideset_elems_loc_inds_and_ords( const std::string &aSetName, Matrix< IndexMat > &aElemIndices, Matrix< IndexMat > &aSidesetOrdinals ) const override;
        moris_id                                    get_max_entity_id( mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const override;
        uint                                        get_node_owner( moris_index aNodeIndex ) const override;
        uint                                        get_element_owner( moris_index aElementIndex ) const override;

        //------------------------------------------------------------------------------
        // end mesh core functions
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // MTK Integration Mesh Functions
        // see base class mtk::Integration_Mesh for documentation
        //------------------------------------------------------------------------------
        mtk::Cell_Cluster const       &get_cell_cluster( mtk::Cell const &aInterpCell ) const override;
        mtk::Cell_Cluster const       &get_cell_cluster( moris_index aInterpCellIndex ) const override;
        Vector< std::string >          get_block_set_names() const override;
        std::string                    get_block_set_label( moris_index aBlockSetOrdinal ) const override;
        moris_index                    get_block_set_index( const std::string &aBlockSetLabel ) const override;
        Vector< mtk::Cluster const * > get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const override;
        Matrix< IndexMat >             get_block_set_colors( moris_index aBlockSetOrdinal ) const;
        Vector< mtk::Cluster const * > get_side_set_cluster( moris_index aSideSetOrdinal ) const override;
        Matrix< IndexMat >             get_side_set_colors( moris_index aSideSetOrdinal ) const;
        uint                           get_num_side_sets() const override;
        std::string                    get_side_set_label( moris_index aSideSetOrdinal ) const override;
        moris_index                    get_side_set_index( std::string aSideSetLabel ) const override;
        uint                           get_num_double_sided_sets() const override;
        std::string                    get_double_sided_set_label( moris_index aSideSetOrdinal ) const override;
        moris_index                    get_double_sided_set_index( const std::string &aDoubleSideSetLabel ) const override;
        Vector< mtk::Cluster const * > get_double_side_set_cluster( moris_index aSideSetOrdinal ) const override;
        Matrix< IndexMat >             get_double_side_set_colors( moris_index aSideSetOrdinal ) const;
        uint                           get_sidesets_num_faces( Vector< moris_index > aSideSetIndex ) const override;
        //------------------------------------------------------------------------------
        // end integration mesh functions
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // Additional Get/Set Functions
        //------------------------------------------------------------------------------

        void                                                        setup_mesh_index_map();
        moris_index                                                 get_local_mesh_index_xtk( moris_index const &aDiscretizationMeshIndex ) const;
        Matrix< IndexMat >                                          get_enriched_mesh_indices() const override;
        uint                                                        get_num_interpolation_types() const;
        uint                                                        get_num_cell_cluster_groups( const moris_index aDiscretizationMeshIndex ) const override;
        uint                                                        get_num_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const override;
        uint                                                        get_num_dbl_side_single_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const override;
        Vector< std::shared_ptr< xtk::Cell_Cluster_Group > > const &get_cell_cluster_groups( const moris_index aDiscretizationMeshIndex ) const;
        Vector< std::shared_ptr< xtk::Side_Cluster_Group > > const &get_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const;

        //------------------------------------------------------------------------------

        Vector< xtk::Cell_Cluster const * > const &
        get_xtk_cell_clusters_in_block_set( moris_index aBlockSetOrdinal ) const;

        /**
         * @return Interface side set name
         */
        std::string
        get_interface_side_set_name(
                moris_index aGeomIndex,
                moris_index aBulkPhaseIndex0,
                moris_index aBulkPhaseIndex1 );

        //------------------------------------------------------------------------------
        /**
         * @return Double sided interface side set name
         */
        std::string
        get_dbl_interface_side_set_name(
                moris_index aBulkPhaseIndex0,
                moris_index aBulkPhaseIndex1 );

        //------------------------------------------------------------------------------

        /**
         * @return Primary cell local indices in a block set
         */
        Matrix< IndexMat >
        get_block_entity_loc_inds( const std::string &aSetName ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Creates a double sided interface between the two bulk phases.
         * This function creates additional dbl sided interfaces. By default,
         * the enriched integration mesh creates only the low-leader to high-follower
         * dbl-sided interfaces. This functions allows creation of low-follower to high-leader
         * interfaces.
         * The creation is done in groups as the finalization of side-sets is computationally expensive, but can be done in bulk.
         * Hence, avoid using this function in loops with a large number of iterations.
         * @param[in] aLeaderBulkPhaseIndex List of leader bulk phase indices
         * @param[in] aFollowerBulkPhaseIndex List of follower bulk phase indices
         */
        void
        create_dbl_sided_interface_sets(
                Vector< moris_index > aLeaderBulkPhaseIndex,
                Vector< moris_index > aFollowerBulkPhaseIndex );

        //------------------------------------------------------------------------------
        // Output/ Viz Functions
        //------------------------------------------------------------------------------
        /**
         * @brief For cleanup when writing to an exodus file (note: in general this should not be used because
         * sets may not be always empty through an optimization run)
         */
        void
        deactivate_empty_sets();

        //------------------------------------------------------------------------------

        /**
         * @brief Deactivate empty side sets in mesh.
         */
        void
        deactivate_empty_side_sets();

        //------------------------------------------------------------------------------

        /**
         * @brief Deactivate empty block sets in mesh.
         */
        void
        deactivate_empty_block_sets();

        //------------------------------------------------------------------------------
        /**
         * @brief Create basis support fields
         * @return Cell of field names for basis support (1xNumBasisFunc)
         */

        Vector< std::string >
        create_basis_support_fields( Matrix< DDRMat > const &aProbeSpheres );

        void
        create_bg_cell_id_field();

        void create_subphase_fields();

        void
        write_subphase_neighborhood( const std::string &aFile );

        void
        create_cell_id_fields();

        void
        add_mesh_field_real_scalar_data_loc_inds(
                std::string const      &aFieldName,
                mtk::EntityRank const  &aFieldEntityRank,
                Matrix< DDRMat > const &aFieldData ) override;
        //------------------------------------------------------------------------------

        /**
         * @brief Write mesh
         */
        void
        write_mesh(
                moris::Parameter_List *aParamList );

        void
        create_union_block(
                Vector< std::string > const &aBlocks,
                std::string                  aNewBlock,
                Matrix< IndexMat > const    &aNewBlockColor );

        void
        create_union_side_set(
                Vector< std::string > const &aSideSets,
                std::string                  aNewSideSet,
                Matrix< IndexMat > const    &aNewSideSetColor );

        void
        deactivate_all_blocks_except_selected(
                Vector< std::string > const &aBlockSetsToKeep );

        void
        deactivate_all_side_sets_except_selected(
                Vector< std::string > const &aSideSetsToKeep );

        //------------------------------------------------------------------------------
        // Memory Map
        //------------------------------------------------------------------------------
        /**
         * @brief Memory map of the enriched integration mesh
         * @return Memory map
         */
        moris::Memory_Map
        get_memory_usage();

        //------------------------------------------------------------------------------
        // Additional Field Functions
        //------------------------------------------------------------------------------

        /**
         * @brief Get the labels of the fields
         *
         * @param aEntityRank
         * @param aSetOrdinal ordinal of the set or bulk-phase of the block for which the field names should be obtained (leave empty for global names)
         * @return Vector< std::string > list of labels of the fields
         */
        Vector< std::string >
        get_field_names(
                mtk::EntityRank          aEntityRank,
                const moris::moris_index aSetOrdinal = MORIS_INDEX_MAX );

        //------------------------------------------------------------------------------

        /**
         * @brief Create a field
         * @param[in] aLabel Field label
         * @param[in] aEntityRank Field entity rank
         * @param[in]aSetOrdinal Bulk phase field defined over
         * aBulkphaseIndex of MORIS_INDEX_MAX results in a field over all phases
         */
        moris::moris_index
        create_field(
                const std::string &aLabel,
                mtk::EntityRank    aEntityRank,
                moris::moris_index aSetOrdinal = MORIS_INDEX_MAX );

        //------------------------------------------------------------------------------

        /**
         * @brief Given a field name and rank, gets the field index(ordinal)
         * @return Field index
         */
        moris::moris_index
        get_field_index(
                const std::string       &aLabel,
                const mtk::EntityRank    aEntityRank,
                const moris::moris_index aSetOrdinal = MORIS_INDEX_MAX );

        //------------------------------------------------------------------------------

        /**
         * @brief Add field data to created field
         * @param[in] aFieldIndex Field index (use fn:get_field_index(...) to get)
         * @param[in] aEntityRank Field entity rank
         * @param[in] aFieldData Field data
         */
        void
        add_field_data(
                moris::moris_index      aFieldIndex,
                mtk::EntityRank         aEntityRank,
                Matrix< DDRMat > const &aFieldData,
                moris::moris_index      aSetOrdinal = MORIS_INDEX_MAX );

        //------------------------------------------------------------------------------

        /**
         * @brief Given a field index and field entity rank, get the field data
         * @param[in] aFieldIndex Field index (use fn:get_field_index(...) to get)
         * @param[in] aEntityRank Field entity rank
         * @return Field data
         */
        Matrix< DDRMat > const &
        get_field_data(
                moris::moris_index aFieldIndex,
                mtk::EntityRank    aEntityRank,
                moris::moris_index aSetOrdinal = MORIS_INDEX_MAX ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Convert entity indices to entity ids
         * @param[in] aIndices Entity indices
         * @param[in] aEntityRank Entity rank
         * @return Entity ids
         */
        Matrix< IdMat > convert_indices_to_ids(
                Matrix< IndexMat > const &aIndices,
                mtk::EntityRank           aEntityRank ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Convert entity ids to entity indices
         * @param[in] aIds Entity ids
         * @param[in] aEntityRank Entity rank
         * @return Entity indices
         */
        Matrix< IndexMat > convert_ids_to_indices(
                Matrix< IdMat > const &aIds,
                mtk::EntityRank        aEntityRank ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Get MTK cells from their indices
         * @param[in] aCellIndices Cell indices
         * @return MTK cells
         */
        Vector< mtk::Cell const * >
        get_mtk_cells_loc_inds( Matrix< IndexMat > const &aCellIndices );

        //------------------------------------------------------------------------------

        /**
         * @brief Get MTK vertices from their indices
         * @param[in] aVertexIndices Vertex indices
         * @return MTK vertices
         */
        Vector< mtk::Vertex const * >
        get_mtk_vertices_loc_inds( Matrix< IndexMat > const &aVertexIndices );
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // Accessor functions of XTK specific data structures
        //------------------------------------------------------------------------------
        /**
         * @brief Get the XTK cell cluster implementation
         * @param[in] aInterpCell Interpolation MTK cell
         * @return XTK cell cluster
         * NOTE: requires there is only one cell cluster associated with the interpolation cell
         */
        xtk::Cell_Cluster const &
        get_xtk_cell_cluster( mtk::Cell const &aInterpCell ) const;

        //------------------------------------------------------------------------------
        /**
         * @brief Get the XTK cell cluster implementation
         * @param[in] aInterpCell Interpolation MTK cell
         * @return XTK cell cluster
         * NOTE: requires there is only one cell cluster associated with the interpolation cell
         */
        xtk::Cell_Cluster const &
        get_xtk_cell_cluster( moris_index aInterpCellIndex ) const;

        //------------------------------------------------------------------------------
        // Printing
        //------------------------------------------------------------------------------
        void print() const;
        void print_general() const;
        void print_vector_clusters( moris::uint aVerbosityLevel = 0 ) const;
        void print_block_sets( moris::uint aVerbosityLevel = 0 ) const;
        void print_side_sets( moris::uint aVerbosityLevel = 0 ) const;
        void print_double_side_sets( moris::uint aVerbosityLevel = 0 ) const;
        void print_double_side_clusters( moris::uint aVerbosityLevel = 0 ) const;

        //--------------------------------------------------------------------------------
        // Utilities for manipulating sets
        //--------------------------------------------------------------------------------
        /**
         * @brief Create a single side set from a double sided side set
         * @param[in] aDblSideSetIndex Double side set index
         * @param[in] aSideSetName New side set name
         * @param[in] aCollectSets Perform collect operation upon completion of function.
         * @return Side set index
         */
        moris_index
        create_side_set_from_dbl_side_set(
                moris_index const &aDblSideSetIndex,
                std::string const &aSideSetName,
                bool               aCollectSets = true );

        /**
         * @brief Create a block set from a single sided side set
         * @param[in] aDblSideSetIndex Double side set index
         * @param[in] aBlockSetName New block set name
         * @param[in] aCellTopo New block set cell topology
         * @param[in] bool set this to true if the blockset's only purpose is for visualization
         * @return Block set index
         */
        moris_index
        create_block_set_from_cells_of_side_set(
                moris_index const       &aSideSetIndex,
                std::string const       &aBlockSetName,
                mtk::CellTopology const &aCellTopo,
                bool                     aCreateOnlyForVis = false );

        //------------------------------------------------------------------------------

        void
        setup_cluster_groups();

        //------------------------------------------------------------------------------

        void
        visualize_cluster_measures();

        //------------------------------------------------------------------------------

        void
        visualize_cluster_group_measures( const bool aWriteBsplineClusterInfo );

        //------------------------------------------------------------------------------

      protected:
        //------------------------------------------------------------------------------
        // Parallel functions
        //------------------------------------------------------------------------------
        moris_id allocate_entity_ids(
                moris::size_t   aNumReqs,
                mtk::EntityRank aEntityRank );

        //------------------------------------------------------------------------------

        void
        commit_double_side_set( moris_index const &aDoubleSideSetIndex );

        //------------------------------------------------------------------------------

        void
        commit_double_side_set( const Vector< moris_index > &aDoubleSideSetIndexList );

        //------------------------------------------------------------------------------

        void
        commit_side_set( moris_index const &aSideSetIndex );

        //------------------------------------------------------------------------------

        void
        commit_side_set( const Vector< moris_index > &aSideSetIndexList );

        //------------------------------------------------------------------------------

        void
        commit_block_set( moris_index const &aBlockSetIndex );

        //------------------------------------------------------------------------------

      public:
        void
        communicate_sets_of_type( const mtk::SetType aSetType );

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------

        uint
        get_num_owned_cells() const;

        //------------------------------------------------------------------------------

        /**
         * @brief Fill enriched integration mesh with Cell clusters, setting their active and void IG cells
         *
         */
        void
        setup_cell_clusters();

        void
        setup_cell_clusters_new();

        //------------------------------------------------------------------------------

        /**
         * @brief Group clusters with the same active bulk-phase into "block-" sets
         *
         */
        void
        setup_blockset_with_cell_clusters();

        //------------------------------------------------------------------------------

        /**
         * @brief Set the up side set clusters object
         *
         */
        void
        setup_side_set_clusters();

        //------------------------------------------------------------------------------

        void
        setup_double_side_set_clusters();

        //------------------------------------------------------------------------------

        void
        setup_color_to_set();

        //------------------------------------------------------------------------------

        // TODO!

        void
        setup_cell_cluster_groups();

        void
        setup_side_cluster_groups();

        void
        setup_side_cluster_groups_old();

        void
        setup_dbl_side_cluster_groups();

        void
        setup_dbl_side_cluster_groups_old();

        //------------------------------------------------------------------------------

        void
        setup_double_sided_interface_sides();

        //------------------------------------------------------------------------------

        void
        declare_interface_double_side_sets();

        //------------------------------------------------------------------------------

        moris_index
        get_dbl_side_set_index(
                moris_index aPhase0,
                moris_index aPhase1 );

        //------------------------------------------------------------------------------

        void
        create_interface_double_side_sets_and_clusters();

        //------------------------------------------------------------------------------

        void
        add_side_to_cluster(
                const std::shared_ptr< xtk::Side_Cluster > &aSideCluster,
                moris_index                                 aCellIndex,
                moris_index                                 aSideOrdinal );

        //------------------------------------------------------------------------------

        Vector< std::string >
        split_set_name_by_bulk_phase( const std::string &aBaseName );

        //------------------------------------------------------------------------------

        Vector< std::string >
        split_set_name_by_child_no_child( const std::string &aBaseName );

        //------------------------------------------------------------------------------

        Vector< moris_index >
        register_vertex_set_names( Vector< std::string > const &aVertexSetNames );

        //------------------------------------------------------------------------------

        Vector< moris_index >
        register_block_set_names_with_cell_topo(
                Vector< std::string > const &aBlockSetNames,
                mtk::CellTopology            aBlockTopology );

        //------------------------------------------------------------------------------

        void
        set_block_set_colors(
                moris_index const        &aBlockSetIndex,
                Matrix< IndexMat > const &aBlockSetColors );

        //------------------------------------------------------------------------------

        Vector< moris_index >
        register_side_set_names( Vector< std::string > const &aSideSetNames );

        //------------------------------------------------------------------------------

        void
        set_side_set_colors(
                moris_index const        &aSideSetIndex,
                Matrix< IndexMat > const &aSideSetColors );

        //------------------------------------------------------------------------------

        Vector< moris_index >
        register_double_side_set_names( Vector< std::string > const &aDblSideSetNames );

        //------------------------------------------------------------------------------

        void
        set_double_side_set_colors(
                moris_index const        &aDblSideSetIndex,
                Matrix< IndexMat > const &aLeaderSideColors,
                Matrix< IndexMat > const &aFollowerSideColors );

        //------------------------------------------------------------------------------

        void
        setup_interface_side_sets();

        //------------------------------------------------------------------------------

        void
        declare_interface_side_sets();

        //------------------------------------------------------------------------------

        void
        create_interface_side_sets_and_clusters();

        //------------------------------------------------------------------------------

        Vector< moris_index >
        declare_interface_vertex_sets();

        //------------------------------------------------------------------------------

        void
        set_vertex_set_color(
                moris_index const        &aVertexSetIndex,
                Matrix< IndexMat > const &aVertexSetColors );

        //------------------------------------------------------------------------------

        void
        construct_color_to_set_relationship(
                Vector< Matrix< IndexMat > > const &aSetColors,
                Vector< Vector< moris_index > >    &aColorToSetIndex );

        //------------------------------------------------------------------------------

        void
        create_interface_side_sets_from_interface_double_side_set(
                moris_index const &aBulkphase0,
                moris_index const &aBulkphase1 );

        //------------------------------------------------------------------------------

        uint
        get_num_elements_in_side_set( const std::string &aSideSetName );

        //------------------------------------------------------------------------------
        // Internal Additional Field Functions
        //------------------------------------------------------------------------------

        /*
         * Returns an index in the data structure for a given entity rank (i.e. NODE = 0)
         */
        moris_index
        get_entity_rank_field_index( const mtk::EntityRank aEntityRank );

        //------------------------------------------------------------------------------
        /**
         * @return  whether a field exists or not
         */
        bool
        field_exists(
                const std::string       &aLabel,
                const mtk::EntityRank    aEntityRank,
                const moris::moris_index aSetOrdinal );

    };    // class Enriched_Integration_Mesh

}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTEGRATION_MESH_HPP_ */
