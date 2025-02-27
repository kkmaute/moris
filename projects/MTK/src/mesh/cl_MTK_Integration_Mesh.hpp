/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Mesh.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_

#include <string>
#include "assert.hpp"

#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Block_Set.hpp"
#include "cl_MTK_Side_Set.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Nonconformal_Side_Set.hpp"

namespace moris::mtk
{
    class Integration_Mesh : public virtual Mesh
    {
        // class Cluster_Group;

      protected:
        Vector< moris::mtk::Block_Set * >             mListOfBlocks;
        Vector< moris::mtk::Side_Set * >              mListOfSideSets;
        Vector< moris::mtk::Double_Side_Set * >       mListOfDoubleSideSets;
        Vector< moris::mtk::Nonconformal_Side_Set * > mListOfNonconformalSideSets;
        Vector< moris::mtk::Set * >                   mListOfAllSets;

        Vector< mtk::Double_Side_Cluster * > mDoubleSideClusters;

        map< std::string, moris_index > mSetNameToIndexMap;

        // set by color
        moris_index                                             mMaxColor{};
        Vector< Vector< moris::mtk::Block_Set * > >             mColorToBlockSet;               // outer cell: color | inner cell: block set index
        Vector< Vector< moris::mtk::Side_Set * > >              mColorToSideSet;                // outer cell: color | inner cell: side set index
        Vector< Vector< moris::mtk::Double_Side_Set * > >       mColorToDoubleSideSet;          // outer cell: color | inner cell: double side set index
        Vector< Vector< moris::mtk::Nonconformal_Side_Set * > > mColorToNonconformalSideSet;    // outer cell: color | inner cell: nonconformal side set index
        Vector< Vector< moris::mtk::Set * > >                   mColorToAllSets;                // outer cell: color | inner cell: set index
        
        // ----------------------------------------------------------------------------

      public:
        // ----------------------------------------------------------------------------

        Integration_Mesh() = default;
        // Functions only valid for integration meshes

        // ----------------------------------------------------------------------------

        ~Integration_Mesh() override;

        // ##############################################
        //  Cell Cluster Access
        // ##############################################

        /*
         * Get a cell cluster related to an interpolation
         * cell
         */
        virtual Cell_Cluster const &get_cell_cluster( Cell const &aInterpCell ) const = 0;

        // ----------------------------------------------------------------------------

        /*
         * Get a cell cluster related to an interpolation
         * cell
         */
        virtual Cell_Cluster const &get_cell_cluster( moris_index aInterpCellIndex ) const = 0;

        // ##############################################
        //  MTK Set Access
        // ##############################################

        // ----------------------------------------------------------------------------

        moris::uint
        get_num_sets() const override;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get mesh set by name
         * @param[in] aSetLabel Set label
         */
        moris::mtk::Set *get_set_by_name( const std::string &aSetLabel ) const override;

        // ----------------------------------------------------------------------------

        /**
         * @brief Get set name by set index
         * @param[in] aIndex Set index
         */
        moris::mtk::Set *get_set_by_index( moris_index aIndex ) const override;

        // ----------------------------------------------------------------------------

        void get_Mesh_GEN_map( Vector< moris_index > &aMesh_GEN_map ) override;

        // ----------------------------------------------------------------------------

        /**
         * @brief Get set index by set name
         * @param[in] aSetLabel Set label
         * @return Set index
         */
        moris_index get_set_index_by_name( const std::string &aSetLabel );

        // ----------------------------------------------------------------------------

        /**
         * @brief Get block sets with color
         * @param[in] aColor Set color
         */
        Vector< moris::mtk::Block_Set * > const &get_block_sets_with_color( moris_index const &aColor );

        // ----------------------------------------------------------------------------

        /**
         * @brief Get block set names with color
         * @param[in] aColor    Set color
         * @param[in] aSetNames List of set names
         */
        void get_block_set_names_with_color(
                moris_index const     &aColor,
                Vector< std::string > &aSetNames );

        // ----------------------------------------------------------------------------

        /**
         * @brief Get side sets with color
         * @param[in] aColor Set color
         */
        Vector< moris::mtk::Side_Set * > const &get_side_sets_with_color( moris_index const &aColor );

        // ----------------------------------------------------------------------------

        /**
         * @brief Get double side sets with color
         * @param[in] aColor Set color
         */
        Vector< moris::mtk::Double_Side_Set * > const &get_double_side_sets_with_color( moris_index const &aColor );

        // ----------------------------------------------------------------------------

        /**
         * @brief Get all sets with color
         * @param[in] aColor Set color
         */
        Vector< moris::mtk::Set * > const &get_all_sets_with_color( moris_index const &aColor );

        // ----------------------------------------------------------------------------

        /**
         * @brief Print sets by colors
         */
        void print_sets_by_colors();

        // ----------------------------------------------------------------------------

        /**
         * Get block set names
         */
        virtual Vector< std::string > get_block_set_names() const = 0;

        // ----------------------------------------------------------------------------

        /**
         * Returns the label
         */
        virtual std::string get_block_set_label( moris_index aBlockSetOrdinal ) const;

        // ----------------------------------------------------------------------------

        /**
         * Returns the index given a label
         */
        virtual moris_index get_block_set_index( const std::string &aBlockSetLabel ) const;

        // ----------------------------------------------------------------------------

        // Check if moment fitting is being used or not.
        virtual bool get_moment_fitting_flag() const;

        // ----------------------------------------------------------------------------

        // Get moment fitting poins associated with the entire mesh, 
        // since the points in IP cell parent element space are the same for all IP cells.
        virtual const Matrix< DDRMat >& get_moment_fitting_points() const;

        // ----------------------------------------------------------------------------

        /*
         * Get number of blocks
         */
        moris::uint get_num_blocks() const;

        // ----------------------------------------------------------------------------

        /**
         * Get number of blocks
         * Sometimes num side set * 2. Ask Keenan
         */
        moris::uint get_num_side_set() const;

        // ----------------------------------------------------------------------------

        /**
         * Get number of blocks
         * Sometimes num side set * 2. Ask Keenan
         */
        moris::uint get_num_double_side_set() const;

        // ##############################################
        //  Cell Cluster Access
        // ##############################################

        /*
         * Get cell clusters within a block set
         */
        virtual Vector< Cluster const * > get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const = 0;

        // ----------------------------------------------------------------------------

        virtual uint get_num_cell_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
        {
            MORIS_ERROR( false, "mtk::Integration_Mesh::get_num_cell_cluster_groups() - Not implemented in base class." );
            return 0;
        }

        // ----------------------------------------------------------------------------

        virtual uint get_num_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
        {
            MORIS_ERROR( false, "mtk::Integration_Mesh::get_num_side_cluster_groups() - Not implemented in base class." );
            return 0;
        }

        // ----------------------------------------------------------------------------

        virtual uint get_num_dbl_side_single_side_cluster_groups( const moris_index aDiscretizationMeshIndex ) const
        {
            MORIS_ERROR( false, "mtk::Integration_Mesh::get_num_dbl_side_cluster_groups() - Not implemented in base class." );
            return 0;
        }

        // ##############################################
        //  Side Cluster Access
        // ##############################################

        /**
         * Get side clusters within a side set
         */
        virtual Vector< Cluster const * > get_side_set_cluster( moris_index aSideSetOrdinal ) const = 0;

        // ----------------------------------------------------------------------------
        /**
         * get number of side sets
         */
        virtual uint get_num_side_sets() const = 0;

        // ----------------------------------------------------------------------------
        /**
         * Returns the label
         */
        virtual std::string get_side_set_label( moris_index aSideSetOrdinal ) const = 0;

        // ----------------------------------------------------------------------------
        /**
         * Returns the index given a label
         */
        virtual moris_index get_side_set_index( std::string aSideSetLabel ) const = 0;

        // ##############################################
        //  Double Side Set Cluster Access
        // ##############################################

        /**
         * Returns the number of double sided side sets in the mesh
         */
        virtual uint get_num_double_sided_sets() const = 0;

        // ----------------------------------------------------------------------------
        /**
         * Returns the label
         */
        virtual std::string get_double_sided_set_label( moris_index aSideSetOrdinal ) const = 0;

        // ----------------------------------------------------------------------------
        /**
         * Returns the index given a label
         */
        virtual moris_index get_double_sided_set_index( const std::string &aDoubleSideSetLabel ) const;

        // ----------------------------------------------------------------------------
        /**
         * Returns the double side clusters in the side set
         */
        virtual Vector< Cluster const * > get_double_side_set_cluster( moris_index aSideSetOrdinal ) const = 0;

        // ----------------------------------------------------------------------------
        /**
         * adds a double sided side set to the mesh
         */
        void add_double_side_set( mtk::Double_Side_Set *aDblSideSet );

        // ----------------------------------------------------------------------------
        /**
         * adds a double sided side set to the mesh
         */
        void add_double_sided_cluster( mtk::Double_Side_Cluster * );

        // ----------------------------------------------------------------------------
        /**
         * Save MPC tp hdf5
         */
        void save_MPC_to_hdf5( std::string aFileName );

        // ----------------------------------------------------------------------------

        void create_MPC_maps(
                Matrix< IdMat > &tBSToIPMap,
                Matrix< IdMat > &tIPToIGMap ) const;

        // ----------------------------------------------------------------------------

        void build_hanging_node_MPC(
                Vector< Matrix< IdMat > >  &tBSToIPIds,
                Vector< Matrix< DDRMat > > &tBSToIPWeights,
                const Matrix< IdMat >      &tBSToIPMap,
                const Matrix< IdMat >      &tIPToIGMap );

        // ----------------------------------------------------------------------------
        /**
         * @brief Save elemental T-matrices for the IG-mesh elements to .hdf5
         *
         * @param aFileName
         * @param aNumBsplineMeshes
         */
        void save_elemental_T_matrices_to_file(
                const std::string &aFileName,
                uint               aNumBsplineMeshes = 1 );

        // ----------------------------------------------------------------------------
        /**
         * Save nodal T-matrices for IG-mesh nodes to .hdf5 or .dat file
         */
        void save_IG_global_T_matrix_to_file( const std::string &tFileName );

        // ----------------------------------------------------------------------------

        // elemental
        void get_IG_to_IP_elemental_T_matrices(
                Vector< moris_id >           &aIgCellIds,
                Vector< Matrix< IndexMat > > &aIgToIpIndices,
                Vector< Matrix< DDRMat > >   &aIgToIpTmatrices );

        // ----------------------------------------------------------------------------

        // global
        void get_IG_to_IP_nodal_T_matrices(
                Vector< Matrix< IdMat > >  &aIGtoIPIds,
                Vector< Matrix< DDRMat > > &aIGtoIPWeights,
                uint                        aSetIndex );

        // ----------------------------------------------------------------------------

        // elemental
        void get_IP_to_BS_nodal_T_matrices(
                Vector< Vector< Matrix< IdMat > > >  &aIPtoBSIds,
                Vector< Vector< Matrix< DDRMat > > > &aIPtoBSWeights,
                uint                                  aNumBsplineMeshes );

        // ----------------------------------------------------------------------------

        // global
        void get_IP_to_BS_nodal_T_matrices(
                Vector< Matrix< IdMat > >  &aIPtoBSIds,
                Vector< Matrix< DDRMat > > &aIPtoBSWeights,
                uint                        aSetIndex );

        // ----------------------------------------------------------------------------

        // elemental
        void get_elemental_IG_to_BS_T_matrices(
                Vector< Matrix< IndexMat > > const         &aIgToIpIndices,      // outer cell: IG cell index
                Vector< Matrix< DDRMat > > const           &aIgToIpTmatrices,    // outer cell: IG cell index
                Vector< Vector< Matrix< IdMat > > > const  &aIPtoBSIds,          // outer cell: B-spline mesh index | inner cell: enr. Lagrange BF index
                Vector< Vector< Matrix< DDRMat > > > const &aIPtoBSWeights,      // outer cell: B-spline mesh index | inner cell: enr. Lagrange BF index
                Vector< Vector< Matrix< IdMat > > >        &aIGtoBSIds,          // outer cell: B-spline mesh index | inner cell: IG cell index
                Vector< Vector< Matrix< DDRMat > > >       &aIGtoBSWeights );          // outer cell: B-spline mesh index | inner cell: IG cell index

        // ----------------------------------------------------------------------------

        // global
        void get_IG_to_BS_nodal_T_matrices(
                Vector< Matrix< IdMat > >  &aIGtoBSIds,
                Vector< Matrix< DDRMat > > &aIGtoBSWeights,
                Vector< Matrix< IdMat > >  &aIPtoBSIds,
                Vector< Matrix< DDRMat > > &aIPtoBSWeights,
                Vector< Matrix< IdMat > >  &aIGtoIPIds,
                Vector< Matrix< DDRMat > > &aIGtoIPWeights );

        // ----------------------------------------------------------------------------

        uint get_max_IP_ID_on_set( uint aSetIndex );

        // ----------------------------------------------------------------------------

        uint get_max_IP_index_in_mesh();

        // ----------------------------------------------------------------------------

        void build_sparse_extraction_operator(
                Vector< Matrix< IdMat > >  &aIGtoBSIds,
                Vector< Matrix< DDRMat > > &aIGtoBSWeights,
                Matrix< DDUMat >           &aSparseIndices,
                Matrix< DDRMat >           &aWeights );

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all sets in the mesh.
         * @return The list of all sets in the mesh.
         */
        Vector< moris::mtk::Set * > const &get_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all block sets in the mesh.
         * @return The list of all block sets in the mesh.
         */
        Vector< moris::mtk::Block_Set * > const &get_block_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all side sets in the mesh.
         * @return The list of all side sets in the mesh.
         */
        Vector< moris::mtk::Side_Set * > const &get_side_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all double side sets in the mesh.
         * @return The list of all double side sets in the mesh.
         */
        Vector< moris::mtk::Double_Side_Set * > const &get_double_side_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all sets with a specific color.
         * @return The outer cell is the color index, the inner cell contains the sets.
         */
        Vector< Vector< moris::mtk::Set * > > const &get_color_to_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all block sets with a specific color.
         * @return The outer cell is the color index, the inner cell contains the block sets.
         */
        Vector< Vector< moris::mtk::Block_Set * > > const &get_color_to_block_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all side sets with a specific color.
         * @return The outer cell is the color index, the inner cell contains the side sets.
         */
        Vector< Vector< moris::mtk::Side_Set * > > const &get_color_to_side_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief Get the list of all double side sets with a specific color.
         * @return The outer cell is the color index, the inner cell contains the double side sets.
         */
        Vector< Vector< moris::mtk::Double_Side_Set * > > const &get_color_to_double_side_sets() const;

        // ----------------------------------------------------------------------------
        /**
         * Returns set name to index map
         */
        map< std::string, moris_index > const &get_set_name_to_index_map() const;

        // ----------------------------------------------------------------------------
        /**
         * Returns max color
         */
        moris_index const &get_max_color() const;

        // ----------------------------------------------------------------------------
        /**
         * @brief deletes the ghost visualization sets including ghost block sets and side sets
         *
         */
        void delete_visualization_sets();

      protected:
        void collect_all_sets( bool aSetShape = true );

        void collect_color_to_set_info();

        // ----------------------------------------------------------------------------
        /**
         * @brief Populates the "color to set list". The maximum number of colors should be known before calling this function (see determine_max_color()).
         * @tparam T The Set-type (Block_Set, Side_Set, Double_Side_Set)
         * @param aListOfSets The list of sets to populate the color to set list from
         * @param aColorToSet The color to set list to populate
         */
        template< class T >
        void populate_color_to_set_list( Vector< T > aListOfSets, Vector< Vector< T > > &aColorToSet )
        {
            // clear old data
            aColorToSet.clear();
            aColorToSet.resize( mMaxColor + 1 );

            for ( moris::uint i = 0; i < aListOfSets.size(); i++ )
            {
                Matrix< IndexMat > const &tSetColors = aListOfSets( i )->get_set_colors();
                // iterate through the colors and add to related color grouping
                for ( moris::uint j = 0; j < tSetColors.numel(); j++ )
                {
                    aColorToSet( tSetColors( j ) ).push_back( aListOfSets( i ) );
                }
            }
        }

        // ----------------------------------------------------------------------------
        /**
         * @brief Determines the maximum color number over all mesh sets
         */
        void determine_max_color();
    };    // class mtk::Integration_Mesh

    // ----------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_INTEGRATION_MESH_HPP_ */
