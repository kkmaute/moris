/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Intersection_Mesh.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_MESH_HPP_

#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris::mtk
{
    class Intersection_Mesh : public mtk::Integration_Mesh
    {
      private:
        moris::mtk::Integration_Mesh*    mBackGroundMesh;     /**< Background mesh that is not intersected */
        moris::mtk::Intersection_Detect* mIntersectionDetect; /**< Holds all the information and data created */

      public:
        // ----------------------------------------------------------------------------
        /*
         * Default constructor
         */
        Intersection_Mesh( moris::mtk::Integration_Mesh* aBackGroundMesh,
                moris::mtk::Intersection_Detect*         aIntersectionDetect );

        // ----------------------------------------------------------------------------
        /*
         * Default destructor
         */
        ~Intersection_Mesh() override;

        // ----------------------------------------------------------------------------
        /**
         * Returns the spatial dimension of the mesh
         */
        uint
        get_spatial_dim() const override
        {
            return mBackGroundMesh->get_spatial_dim();
        }

        // ----------------------------------------------------------------------------
        /*
         * Get number of nodes/elements in the mesh
         * @param[ in ] aEntityRank Rank of the entity
         * @param[ in ] aIndex A mesh index
         */
        moris::uint
        get_num_entities( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const override;

        // ----------------------------------------------------------------------------
        /*
         * adds side sets and block sets created in the process of intersection to the mesh
         * @param[ in ] aEntityRank Rank of the entity
         */
        Vector< std::string >
        get_set_names( enum EntityRank aSetEntityRank ) const override;

        // ----------------------------------------------------------------------------
        /*
         * Determines if a given side set name belongs to intersection data or background mesh
         * @param[ in ] aSetName Name of the set
         */
        bool
        is_intersection_side_set( const std::string& aSetName ) const;

        // ----------------------------------------------------------------------------

        /**
         *Determines if a given block set name belongs to intersection data or background mesh
         * @param[ in ] aSetName Name of the set
         */
        bool
        is_intersection_block_set( const std::string& aSetName ) const;

        // ----------------------------------------------------------------------------
        /*
         * get side set elements indices and their respective ordinal
         * @param[ in ] aSetName Name of the set
         * @param[ out ] aElemIndices Indices of the element in the set
         * @param[ out ] aSidesetOrdinals Side ordinals of the element that make the side set
         */
        void
        get_sideset_elems_loc_inds_and_ords(
                const std::string&  aSetName,
                Matrix< IndexMat >& aElemIndices,
                Matrix< IndexMat >& aSidesetOrdinals ) const override;

        // ----------------------------------------------------------------------------
        /*
         * Get element indices in a block set
         * @param[ in ] aSetIndex Index of the set
         */
        Matrix< IndexMat >
        get_element_indices_in_block_set( uint aSetIndex ) override;

        // ----------------------------------------------------------------------------
        /*
         * Get element ids in a block set
         * @param[ in ] aSetIndex Index of the set
         */
        Matrix< IdMat >
        get_element_ids_in_block_set( uint aSetIndex ) override;

        // ----------------------------------------------------------------------------
        /*
         * Get the node coordinates for an index
         * @param[ in ] aNodeIndex Index of the node
         */
        moris::Matrix< moris::DDRMat >
        get_node_coordinate( moris_index aNodeIndex ) const override;

        // ----------------------------------------------------------------------------
        /*
         * Get global entity id from entity local index
         * @param[ in ] aEntityIndex Index of the entity
         * @param[ in ] aEntityRank Rank of the entity
         * @param[ in ] aMeshIndex Index of the mesh
         */
        moris::moris_id
        get_glb_entity_id_from_entity_loc_index( moris_index aEntityIndex,
                enum EntityRank                              aEntityRank,
                const moris_index                            aMeshIndex ) const override;

        // ----------------------------------------------------------------------------
        /**
         * Returns block set toplogy
         * @param[ in ] aSetName Name of the set
         */
        enum CellTopology get_blockset_topology( const std::string& aSetName ) override;

        // ----------------------------------------------------------------------------
        /**
         * Get nodes connected to an element
         * @param[ in ] aElementIndex An element index
         */
        Matrix< IndexMat >
        get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const override;

        // ##############################################
        //  Pure virtual function ( not implemented )
        // ##############################################

        // ----------------------------------------------------------------------------
        /**
         * Get nodes connected to an element
         * @param[ in ] aInterpCell An Interpolation cell
         */
        Cell_Cluster const &
        get_cell_cluster( Cell const & aInterpCell ) const override
        {
            MORIS_ERROR( false, "Not Implemented for Intersection Mesh" );
            mtk::Cell_Cluster* tDummyCells = nullptr;
            return *tDummyCells;
        }

        // ----------------------------------------------------------------------------
        /**
         * Get nodes connected to an element
         * @param[ in ] aInterpCell An Interpolation cell
         */
        Cell_Cluster const &
        get_cell_cluster( moris_index aInterpCellIndex ) const override
        {
            MORIS_ERROR( false, "Not Implemented for Intersection Mesh" );
            mtk::Cell_Cluster* tDummyCells = nullptr;
            return *tDummyCells;
        }

        // ----------------------------------------------------------------------------

        /**
         * Get block set names
         */
        Vector< std::string >
        get_block_set_names() const override
        {
            MORIS_ERROR( false, "Not Implemented for Intersection Mesh" );
            return Vector< std::string >( 0 );
        }

        // ----------------------------------------------------------------------------

        /**
         * Get cell clusters in set
         * @param[ in ] aBlockSetOrdinal Ordinal of the block set
         */
        Vector< Cluster const * >
        get_cell_clusters_in_set( moris_index aBlockSetOrdinal ) const override
        {
            MORIS_ERROR( false, "Not Implemented for Intersection Mesh" );
            return Vector< Cluster const * >( 0 );
        }

        // ----------------------------------------------------------------------------
        /**
         * Get side clusters within a side set
         * @param[ in ] aSideSetOrdinal Ordinal of the side set
         */
        Vector< Cluster const * >
        get_side_set_cluster( moris_index aSideSetOrdinal ) const override
        {
            MORIS_ERROR( false, "Not Implemented for Intersection Mesh" );
            return Vector< Cluster const * >( 0 );
        }

        // ----------------------------------------------------------------------------
        /**
         * get number of side sets
         */
        uint
        get_num_side_sets() const override
        {
            MORIS_ERROR( false, "get_num_side_sets(), Not Implemented for Intersection Mesh" );
            return 0;
        }

        // ----------------------------------------------------------------------------
        /**
         * Returns the label
         * @param[ in ] aSideSetOrdinal Ordinal of the side set
         */
        std::string
        get_side_set_label( moris_index aSideSetOrdinal ) const override
        {
            MORIS_ERROR( false, " get_side_set_label, Not Implemented for Intersection Mesh" );
            return std::string();
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the index given a label
         * @param[ in ] aSideSetLabel Label of the side set
         */
        moris_index
        get_side_set_index( std::string aSideSetLabel ) const override
        {
            MORIS_ERROR( false, "get_side_set_index, Not Implemented for Intersection Mesh" );
            return 0;
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the number of double sided side sets in the mesh
         */
        uint
        get_num_double_sided_sets() const override
        {
            MORIS_ERROR( false, "get_num_double_sided_sets, Not Implemented for Intersection Mesh" );
            return 0;
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the label
         * @param[ in ] aSideSetOrdinal Ordinal of the side set
         */
        std::string
        get_double_sided_set_label( moris_index aSideSetOrdinal ) const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return std::string();
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the label
         */
        MeshType
        get_mesh_type() const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return MeshType::UNDEFINED;
        }

        // ----------------------------------------------------------------------------
        /*!
         * Get local indices of the entities connected to a specific entity
         * @param[ in ] aEntityIndex Entity index
         * @param[ in ] aInputEntityRank Entity rank
         * @param[ in ] aOutputEntityRank Entity rank of the output
         * @param[ in ] aDiscretizationIndex A discretization index
         */
        Matrix< IndexMat >
        get_entity_connected_to_entity_loc_inds(
                moris_index       aEntityIndex,
                enum EntityRank   aInputEntityRank,
                enum EntityRank   aOutputEntityRank,
                const moris_index aDiscretizationIndex = 0 ) const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return Matrix< IndexMat >( 0, 0 );
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns node owner
         * @param[ in ] aNodeIndex Index of the node
         */
        uint get_node_owner( moris_index aNodeIndex ) const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return 0;
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the element owner
         * @param[ in ] aNodeIndex Index of the node
         */
        uint get_element_owner( moris_index aNodeIndex ) const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return 0;
        }

        // ----------------------------------------------------------------------------
        /*!
         * get communication table
         */
        Matrix< IdMat >
        get_communication_table() const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return Matrix< IndexMat >( 0, 0 );
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns integration block set shape
         * @param[ in ] aSetName Name of the set
         */
        enum CellShape
        get_IG_blockset_shape( const std::string& aSetName ) override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return CellShape::UNDEFINED;
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns interpolation block set shape
         * @param[ in ] aSetName Name of the set
         */
        enum CellShape
        get_IP_blockset_shape( const std::string& aSetName ) override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return CellShape::UNDEFINED;
        }

        // ----------------------------------------------------------------------------
        /*!
         * Get double sided cluster
         * @param[ in ] aSideSetOrdinal Ordinal of the side set
         */
        Vector< Cluster const * >
        get_double_side_set_cluster( moris_index aSideSetOrdinal ) const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return Vector< Cluster const * >( 0 );
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the label
         * @param[ in ] aElementIndex Index of the element
         */
        Matrix< IndexMat >
        get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const override
        {
            MORIS_ERROR( false, "get_double_sided_set_label, Not Implemented for Intersection Mesh" );
            return Matrix< IndexMat >( 0, 0 );
        }

        // ----------------------------------------------------------------------------
        /*!
         * Returns the label
         * @param[ in ]
         */
        Matrix< IndexMat >
        get_set_entity_loc_inds(
                enum EntityRank    aSetEntityRank,
                const std::string& aSetName ) const override
        {
            return Matrix< IndexMat >( 0, 0 );
        }
    };
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERSECTION_MESH_HPP_ */
