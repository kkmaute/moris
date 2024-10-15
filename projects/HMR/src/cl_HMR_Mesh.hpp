/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh.hpp
 *
 */

#pragma once

#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "cl_Vector.hpp"    //CNT/src
#include "cl_MTK_Enums.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Mesh_Core.hpp"    //MTK/src

namespace moris::hmr
{
    //-------------------------------------------------------------------------------

    // forward declaration of HMR Field object
    class Field;

    class Database;

    //-------------------------------------------------------------------------------
    /**
     * \brief mesh interface class
     */
    class Mesh : public virtual mtk::Mesh
    {
        //! describing label
        std::string mLabel;

        // Dummy B-spline meshes
        BSpline_Mesh_Base* mDummyBSplineMesh = nullptr;

        // Stored Lagrange mesh
        Lagrange_Mesh_Base* mMesh = nullptr;

        // give access to the member data to stk based interpolation and integration meshes
        friend class Interpolation_Mesh_HMR;
        friend class Integration_Mesh_HMR;

        // Global to Local Entity Map
        // Vector(0) - Node Global to Local
        // Vector(1) - Edge Global to Local
        // Vector(2) - Face Global to Local
        // Vector(3) - Cell Global to Local
        // Keep in mind not all of these are created
        Vector< std::unordered_map< moris_id, moris_index > > mEntityGlobalToLocalMap;

        //-------------------------------------------------------------------------------

      public:
        //-------------------------------------------------------------------------------

        /**
         * mesh constructor, to be called from HMR
         */
        Mesh(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeOrder,
                uint                        aLagrangePattern );

        Mesh(
                std::shared_ptr< Database > aDatabase,
                uint                        aOrder,
                uint                        aLagrangePattern,
                BSpline_Mesh_Base*          aDummyBSplineMesh );

        Mesh(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeOrder,
                uint                        aLagrangePattern,
                uint                        aBSplineOrder,
                uint                        aBSplinePattern );

        Mesh(
                std::shared_ptr< Database > aDatabase,
                uint                        aLagrangeMeshIndex );

        //-------------------------------------------------------------------------------

        /**
         * destructor
         */
        ~Mesh() override;

        //-------------------------------------------------------------------------------

        /**
         * return the type of this mesh
         */
        mtk::MeshType
        get_mesh_type() const override
        {
            return mtk::MeshType::HMR;
        }

        //-------------------------------------------------------------------------------

        /**
         * provides a moris::Matrix< DDUMat > containing the IDs this mesh has
         * to communicate with
         */
        Matrix< IdMat > get_communication_table() const override;

        //-------------------------------------------------------------------------------

        /**
         * creates a new field pointer that is linked to this mesh
         */
        std::shared_ptr< Field > create_field(
                const std::string& aLabel,
                uint               aBSplineOrder );

        //-------------------------------------------------------------------------------

        /**
         * returns a pointer to the underlying lagrange mesh
         */
        Lagrange_Mesh_Base*
        get_lagrange_mesh()
        {
            return mMesh;
        }

        //-------------------------------------------------------------------------------

        /**
         * returns a the underlying background mesh pattern on which this mesh was build on
         */
        uint
        get_lagrange_mesh_pattern()
        {
            return mMesh->get_activation_pattern();
        }

        //-------------------------------------------------------------------------------

        /**
         * return a pointer to the database
         */
        std::shared_ptr< Database >
        get_database()
        {
            return mDatabase;
        }

        //-------------------------------------------------------------------------------
        // Functions for MTK
        //-------------------------------------------------------------------------------

        /**
         * returns the number of dimensions of this mesh
         */
        uint get_spatial_dim() const override;

        //-------------------------------------------------------------------------------

        /**
         * Gets the polynomial order of this mesh
         *
         * @return Polynomial degree
         */
        uint get_order() override;

        //-------------------------------------------------------------------------------

        /**
         * Gets the polynomial order of an underlying B-spline mesh
         *
         * @param aDiscretizationIndex B-spline mesh index
         * @return Polynomial degree
         */
        uint get_discretization_order( uint aDiscretizationIndex = 0 ) override;

        //-------------------------------------------------------------------------------

        uint get_num_entities(
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        uint get_num_elements_including_aura() const;

        //-------------------------------------------------------------------------------

        uint get_num_nodes() const override;

        //-------------------------------------------------------------------------------

        uint get_num_edges() const override;

        //-------------------------------------------------------------------------------

        uint get_num_faces() const override;

        //-------------------------------------------------------------------------------

        uint get_num_elems() const override;

        /**
         * Gets element indices in a block set.
         *
         * @param aSetIndex Block set index
         * @return Element indices in the set
         */
        Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex ) override;

        /**
         * Gets the element IDs in a block set.
         *
         * @param aSetIndex Block set index
         * @return Element IDs in the set
         */
        Matrix< IdMat > get_element_ids_in_block_set( uint aSetIndex ) override;

        //-------------------------------------------------------------------------------

        uint get_max_num_coeffs_on_proc( uint aBSplineMeshIndex ) const override;

        /**
         * Gets the indices of the B-splines which form the basis of the given node.
         *
         * @param aNodeIndex Node index
         * @param aBSplineMeshIndex B-spline mesh index
         * @return B-spline indices
         */
        Matrix< IndexMat > get_coefficient_indices_of_node(
                uint aNodeIndex,
                uint aBSplineMeshIndex ) override;

        /**
         * Gets the IDs of the B-splines which form the basis of a given node.
         *
         * @param aNodeIndex Node index
         * @param aBSplineMeshIndex B-spline mesh index
         * @return B-spline IDs
         */
        Matrix< IdMat > get_coefficient_IDs_of_node(
                uint aNodeIndex,
                uint aBSplineMeshIndex ) override;

        Matrix< IdMat > get_coefficient_owners_of_node(
                uint aNodeIndex,
                uint aBSplineMeshIndex ) override;

        Matrix< IdMat > get_coefficient_ijkl_IDs_of_node(
                uint aNodeIndex,
                uint aBSplineMeshIndex ) override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_entity_connected_to_entity_loc_inds(
                moris_index     aEntityIndex,
                mtk::EntityRank aInputEntityRank,
                mtk::EntityRank aOutputEntityRank,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_entity_connected_to_entity_glob_ids(
                moris_index     aEntityId,
                mtk::EntityRank aInputEntityRank,
                mtk::EntityRank aOutputEntityRank,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_nodes_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_nodes_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_nodes_connected_to_face_loc_inds( moris_index aFaceIndex ) const;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const override;

        //-------------------------------------------------------------------------------

        void get_elements_in_support_of_basis(
                uint                aMeshIndex,
                uint                aBasisIndex,
                Matrix< IndexMat >& aElementIndices ) override;

        //-------------------------------------------------------------------------------

        void get_nodes_indices_in_bounding_box(
                const moris::Matrix< DDRMat >& aPoint,
                const moris::Matrix< DDRMat >& aBoundingBoxSize,
                moris::Matrix< IndexMat >&     aNodeIndices ) override;

        //-------------------------------------------------------------------------------

        luint
        get_num_active_bg_elements_on_discretization_mesh_index_including_aura( moris_index const aDiscretizationMeshIndex ) override
        {
            return mMesh->get_num_active_bg_elements_on_discretization_mesh_index_including_aura( aDiscretizationMeshIndex );
        }

        // ----------------------------------------------------------------------------

        void
        get_active_bg_element_indices_on_discretization_mesh_index_including_aura(
                moris_index const  aDiscretizationMeshIndex,
                Matrix< DDLUMat >& aElementIDs ) override
        {
            mMesh->get_active_bg_element_indices_on_discretization_mesh_index_including_aura( aDiscretizationMeshIndex, aElementIDs );
        }

        // ----------------------------------------------------------------------------

        void get_elements_in_bspline_element(
                moris_index           aBspElementIndex,
                moris_index           aDiscretizationMeshIndex,
                Vector< mtk::Cell* >& aCells ) override;

        // ----------------------------------------------------------------------------

        void
        get_lagrange_elements_in_bspline_elements(
                moris_index                      aDiscretizationMeshIndex,
                Vector< Vector< mtk::Cell* > >&  aCells,
                Vector< Vector< moris_index > >& aCellIndices,
                Vector< moris_index >&           aLagToBspCellIndices,
                Vector< uint >&                  aBspCellRefineLevels,
                Vector< mtk::Cell* >&            aBspCells ) override;

        // ----------------------------------------------------------------------------

        void
        get_extended_t_matrix(
                moris_index                       aDiscretizationMeshIndex,
                moris_index                       aBSplineCellIndex,
                moris::mtk::Cell&                 aLagrangeCell,
                Vector< Vector< mtk::Vertex* > >& tBsplineBasis,
                Vector< Matrix< DDRMat > >&       tWeights ) override;

        // ----------------------------------------------------------------------------

        void
        get_L2_projection_matrix(
                moris_index                             aDiscretizationMeshIndex,
                const mtk::Cell*                        aRootBSplineCell,
                const mtk::Cell*                        aExtendedBSplineCell,
                Vector< Vector< const mtk::Vertex* > >& tRootBsplineBasis,
                Vector< const mtk::Vertex* >&           tExtendedBsplineBasis,
                Vector< Matrix< DDRMat > >&             tWeights ) override;

        // ----------------------------------------------------------------------------

        const luint*
        get_bspline_element_ijk_level(
                moris_index      aDiscretizationMeshIndex,
                const mtk::Cell* aBsplineElement,
                uint&            aLevel ) override;

        // ----------------------------------------------------------------------------

        void get_elements_in_interpolation_cluster(
                moris_index           aElementIndex,
                moris_index           aDiscretizationMeshIndex,
                Vector< mtk::Cell* >& tCells ) override;

        //-------------------------------------------------------------------------------

        void
        get_elements_in_bspline_element_and_side_ordinal(
                moris_index           aBsplineElementIndex,
                moris_index           aDiscretizationMeshIndex,
                moris_index           aSideOrdinal,
                Vector< mtk::Cell* >& aCells ) override;

        //-------------------------------------------------------------------------------

        void get_elements_in_interpolation_cluster_and_side_ordinal(
                moris_index           aElementIndex,
                moris_index           aDiscretizationMeshIndex,
                moris_index           aSideOrdinal,
                Vector< mtk::Cell* >& aCells ) override;

        //-------------------------------------------------------------------------------

        uint get_num_basis_functions( uint aMeshIndex ) override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat >
        get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const override;

        //-------------------------------------------------------------------------------

        bool
        get_elements_connected_to_element_through_face_ord(
                moris_index            aBaseElementIndex,
                moris_index            aMySideOrdinal,
                moris_index&           aMyRefineLevel,
                Vector< moris_index >& aNeighborElements,
                Vector< moris_index >& aNeighborElementSideOrdinals,
                Vector< moris_index >& aTransitionLocations,
                Vector< moris_index >& aNeighborRefinementLevels ) const override;

        //-------------------------------------------------------------------------------
        //          Global ID Functions
        //-------------------------------------------------------------------------------

        moris_id get_glb_entity_id_from_entity_loc_index(
                moris_index     aEntityIndex,
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        moris_index get_loc_entity_ind_from_entity_glb_id(
                moris_id        aEntityId,
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        moris_id get_max_entity_id(
                mtk::EntityRank aEntityRank,
                moris_index     aIndex ) const override;

        //-------------------------------------------------------------------------------
        //          Coordinate Field Functions
        //-------------------------------------------------------------------------------

        Matrix< DDRMat > get_node_coordinate( moris_index aNodeIndex ) const override;

        //-------------------------------------------------------------------------------
        //           Entity Ownership Functions
        //-------------------------------------------------------------------------------

        /**
         * Gets the owner of a node.
         *
         * @param aNodeIndex Node index
         * @return Node owner
         */
        uint get_node_owner( moris_index aNodeIndex ) const override;

        /**
         * Gets the owner of an element.
         *
         * @param aElementIndex Element index
         * @return Element owner
         */
        uint get_element_owner( moris_index aElementIndex ) const override;

        uint get_entity_owner(
                moris_index     aEntityIndex,
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) const override;

        // FIXME Needs parallel implementation
        void
        get_processors_whom_share_entity(
                moris_index      aEntityIndex,
                mtk::EntityRank  aEntityRank,
                Matrix< IdMat >& aProcsWhomShareEntity ) const override;

        mtk::EntityRank
        get_facet_rank() const override;

        //-------------------------------------------------------------------------------
        //           Set Functions
        //-------------------------------------------------------------------------------

        void get_sideset_elems_loc_inds_and_ords(
                const std::string&  aSetName,
                Matrix< IndexMat >& aElemIndices,
                Matrix< IndexMat >& aSidesetOrdinals ) const override;

        //-------------------------------------------------------------------------------

        Vector< std::string > get_set_names( mtk::EntityRank aSetEntityRank ) const override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_set_entity_loc_inds(
                mtk::EntityRank    aSetEntityRank,
                const std::string& aSetName ) const override;

        //-------------------------------------------------------------------------------
        //           Pointer Functions for FEM
        //-------------------------------------------------------------------------------

        mtk::Vertex&
        get_mtk_vertex( moris_index aVertexIndex ) override
        {
            return *mMesh->get_node_by_index( aVertexIndex );
        }

        //-------------------------------------------------------------------------------

        mtk::Vertex const &
        get_mtk_vertex( moris_index aVertexIndex ) const override
        {
            return *mMesh->get_node_by_index( aVertexIndex );
        }

        //-------------------------------------------------------------------------------

        mtk::Cell& get_mtk_cell( moris_index aElementIndex ) override;

        //-------------------------------------------------------------------------------

        mtk::Cell const & get_mtk_cell( moris_index aElementIndex ) const override;

        //-------------------------------------------------------------------------------

        mtk::Cell& get_writable_mtk_cell( moris_index aElementIndex ) override;

        //-------------------------------------------------------------------------------

        moris::mtk::Facet*
        get_facet( moris_index aFacetIndex ) override
        {
            return mMesh->get_facet( aFacetIndex );
        }

        //-------------------------------------------------------------------------------

        Vector< moris::mtk::Vertex const * >
        get_all_vertices() const override;

        //-------------------------------------------------------------------------------

        //            Vector<moris::mtk::Vertex const *>
        //            get_all_vertices_including_aura() const;

        //-------------------------------------------------------------------------------

        void get_adof_map(
                uint                          aBSplineIndex,
                map< moris_id, moris_index >& aAdofMap ) const override;

        //-------------------------------------------------------------------------------

        moris_index get_field_ind(
                const std::string& aFieldLabel,
                mtk::EntityRank    aEntityRank ) const override;

        //-------------------------------------------------------------------------------

        uint get_num_fields(
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        real& get_value_of_scalar_field(
                moris_index     aFieldIndex,
                mtk::EntityRank aEntityRank,
                uint            aEntityIndex,
                moris_index     aIndex = 0 ) override;

        //-------------------------------------------------------------------------------

        const real& get_value_of_scalar_field(
                moris_index     aFieldIndex,
                mtk::EntityRank aEntityRank,
                uint            aEntityIndex,
                moris_index     aIndex = 0 ) const override;

        //-------------------------------------------------------------------------------

        Matrix< DDRMat >& get_field(
                moris_index     aFieldIndex,
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) override;

        //-------------------------------------------------------------------------------

        mtk::CellTopology get_blockset_topology( const std::string& aSetName ) override;

        //-------------------------------------------------------------------------------

        mtk::CellShape get_IG_blockset_shape( const std::string& aSetName ) override;

        //-------------------------------------------------------------------------------

        mtk::CellShape get_IP_blockset_shape( const std::string& aSetName ) override;

        //-------------------------------------------------------------------------------

      private:
        //-------------------------------------------------------------------------------

        void get_element_indices_from_memory_indices(
                const Matrix< DDLUMat >& aMemoryIndices,
                Matrix< IndexMat >&      aIndices ) const;

        //-------------------------------------------------------------------------------
        /**
         * subroutine for get_elements_connected_to_element_loc_inds(
         */
        void collect_memory_indices_of_active_element_neighbors(
                moris_index        aElementIndex,
                Matrix< DDLUMat >& aMemoryIndices,
                Matrix< DDLUMat >& aThisCellFacetOrds,
                Matrix< DDLUMat >& aNeighborCellFacetOrds,
                Matrix< DDLUMat >& aTransitionNeighborCellLocation,
                luint&             aCounter ) const;

        //-------------------------------------------------------------------------------

        /**
         * subroutine for get_elements_connected_to_node_loc_inds
         */
        void collect_memory_indices_of_active_elements_connected_to_node(
                moris_index        aNodeIndex,
                Matrix< DDLUMat >& aMemoryIndices ) const;

      public:
        /**
         * Get the T-matrix of a node.
         *
         * @param aNodeIndex Node index
         * @param aBSplineMeshIndex B-spline mesh index
         * @return T-matrix
         */
        const Matrix< DDRMat >& get_t_matrix_of_node_loc_ind(
                uint aNodeIndex,
                uint aBSplineMeshIndex ) override;

      private:
        //-------------------------------------------------------------------------------

        uint get_level_of_entity_loc_ind(
                mtk::EntityRank aEntityRank,
                uint            aEntityIndex,
                moris_index     aIndex = 0 ) override;

        //-------------------------------------------------------------------------------

        uint get_max_level_of_entity(
                mtk::EntityRank aEntityRank,
                moris_index     aIndex = 0 ) override;

        //-------------------------------------------------------------------------------

        Matrix< IndexMat > get_inds_of_active_elements_connected_to_basis( const Basis* aBasis ) const;

        //-------------------------------------------------------------------------------

        void setup_glb_to_local_maps();

        //-------------------------------------------------------------------------------

        void setup_entity_global_to_local_map(
                mtk::EntityRank aEntityRank,
                uint&           aCounter,
                moris_index     aIndex = 0 );

        //------------------------------------------------------------------------------
        // ##############################################
        // Multi-grid accessor functions
        // ##############################################
        //-------------------------------------------------------------------------------

        uint
        get_num_interpolations() override
        {
            return mMesh->get_number_of_bspline_meshes();
        };

        //-------------------------------------------------------------------------------

        uint
        get_max_level( const moris_index aInterpolationIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )->get_max_level();
        };

        //-------------------------------------------------------------------------------

        uint
        get_num_basis( const moris_index aInterpolationIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )->get_number_of_indexed_basis();
        }

        //-------------------------------------------------------------------------------

        uint
        get_basis_level(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )
                    ->get_basis_by_index( aBasisIndex )
                    ->get_level();
        }

        //-------------------------------------------------------------------------------

        uint
        get_num_coarse_basis_of_basis(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )
                    ->get_basis_by_index( aBasisIndex )
                    ->get_number_of_parents();
        }

        //-------------------------------------------------------------------------------

        uint
        get_coarse_basis_index_of_basis(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex,
                moris_index aCoarseParentIndexForBasis ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )
                    ->get_basis_by_index( aBasisIndex )
                    ->get_parent( aCoarseParentIndexForBasis )
                    ->get_index();
        }

        //-------------------------------------------------------------------------------

        moris::Matrix< DDSMat >
        get_fine_basis_inds_of_basis(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )
                    ->get_children_ind_for_basis( aBasisIndex );
        }

        //-------------------------------------------------------------------------------

        moris::Matrix< DDRMat >
        get_fine_basis_weights_of_basis(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )
                    ->get_children_weights_for_parent( aBasisIndex );
        }

        //-------------------------------------------------------------------------------

#ifdef MORIS_HAVE_DEBUG
        Matrix< DDRMat >
        get_basis_coords(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex ) override
        {
            return mMesh->get_bspline_mesh( aInterpolationIndex )
                    ->get_basis_by_index( aBasisIndex )
                    ->get_coords();
        }

        //-------------------------------------------------------------------------------

        sint
        get_basis_status(
                moris_index aInterpolationIndex,
                moris_index aBasisIndex ) override
        {
            sint tStatus = -1;

            if ( mMesh->get_bspline_mesh( aInterpolationIndex )
                            ->get_basis_by_index( aBasisIndex )
                            ->is_active() )
            {
                tStatus = 1;
            }
            else if ( mMesh->get_bspline_mesh( aInterpolationIndex )
                              ->get_basis_by_index( aBasisIndex )
                              ->is_refined() )
            {
                tStatus = 0;
            }
            return tStatus;
        }
#endif

        //-------------------------------------------------------------------------------

    };    // end class: hmr::Mesh

    //-------------------------------------------------------------------------------

}    // namespace moris::hmr
