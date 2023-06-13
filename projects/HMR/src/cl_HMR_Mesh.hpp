/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Mesh.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_INTERFACE_HPP_
#define SRC_HMR_CL_HMR_INTERFACE_HPP_

#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "cl_Cell.hpp"    //CNT/src
#include "cl_Mesh_Enums.hpp"
#include "MTK_Tools.hpp"
#include "cl_MTK_Mesh_Core.hpp"    //MTK/src

namespace moris
{
    namespace hmr
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

            moris::Cell< moris::hmr::BSpline_Mesh_Base* > mDummyBSplineMeshes;

            // give access to the member data to stk based interpolation and integration meshes
            friend class Interpolation_Mesh_HMR;
            friend class Integration_Mesh_HMR;

            // Global to Local Entity Map
            // moris::Cell(0) - Node Global to Local
            // moris::Cell(1) - Edge Global to Local
            // moris::Cell(2) - Face Global to Local
            // moris::Cell(3) - Cell Global to Local
            // Keep in mind not all of these are created

            moris::Cell< std::unordered_map< moris_id, moris_index > > mEntityGlobalToLocalMap;

            //-------------------------------------------------------------------------------

          public:
            //-------------------------------------------------------------------------------

            /**
             * mesh constructor, to be called from HMR
             */
            Mesh(
                    std::shared_ptr< Database > aDatabase,
                    uint                 aLagrangeOrder,
                    uint                 aLagrangePattern );

            Mesh(
                    std::shared_ptr< Database > aDatabase,
                    uint                 aOrder,
                    uint                 aLagrangePattern,
                    uint                 aBsplinePattern );

            Mesh(
                    std::shared_ptr< Database > aDatabase,
                    uint                 aLagrangeOrder,
                    uint                 aLagrangePattern,
                    uint                 aBSplineOrder,
                    uint                 aBSplinePattern );

            Mesh(
                    std::shared_ptr< Database > aDatabase,
                    uint                 aLagrangeMeshIndex );

            //-------------------------------------------------------------------------------

            /**
             * destructor
             */
            virtual ~Mesh();

            //-------------------------------------------------------------------------------

            /**
             * return the type of this mesh
             */
            MeshType
            get_mesh_type() const
            {
                return MeshType::HMR;
            }

            //-------------------------------------------------------------------------------

            /**
             * provides a moris::Matrix< DDUMat > containing the IDs this mesh has
             * to communicate with
             */
            Matrix< IdMat > get_communication_table() const;

            // -----------------------------------------------------------------------------

            /**
             * returns the proc neighbors for this proc. This function return 9 entries for 2D and 27 for 3D.
             * Some of them can be gNoId if this proc does not exist.
             */
            Matrix< IdMat > get_proc_neighbors() const;

            //-------------------------------------------------------------------------------

            /**
             * creates a new field pointer that is linked to this mesh
             */
            std::shared_ptr< Field > create_field(
                    const std::string& aLabel,
                    uint        aBSplineOrder );

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
            uint get_spatial_dim() const;

            //-------------------------------------------------------------------------------

            uint get_num_entities( const enum EntityRank aEntityRank,
                    const moris_index                    aIndex = 0 ) const;

            //-------------------------------------------------------------------------------

            uint get_num_elements_including_aura() const;

            //-------------------------------------------------------------------------------

            uint get_num_nodes() const;

            //-------------------------------------------------------------------------------

            uint get_num_edges() const;

            //-------------------------------------------------------------------------------

            uint get_num_faces() const;

            //-------------------------------------------------------------------------------

            uint get_num_elems() const;

            /**
             * Gets element indices in a block set.
             *
             * @param aSetIndex Block set index
             * @return Element indices in the set
             */
            Matrix< IndexMat > get_element_indices_in_block_set( uint aSetIndex );

            /**
             * Gets the element IDs in a block set.
             *
             * @param aSetIndex Block set index
             * @return Element IDs in the set
             */
            Matrix< IdMat > get_element_ids_in_block_set( uint aSetIndex );

            //-------------------------------------------------------------------------------

            uint get_max_num_coeffs_on_proc( const uint aBSplineMeshIndex ) const;

            /**
             * Gets the indices of the B-splines which form the basis of the given node.
             *
             * @param aNodeIndex Node index
             * @param aBSplineMeshIndex B-spline mesh index
             * @return B-spline indices
             */
            Matrix< IndexMat > get_coefficient_indices_of_node(
                    uint aNodeIndex,
                    uint aBSplineMeshIndex );

            /**
             * Gets the IDs of the B-splines which form the basis of a given node.
             *
             * @param aNodeIndex Node index
             * @param aBSplineMeshIndex B-spline mesh index
             * @return B-spline IDs
             */
            Matrix< IdMat > get_coefficient_IDs_of_node(
                    uint aNodeIndex,
                    uint aBSplineMeshIndex );

            Matrix< IdMat > get_coefficient_owners_of_node(
                    uint aNodeIndex,
                    uint aBSplineMeshIndex );

            Matrix< IdMat > get_coefficient_ijkl_IDs_of_node(
                    uint aNodeIndex,
                    uint aBSplineMeshIndex );

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_entity_connected_to_entity_loc_inds(
                    moris_index       aEntityIndex,
                    enum EntityRank   aInputEntityRank,
                    enum EntityRank   aOutputEntityRank,
                    const moris_index aIndex = 0 ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_entity_connected_to_entity_glob_ids(
                    moris_index       aEntityId,
                    enum EntityRank   aInputEntityRank,
                    enum EntityRank   aOutputEntityRank,
                    const moris_index aIndex = 0 ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_edge_loc_inds( moris_index aEdgeIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_face_loc_inds( moris_index aFaceIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_edges_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_edges_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_faces_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_faces_connected_to_element_loc_inds( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_elements_connected_to_node_loc_inds( moris_index aNodeIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_elements_connected_to_face_loc_inds( moris_index aFaceIndex ) const;

            //-------------------------------------------------------------------------------

            void get_elements_in_support_of_basis( const uint aMeshIndex,
                    const uint                                aBasisIndex,
                    Matrix< IndexMat >&                       aElementIndices );

            //-------------------------------------------------------------------------------

            void get_nodes_indices_in_bounding_box(
                    const moris::Matrix< DDRMat >& aPoint,
                    const moris::Matrix< DDRMat >& aBoundingBoxSize,
                    moris::Matrix< IndexMat >&     aNodeIndices );

            //-------------------------------------------------------------------------------

            luint
            get_num_active_bg_elements_on_discretization_mesh_index_including_aura( moris_index const aDiscretizationMeshIndex )
            {
                return mMesh->get_num_active_bg_elements_on_discretization_mesh_index_including_aura( aDiscretizationMeshIndex );
            }

            // ----------------------------------------------------------------------------

            void
            get_active_bg_element_indices_on_discretization_mesh_index_including_aura(
                    moris_index const  aDiscretizationMeshIndex,
                    Matrix< DDLUMat >& aElementIDs )
            {
                mMesh->get_active_bg_element_indices_on_discretization_mesh_index_including_aura( aDiscretizationMeshIndex, aElementIDs );
            }

            // ----------------------------------------------------------------------------

            void get_elements_in_bspline_element(
                    moris_index const          aBspElementIndex,
                    moris_index const          aDiscretizationMeshIndex,
                    moris::Cell< mtk::Cell* >& aCells );

            // ----------------------------------------------------------------------------

            void
            get_lagrange_elements_in_bspline_elements(
                    moris_index const                          aDiscretizationMeshIndex,
                    moris::Cell< moris::Cell< mtk::Cell* > >&  aCells,
                    moris::Cell< moris::Cell< moris_index > >& aCellIndices,
                    moris::Cell< moris_index >&                aLagToBspCellIndices,
                    moris::Cell< uint >&                       aBspCellRefineLevels );

        // ----------------------------------------------------------------------------

        virtual void
        get_extended_t_matrix(
                moris_index const &                  aDiscretizationMeshIndex,
                moris_index const &                  aBSplineCellIndex,
                moris::mtk::Cell&                    aLagrangeCell,
                moris::Cell< Cell< mtk::Vertex* > >& tBsplineBasis,
                moris::Cell< Matrix< DDRMat > >&     tWeights ) override;

        // ----------------------------------------------------------------------------

        virtual void
        get_L2_projection_matrix(
                moris_index const &                         aDiscretizationMeshIndex,
                moris_index const &                         aRootBSplineCellIndex,
                moris_index const &                         aExtendedBSplineCellIndex,
                moris::Cell< moris::Cell< mtk::Vertex* > >& tRootBsplineBasis,
                moris::Cell< mtk::Vertex* >&                tExtendedBsplineBasis,
                moris::Cell< Matrix< DDRMat > >&            tWeights ) override;

        // ----------------------------------------------------------------------------

        virtual const luint*
        get_bspline_element_ijk_level(
                moris_index const & aDiscretizationMeshIndex,
                moris_index const & aBsplineElementIndex,
                uint                aLevel ) override;

        // ----------------------------------------------------------------------------

        void get_elements_in_interpolation_cluster(
                moris_index                aElementIndex,
                moris_index                aDiscretizationMeshIndex,
                moris::Cell< mtk::Cell* >& tCells );

        //-------------------------------------------------------------------------------

        void
        get_elements_in_bspline_element_and_side_ordinal(
                moris_index const          aBsplineElementIndex,
                moris_index const          aDiscretizationMeshIndex,
                moris_index const          aSideOrdinal,
                moris::Cell< mtk::Cell* >& aCells );

        //-------------------------------------------------------------------------------

        void get_elements_in_interpolation_cluster_and_side_ordinal(
                moris_index const          aElementIndex,
                moris_index const          aDiscretizationMeshIndex,
                moris_index const          aSideOrdinal,
                moris::Cell< mtk::Cell* >& aCells );

            //-------------------------------------------------------------------------------

            uint get_num_basis_functions( const uint aMeshIndex );

            //-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat >
            get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

            bool
            get_elements_connected_to_element_through_face_ord(
                    moris_index                 aBaseElementIndex,
                    moris_index                 aMySideOrdinal,
                    moris_index&                aMyRefineLevel,
                    moris::Cell< moris_index >& aNeighborElements,
                    moris::Cell< moris_index >& aNeighborElementSideOrdinals,
                    moris::Cell< moris_index >& aTransitionLocations,
                    moris::Cell< moris_index >& aNeighborRefinementLevels ) const;

            //-------------------------------------------------------------------------------
            //          Global ID Functions
            //-------------------------------------------------------------------------------

            moris_id get_glb_entity_id_from_entity_loc_index(
                    moris_index       aEntityIndex,
                    enum EntityRank   aEntityRank,
                    const moris_index aIndex = 0 ) const;

            moris_id get_glb_element_id_from_element_loc_index( moris_index aEntityIndex ) const;

            //-------------------------------------------------------------------------------
            moris_index get_loc_entity_ind_from_entity_glb_id(
                    moris_id          aEntityId,
                    enum EntityRank   aEntityRank,
                    const moris_index aIndex = 0 ) const;

            //-------------------------------------------------------------------------------

            moris_id get_max_entity_id( enum EntityRank aEntityRank,
                    const moris_index                   aIndex ) const;

            //-------------------------------------------------------------------------------
            //          Coordinate Field Functions
            //-------------------------------------------------------------------------------

            Matrix< DDRMat > get_node_coordinate( moris_index aNodeIndex ) const;

            //-------------------------------------------------------------------------------
            //           Entity Ownership Functions
            //-------------------------------------------------------------------------------

            /**
             * Gets the owner of a node.
             *
             * @param aNodeIndex Node index
             * @return Node owner
             */
            uint get_node_owner( moris_index aNodeIndex ) const;

            /**
             * Gets the owner of an element.
             *
             * @param aElementIndex Element index
             * @return Element owner
             */
            uint get_element_owner( moris_index aElementIndex ) const;

            uint get_entity_owner(
                    moris_index       aEntityIndex,
                    enum EntityRank   aEntityRank,
                    const moris_index aIndex = 0 ) const;

            // FIXME Needs parallel implementation
            void
            get_processors_whom_share_entity(
                    moris_index      aEntityIndex,
                    enum EntityRank  aEntityRank,
                    Matrix< IdMat >& aProcsWhomShareEntity ) const;

            enum EntityRank
            get_facet_rank() const;

            //-------------------------------------------------------------------------------
            //           Set Functions
            //-------------------------------------------------------------------------------

            void get_sideset_elems_loc_inds_and_ords(
                    const std::string&  aSetName,
                    Matrix< IndexMat >& aElemIndices,
                    Matrix< IndexMat >& aSidesetOrdinals ) const;

            //-------------------------------------------------------------------------------

            moris::Cell< std::string > get_set_names( enum EntityRank aSetEntityRank ) const;

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_set_entity_loc_inds(
                    enum EntityRank aSetEntityRank,
                    std::string     aSetName ) const;

            //-------------------------------------------------------------------------------
            //           Pointer Functions for FEM
            //-------------------------------------------------------------------------------

            mtk::Vertex&
            get_mtk_vertex( moris_index aVertexIndex )
            {
                return *mMesh->get_node_by_index( aVertexIndex );
            }

            //-------------------------------------------------------------------------------

            mtk::Vertex const &
            get_mtk_vertex( moris_index aVertexIndex ) const
            {
                return *mMesh->get_node_by_index( aVertexIndex );
            }

            //-------------------------------------------------------------------------------

            //            mtk::Vertex & get_mtk_vertex_including_aura( moris_index aVertexIndex )
            //            {
            //                return *mMesh->get_node_by_index_including_aura( aVertexIndex );
            //            }

            //-------------------------------------------------------------------------------

            //            mtk::Vertex const & get_mtk_vertex_including_aura( moris_index aVertexIndex ) const
            //            {
            //                return *mMesh->get_node_by_index_including_aura( aVertexIndex );
            //            }

            //-------------------------------------------------------------------------------

            mtk::Cell& get_mtk_cell( moris_index aElementIndex );

            //-------------------------------------------------------------------------------

            mtk::Cell const & get_mtk_cell( moris_index aElementIndex ) const;

            //-------------------------------------------------------------------------------

            mtk::Cell& get_writable_mtk_cell( moris_index aElementIndex );
            //-------------------------------------------------------------------------------

            moris::mtk::Facet*
            get_facet( moris_index aFacetIndex )
            {
                return mMesh->get_facet( aFacetIndex );
            }

            //-------------------------------------------------------------------------------

            moris::Cell< moris::mtk::Vertex const * >
            get_all_vertices() const;

            //-------------------------------------------------------------------------------

            //            moris::Cell<moris::mtk::Vertex const *>
            //            get_all_vertices_including_aura() const;

            //-------------------------------------------------------------------------------

            void get_adof_map(
                    const uint                    aBSplineIndex,
                    map< moris_id, moris_index >& aAdofMap ) const;

            //-------------------------------------------------------------------------------

            moris_index get_field_ind(
                    const std::string&    aFieldLabel,
                    const enum EntityRank aEntityRank ) const;

            //-------------------------------------------------------------------------------

            uint get_num_fields(
                    const enum EntityRank aEntityRank,
                    const moris_index     aIndex = 0 ) const;

            //-------------------------------------------------------------------------------

            real& get_value_of_scalar_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex,
                    const moris_index     aIndex = 0 );

            //-------------------------------------------------------------------------------

            const real& get_value_of_scalar_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex,
                    const moris_index     aIndex = 0 ) const;

            //-------------------------------------------------------------------------------

            Matrix< DDRMat >& get_field(
                    const moris_index     aFieldIndex,
                    const enum EntityRank aEntityRank,
                    const moris_index     aIndex = 0 );

            //-------------------------------------------------------------------------------

            enum CellTopology get_blockset_topology( const std::string& aSetName );

            //-------------------------------------------------------------------------------

            enum CellShape get_IG_blockset_shape( const std::string& aSetName );

            //-------------------------------------------------------------------------------

            enum CellShape get_IP_blockset_shape( const std::string& aSetName );

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
                    const moris_index  aElementIndex,
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
                    const moris_index  aNodeIndex,
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
                    uint aBSplineMeshIndex );

          private:
            //-------------------------------------------------------------------------------

            uint get_level_of_entity_loc_ind(
                    const enum EntityRank aEntityRank,
                    const uint            aEntityIndex,
                    const moris_index     aIndex = 0 );

            //-------------------------------------------------------------------------------

            uint get_max_level_of_entity( const enum EntityRank aEntityRank, const moris_index aIndex = 0 );

            //-------------------------------------------------------------------------------

            Matrix< IndexMat > get_inds_of_active_elements_connected_to_basis( const Basis* aBasis ) const;

            //-------------------------------------------------------------------------------

            void setup_glb_to_local_maps();

            //-------------------------------------------------------------------------------

            void setup_entity_global_to_local_map(
                    enum EntityRank   aEntityRank,
                    uint&             aCounter,
                    const moris_index aIndex = 0 );

            //------------------------------------------------------------------------------
            // ##############################################
            // Multi-grid accessor functions
            // ##############################################
            //-------------------------------------------------------------------------------

            uint
            get_num_interpolations()
            {
                return mMesh->get_number_of_bspline_meshes();
            };

            //-------------------------------------------------------------------------------

            uint
            get_max_level( const moris_index aInterpolationIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )->get_max_level();
            };

            //-------------------------------------------------------------------------------

            uint
            get_num_basis( const moris_index aInterpolationIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )->get_number_of_indexed_basis();
            }

            //-------------------------------------------------------------------------------

            uint
            get_basis_level(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )
                        ->get_basis_by_index( aBasisIndex )
                        ->get_level();
            }

            //-------------------------------------------------------------------------------

            uint
            get_num_coarse_basis_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )
                        ->get_basis_by_index( aBasisIndex )
                        ->get_number_of_parents();
            }

            //-------------------------------------------------------------------------------

            uint
            get_coarse_basis_index_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex,
                    const moris_index aCoarseParentIndexForBasis )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )
                        ->get_basis_by_index( aBasisIndex )
                        ->get_parent( aCoarseParentIndexForBasis )
                        ->get_index();
            }

            //-------------------------------------------------------------------------------

            moris::Matrix< DDSMat >
            get_fine_basis_inds_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )
                        ->get_children_ind_for_basis( aBasisIndex );
            }

            //-------------------------------------------------------------------------------

            moris::Matrix< DDRMat >
            get_fine_basis_weights_of_basis(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )
                        ->get_children_weights_for_parent( aBasisIndex );
            }

            //-------------------------------------------------------------------------------

#ifdef MORIS_HAVE_DEBUG
            Matrix< DDRMat >
            get_basis_coords( const moris_index aInterpolationIndex,
                    const moris_index           aBasisIndex )
            {
                return mMesh->get_bspline_mesh( aInterpolationIndex )
                        ->get_basis_by_index( aBasisIndex )
                        ->get_coords();
            }

            //-------------------------------------------------------------------------------

            sint
            get_basis_status(
                    const moris_index aInterpolationIndex,
                    const moris_index aBasisIndex )
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
        };

    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_INTERFACE_HPP_ */
