/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Enriched_Interpolation_Mesh.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_

#include "cl_Vector.hpp"

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_XTK_Interpolation_Vertex_Unzipped.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Vertex_Enrichment.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_Parameter_List.hpp"

using namespace moris;


namespace moris::xtk
{

    class Model;
    class Field;
    class Enriched_Interpolation_Mesh : public mtk::Interpolation_Mesh
    {
        // friend classes
        friend class Enrichment;
        friend class Enriched_Integration_Mesh;
        friend class Ghost_Stabilization;

      protected:
        // Model pointer
        Model* mXTKModel;

        // basis rank
        mtk::EntityRank mBasisRank = mtk::EntityRank::INVALID;

        // mesh index
        Matrix< IndexMat >                             mMeshIndices;
        Matrix< IndexMat >                             mUnenrichedMeshIndices;      // lis of the meshes that will be unenriched by overriding t-matrices
        std::unordered_map< moris_index, moris_index > mMeshIndexToLocMeshIndex;    // over allocated

        // enriched interpolation vertices
        moris::uint                              mNumVerts;
        Vector< Interpolation_Vertex_Unzipped* > mEnrichedInterpVerts;    // over allocated
        Vector< moris_index >                    mOwnedUnzippedVertices;
        Vector< moris_index >                    mNotOwnedUnzippedVertices;

        // enriched interpolation cells
        moris::uint                            mNumVertsPerInterpCell;
        Vector< Interpolation_Cell_Unzipped* > mEnrichedInterpCells;    // over allocated
        Vector< moris_index >                  mOwnedEnrichedInterpCells;
        Vector< moris_index >                  mNotOwnedEnrichedInterpCells;

        // for each outer cell (base interpolation vertex), indices of enriched vertices
        Vector< Vector< Vector< moris_index > > > mBaseInterpVertToVertEnrichmentIndex;
        // input: Discretization mesh index (DMI), BG mesh vertex index || output: list of enriched vertex indices

        // vertex enrichments or t-matrix
        // outer cell - mesh index (i.e. linear or quadratic b-spline enrichment)
        // inner cell
        Vector< Vector< Vertex_Enrichment* > > mInterpVertEnrichment;    // input: DMI || output: list of T-Matrices ("Vertex Enrichments")

        // vertex enrichment to parent vertex index (these are enriched interpolation vertex indices)
        Vector< Vector< moris_index > > mVertexEnrichmentParentVertexIndex;    // input: DMI, Vertex Enrichment index || output: index parent node on BG mesh

        // Bulk phase of each interpolation vertex
        Matrix< IndexMat > mVertexBulkPhase;      // input: index of unzipped vertex || output: bulk-phase index
        Matrix< IndexMat > mVertexMaxSubphase;    // input: index of unzipped vertex || output: maximum subphase ID on vertex // TODO: why do we need this?

        // basis coefficient to enriched basis coefficient
        Vector< Vector< Matrix< IndexMat > > > mCoeffToEnrichCoeffs;
        // input: DMI, index of non-enriched BF coefficient || ouput: list of enriched BF indices associated with it

        // local to global enriched basis vector
        Vector< Matrix< IdMat > >                             mEnrichCoeffLocToGlob;      // input: DMI, enriched BF index || output: global ID of that enriched BF
        Vector< std::unordered_map< moris_id, moris_index > > mGlobalToLocalBasisMaps;    // input: DMI || output: map ordered by global BF IDs with corresponding local BF index

        // basis ownership
        Vector< Matrix< IdMat > > mEnrichCoeffOwnership;    // input: DMI, enriched BF index || output: ?

        // basis bulk phase
        Vector< Matrix< IdMat > > mEnrichCoeffBulkPhase;    // input: DMI, enriched BF index || output: bulk phase the BF interpolates into

        // Entity maps
        Vector< Matrix< IdMat > >                             mLocalToGlobalMaps;    // TODO input: DMI || output:
        Vector< std::unordered_map< moris_id, moris_index > > mGlobalToLocalMaps;    // TODO input: DMI || output:

        // base interpolation cells to their enriched interpolation cells
        Vector< Vector< Interpolation_Cell_Unzipped* > > mBaseCellToEnrichedCell;    // TODO input: DMI || output:

        // a connectivity pointer that all the enriched interpolation cells use
        std::shared_ptr< moris::mtk::Cell_Info > mCellInfo;

        // not owned basis functions
        Vector< moris_index > mNotOwnedBasis;    // TODO input:  || output:
        Vector< moris_index > mOwnedBasis;       // TODO input:  || output:

        // block set information in an implicit form
        Vector< std::string >           mBlockSetNames;
        Vector< Vector< moris_index > > mElementIndicesInBlock;

        // Fields
        Vector< xtk::Field >                                     mFields;               /*Structure Node (0), Cell(1)*/
        Vector< std::unordered_map< std::string, moris_index > > mFieldLabelToIndex;    // input: Field rank (0: Node, 1:Element) , field label || output: field index

      public:
        Enriched_Interpolation_Mesh( Model* aXTKModel );

        ~Enriched_Interpolation_Mesh();

        //------------------------------------------------------------------------------
        // MTK Mesh Core Impl (see base class mtk::Mesh for function descriptions)
        //------------------------------------------------------------------------------
        mtk::MeshType                               get_mesh_type() const;
        moris::uint                                 get_spatial_dim() const;
        uint                                        get_num_entities( mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
        uint                                        get_max_num_coeffs_on_proc( const uint aBSplineMeshIndex ) const;
        Matrix< IndexMat >                          get_entity_connected_to_entity_loc_inds( moris_index aEntityIndex, mtk::EntityRank aInputEntityRank, mtk::EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const;
        Matrix< IndexMat >                          get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const;
        moris_id                                    get_glb_entity_id_from_entity_loc_index( moris_index aEntityIndex, mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
        moris_index                                 get_loc_entity_ind_from_entity_glb_id( moris_id aEntityId, mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
        const mtk::Interpolation_Mesh&              get_background_mesh() override;
        std::unordered_map< moris_id, moris_index > get_vertex_glb_id_to_loc_vertex_ind_map() const;
        Vector< mtk::Vertex const * >               get_all_vertices() const;
        Matrix< IdMat >                             get_entity_connected_to_entity_glob_ids( moris_id aEntityId, mtk::EntityRank aInputEntityRank, mtk::EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const;
        Matrix< DDRMat >                            get_node_coordinate( moris_index aNodeIndex ) const;
        Matrix< DDRMat >                            get_base_node_coordinate( moris_index aBaseNodeIndex ) const;
        mtk::Vertex&                                get_mtk_vertex( moris_index aVertexIndex );
        mtk::Vertex const &                         get_mtk_vertex( moris_index aVertexIndex ) const;
        mtk::Cell const &                           get_mtk_cell( moris_index aElementIndex ) const;
        mtk::Cell&                                  get_mtk_cell( moris_index aElementIndex );
        mtk::Cell&                                  get_writable_mtk_cell( moris_index aElementIndex );
        Matrix< IdMat >                             get_communication_table() const;
        uint                                        get_num_elements();
        Matrix< IndexMat >                          get_element_indices_in_block_set( uint aSetIndex );
        moris_id                                    get_max_entity_id( mtk::EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
        void                                        get_adof_map( const uint aBSplineIndex, map< moris_id, moris_index >& aAdofMap ) const;
        mtk::CellTopology                           get_blockset_topology( const std::string& aSetName );
        mtk::CellShape                              get_IG_blockset_shape( const std::string& aSetName );
        mtk::CellShape                              get_IP_blockset_shape( const std::string& aSetName );
        uint                                        get_node_owner( moris_index aNodeIndex ) const;
        uint                                        get_element_owner( moris_index aElementIndex ) const;

        //------------------------------------------------------------------------------
        // end mesh core functions
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // multi-grid accessor functions
        //------------------------------------------------------------------------------

        uint             get_num_interpolations();
        uint             get_max_level( const moris_index aInterpolationIndex );
        uint             get_num_basis( const moris_index aInterpolationIndex );
        uint             get_basis_level( const moris_index aInterpolationIndex, const moris_index aBasisIndex );
        uint             get_num_coarse_basis_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex );
        uint             get_coarse_basis_index_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex, const moris_index aCoarseParentIndex );
        Matrix< DDSMat > get_fine_basis_inds_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex );
        Matrix< DDRMat > get_fine_basis_weights_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex );

        //------------------------------------------------------------------------------
        // Accessor functions of XTK specific data structures
        //------------------------------------------------------------------------------
        /*!
         * Get the enriched interpolation coefficients associated with a background coefficient
         */
        Matrix< IndexMat > const &
        get_enriched_coefficients_at_background_coefficient(
                moris_index const & aMeshIndex,
                moris_index         aBackgroundCoeffIndex ) const;
        //------------------------------------------------------------------------------

        /*!
         * get the enriched coefficients at all the background mesh coefficients
         */
        Vector< Matrix< IndexMat > > const &
        get_enriched_coefficients_to_background_coefficients( moris_index const & aMeshIndex ) const;

        //------------------------------------------------------------------------------

        /*!
         * get the local enriched coefficient to global map
         */
        Matrix< IndexMat > const &
        get_enriched_coefficient_local_to_global_map( moris_index const & aMeshIndex ) const;

        //------------------------------------------------------------------------------

        /*!
         * Return the vector of background coefficient local to global
         */
        Matrix< IndexMat >
        get_background_coefficient_local_to_global_map() const;

        //------------------------------------------------------------------------------

        uint
        get_num_background_coefficients( moris_index const & aMeshIndex ) const;

        //------------------------------------------------------------------------------

        /*!
         * Returns the number of vertices per interpolation cell
         */
        uint
        get_num_verts_per_interp_cell();

        //------------------------------------------------------------------------------

        /*
         * Returns the interpolation vertex unzipped for provided vertex index
         */
        Interpolation_Vertex_Unzipped*
        get_unzipped_vertex_pointer( moris_index aVertexIndex );

        //------------------------------------------------------------------------------

        /*!
         * Return the enriched interpolation cells
         */
        Vector< Interpolation_Cell_Unzipped* > const &
        get_enriched_interpolation_cells() const;

        //------------------------------------------------------------------------------

        /*!
         * Get the number of interpolation (t-matrices) defined on this mesh
         */
        uint
        get_num_interpolation_types() const;

        //------------------------------------------------------------------------------

        /*
         * Get the interpolation index for a local index
         */
        moris_index
        get_interpolation_index( moris_index const & aLocalInterpIndex ) const;

        //------------------------------------------------------------------------------

        /*!
         * get basis owner
         */
        moris_index
        get_basis_owner(
                moris_index aBasisIndex,
                moris_index aMeshIndex );

        //------------------------------------------------------------------------------

        /*!
         * get basis bulk phase
         */
        moris_index
        get_basis_bulk_phase(
                moris_index const & aBasisIndex,
                moris_index const & aMeshIndex ) const;

      public:
        Vector< Interpolation_Cell_Unzipped* >&
        get_enriched_interpolation_cells();

        //------------------------------------------------------------------------------

        /**
         * @brief Return the owned interpolation cells, shared cells sorted by owning proc, owning procs
         *
         * @param aOwnedInterpCells
         * @param aNotOwnedInterpCells
         * @param aProcRanks
         */
        void
        get_owned_and_not_owned_enriched_interpolation_cells(
                Vector< Interpolation_Cell_Unzipped* >&           aOwnedInterpCells,
                Vector< Vector< Interpolation_Cell_Unzipped* > >& aNotOwnedInterpCells,
                Vector< uint >&                                   aProcRanks );

        //------------------------------------------------------------------------------

        Interpolation_Vertex_Unzipped&
        get_xtk_interp_vertex( moris::uint aVertexIndex );

        //------------------------------------------------------------------------------

        void
        add_proc_to_comm_table( moris_index aProcRank );

        //------------------------------------------------------------------------------

        /**
         * @brief Provided pointers to base interpolation cells, return all the enriched interpolation cells attached to these cells
         *
         * @param aBaseCells
         * @return Vector< Interpolation_Cell_Unzipped const* >
         */
        Vector< Interpolation_Cell_Unzipped const * >
        get_enriched_cells_from_base_cells( Vector< moris::mtk::Cell const * > const & aBaseCells ) const;

        //------------------------------------------------------------------------------

        /*!
         *  Single cell version of above
         */
        Vector< Interpolation_Cell_Unzipped const * >
        get_enriched_cells_from_base_cell( moris::mtk::Cell const * aBaseCells ) const;

        //------------------------------------------------------------------------------

        Interpolation_Vertex_Unzipped const &
        get_xtk_interp_vertex( moris::uint aVertexIndex ) const;

        //------------------------------------------------------------------------------

        moris_index
        get_enr_basis_index_from_enr_basis_id(
                moris_index const & aMeshIndex,
                moris_index const & aBasisId ) const;

        //------------------------------------------------------------------------------
        /**
         * @brief Get the enr basis id from enr basis index
         *
         * @param aMeshIndex Bspline mesh index
         * @param aBasisIndex index of the basis
         * @return moris_index Id of the basis
         */
        moris_index
        get_enr_basis_id_from_enr_basis_index(
                moris_index const & aMeshIndex,
                moris_index const & aBasisIndex ) const;
        //------------------------------------------------------------------------------

        /**
         * @brief Convert a entity indices to entity ids
         *
         * @param aIndices
         * @param aEntityRank
         * @return Matrix< IdMat >
         */
        Matrix< IdMat > convert_indices_to_ids(
                Matrix< IndexMat > const & aIndices,
                mtk::EntityRank            aEntityRank ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Convert a entity ids to entity indices
         *
         * @param aIds
         * @param aEntityRank
         * @return Matrix< IndexMat >
         */
        Matrix< IndexMat > convert_ids_to_indices(
                Matrix< IdMat > const & aIds,
                mtk::EntityRank         aEntityRank ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Return the mesh indices that will be enriched
         *
         * @return Matrix< IndexMat >
         */
        Matrix< IndexMat >
        get_enriched_mesh_indices() const;

        //------------------------------------------------------------------------------

        /**
         * @brief Set the enriched basis to bulkphase map from the outside
         *
         * @param aMeshListIndex
         * @param aBulkPhaseInEnrichedBasis
         */
        void
        set_enriched_basis_to_bulkphase_map(
                const moris_index aMeshIndex,
                Matrix< IdMat >   aBulkPhaseInEnrichedBasis );

        //------------------------------------------------------------------------------

        /**
         * @brief convert enriched //?
         *
         * @param aMeshIndex
         * @param aEnrichedIndices
         * @param aEnrichedIds
         */
        void
        convert_enriched_basis_indices_to_ids( moris_index const & aMeshIndex,
                Matrix< IndexMat > const &                         aEnrichedIndices,
                Matrix< IdMat >&                                   aEnrichedIds ) const;

        void
        convert_enriched_basis_indices_to_ids( moris_index const & aMeshIndex,
                Vector< Matrix< IndexMat > > const &               aEnrichedIndices,
                Vector< Matrix< IdMat > >&                         aEnrichedIds ) const;

        //------------------------------------------------------------------------------
        // Memory Map
        //------------------------------------------------------------------------------

        /*!
         * @brief Memory map of the enriched integration mesh
         * @return Memory map
         */
        moris::Memory_Map
        get_memory_usage();

        //------------------------------------------------------------------------------

        // Print functions
        void write_diagnostics();
        void print_enriched_cells( std::string aFile );
        void print_enriched_verts( std::string aFile );
        void print_enriched_verts_interpolation( const moris_index& aMeshIndex, std::string aFileName );
        void print_vertex_maps() const;
        void print_enriched_cell_maps() const;
        void print_basis_to_enriched_basis() const;
        void print_vertex_interpolation() const;
        void print_basis_information() const;

        // verification functions
        bool
        verify_basis_interpolating_into_cluster(
                mtk::Cluster const &       aCluster,
                moris_index const &        aMeshIndex,
                const mtk::Leader_Follower aIsLeader = mtk::Leader_Follower::LEADER );

        //------------------------------------------------------------------------------

        /**
         * @brief Checks whether a basis exists on a partition of the mesh
         *
         * @param aMeshIndex
         * @param aBasisId
         * @return true
         * @return false
         */
        bool
        basis_exists_on_partition( moris_index const & aMeshIndex,
                moris_index const &                    aBasisId );

        //------------------------------------------------------------------------------

        /**
         * @brief Add a basis function to the mesh. this is used by ghost to add basis functions in the aura returns the index
         *
         * @param aMeshIndex
         * @param aBasisIdToAdd
         * @param aBasisOwner
         * @param aBasisBulkPhase
         * @return moris_index
         */
        moris_index
        add_basis_function(
                moris_index const & aMeshIndex,
                moris_index const & aBasisIdToAdd,
                moris_index const & aBasisOwner,
                moris_index const & aBasisBulkPhase );

        //------------------------------------------------------------------------------

        /**
         * @brief Add a set of basis functions to the mesh. this is used by ghost to add basis functions in the aura returns the index
         *
         * @param aMeshIndex 
         * @param aBfIdsToAdd 
         * @param aBfOwners 
         * @param aBfBulkPhases 
         */
        void
        add_basis_functions(
                moris_index const & aMeshIndex,
                Vector< moris_id > const & aBfIdsToAdd,
                Vector< moris_id > const & aBfOwners,
                Vector< moris_index > const & aBfBulkPhases );


      protected:
        //------------------------------------------------------------------------------
        // functions used by enrichment for construction of the mesh
        //------------------------------------------------------------------------------

        /**
         * @brief Add a vertex enrichment to the member data. returns the index of the vertex enrichment.
         *
         * @param aMeshIndex
         * @param aBaseInterpVertex
         * @param aVertexEnrichment
         * @param aNewVertex
         * @return moris_index index of the vertex enrichment
         */
        moris_index
        add_vertex_enrichment( moris_index const & aMeshIndex,
                mtk::Vertex*                       aBaseInterpVertex,
                Vertex_Enrichment&                 aVertexEnrichment,
                bool&                              aNewVertex );

        /**
         * @brief Add a null vertex enrichment to the member data. returns the index of the vertex enrichment.
         *
         * @param aMeshIndex
         * @param aBaseInterpVertex
         * @param aNewVertex
         * @return moris_index index of the vertex enrichment
         */
        moris_index
        add_empty_vertex_enrichment( moris_index const & aMeshIndex,
                mtk::Vertex*                             aBaseInterpVertex,
                bool&                                    aNewVertex );

        //------------------------------------------------------------------------------

        void
        merge_duplicate_interpolation_vertices();

        //------------------------------------------------------------------------------

        void
        collect_base_vertex_to_enriched_vertex_connectivity( Vector< Vector< Interpolation_Vertex_Unzipped* > >& aBaseVertexToEnrichedVertex );

        //------------------------------------------------------------------------------

        /**
         * @brief Get the pointer to the vertex enrichment provided the vertex enrichment index.
         *
         * @param aMeshIndex
         * @param aVertexEnrichmentIndex
         * @return Vertex_Enrichment*
         */
        Vertex_Enrichment*
        get_vertex_enrichment( moris_index const & aMeshIndex,
                moris_index const &                aVertexEnrichmentIndex );

        //------------------------------------------------------------------------------

        /**
         * @brief Returns the vertex index corresponding to the vertex enrichment
         *
         * @param aMeshIndex
         * @param aVertexEnrichmentIndex
         * @return moris_index
         */
        moris_index
        get_vertex_related_to_vertex_enrichment( moris_index const & aMeshIndex,
                moris_index                                          aVertexEnrichmentIndex ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Converts a mesh index from an external mesh into the local mesh index
         *
         * @param aMeshIndex
         * @return moris_index
         */
        moris_index
        get_local_mesh_index_xtk( moris_index const & aMeshIndex ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Returns a cell of not owned vertex indices. These vertices do not have vertex interpolations and
         * need to be handled differently for ghost stabilization
         *
         * @return Cell< moris_index > const&
         */
        Vector< moris_index > const &
        get_not_owned_vertex_indices() const;

        //------------------------------------------------------------------------------

        void
        finalize_setup();

        void
        finalize_setup_new();

        //------------------------------------------------------------------------------

        bool
        verify_basis_support();

        //------------------------------------------------------------------------------

        // map setup
        void setup_local_to_global_maps();
        void setup_vertex_maps();
        void setup_vertex_to_bulk_phase();
        void setup_cell_maps();
        void setup_basis_maps();
        void setup_basis_ownership();
        void setup_basis_to_bulk_phase();
        void setup_mesh_index_map();

        //------------------------------------------------------------------------------

        // not owned vertex functions
        void setup_not_owned_vertices();


        //------------------------------------------------------------------------------
        // Parallel functions
        //------------------------------------------------------------------------------

        moris_id allocate_entity_ids(
                moris::size_t   aNumReqs,
                mtk::EntityRank aEntityRank,
                bool            aStartFresh );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        /**
         * @brief Assign IDs for all owned unzipped IP vertices and
         */
        void
        assign_ip_vertex_ids();

        //------------------------------------------------------------------------------

        /**
         * @brief Determine which UIPVs are owned by the executing processor and store them.
         * Also determine which UIPCs the not owned UIPVs are attached to.
         *
         * @param aUipcsAssociatedWithNotOwnedUipvs List of UIPCs the not owned UIPVs (as stored in mNotOwnedUnzippedVertices) are attached to.
         */
        void
        sort_unzipped_vertices_into_owned_and_not_owned(
                Vector< moris_index >& aUipcsAssociatedWithNotOwnedUipvs );

        //------------------------------------------------------------------------------

        /**
         * @brief Compile information by which not owned UIPVs can be identified by the owning processors
         *
         * @param aUipcsAssociatedWithNotOwnedUipvs List of UIPCs the not owned UIPVs (as stored in mNotOwnedUnzippedVertices) are attached to.
         * @param aNotOwnedUIPVsToProcs output: outer cell: index of the owning processor in comm-table || inner cell: index of the UIPV to be communicated
         * @param aBaseVertexIds output: outer cell: index of the owning processor in comm-table || inner cell: ID of the base vertex of the UIPV communicated
         * @param aUnzippedIpCellIds output: outer cell: index of the owning processor in comm-table || inner cell: ID of the IP cell the communicated UIPV is attached to
         */
        void
        prepare_requests_for_not_owned_unzipped_vertex_IDs(
                Vector< moris_index > const &    aUipcsAssociatedWithNotOwnedUipvs,
                Vector< Vector< moris_index > >& aNotOwnedUIPVsToProcs,
                Vector< Matrix< IdMat > >&       aBaseVertexIds,
                Vector< Matrix< IdMat > >&       aUnzippedIpCellIds );


        //------------------------------------------------------------------------------

        /**
         * @brief Find the IDs for UIPVs requested by other processors
         *
         * @param aVertIds output: outer cell: index of the processor received from in comm-table || inner cell: IDs of the UIPVs requested
         * @param aReceivedBaseVertexIds input: outer cell: index of the processor received from in comm-table || inner cell: ID of the base vertex of the UIPV communicated
         * @param aReceivedUnzippedIpCellIds input: outer cell: index of the processor received from in comm-table || inner cell: ID of the IP cell the communicated UIPV is attached to
         */
        void
        prepare_answers_for_owned_unzipped_vertex_IDs(
                Vector< Matrix< IdMat > >&        aVertIds,
                Vector< Matrix< IdMat > > const & aReceivedBaseVertexIds,
                Vector< Matrix< IdMat > > const & aReceivedUnzippedIpCellIds );

        //------------------------------------------------------------------------------

        /**
         * @brief Store the IDs received for not owned UIPVs
         *
         * @param aNotOwnedUIPVsToProcs input: outer cell: index of the owning processor in comm-table || inner cell: index of the UIPV communicated
         * @param aReceivedVertIds input: outer cell: index of the owning processor in comm-table || inner cell: IDs of the UIPVs communicated
         */
        void
        handle_requested_unzipped_vertex_ID_answers(
                Vector< Vector< moris_index > > const & aNotOwnedUIPVsToProcs,
                Vector< Matrix< IdMat > > const &       aReceivedVertIds );

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

      public:
        /**
         * @brief This function is for the debug purpose and to make sure that unenriched mesh has been enriched before
         *
         */
        void
        determine_unenriched_meshes_are_enriched_beforehand() const;

        //------------------------------------------------------------------------------

        /**
         * @brief Set the unenriched mesh indices object set function to populate unenriched mesh indices
         *
         * @param aMeshIndices
         */
        void
        set_unenriched_mesh_indices( Matrix< IndexMat > const & aMeshIndices );

        //------------------------------------------------------------------------------

        /**
         * @brief function that overwrites "some" maps to make the mesh unenriched, this function should be used with caution as it does not overwrite all
         * information
         *
         */
        void
        override_maps();

        //------------------------------------------------------------------------------

        /**
         * @brief this function replaces the t-matrices of the coefficients with their non-enriched version
         *
         */
        void
        override_vertex_enrichment_id_index();

        //------------------------------------------------------------------------------

        /**
         * @brief Get the set names object
         *
         * @param aSetEntityRank
         * @return Vector< std::string >
         */

        Vector< std::string >
        get_set_names( mtk::EntityRank aSetEntityRank ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Create a set names object
         *
         */

        void
        create_set_names();

        //------------------------------------------------------------------------------

        /**
         * @brief
         *
         * @param tFieldNames
         * @param tFieldData
         */

        void
        create_basis_support_fields( Matrix< DDRMat > const & aProbeSpheres );

        //------------------------------------------------------------------------------

        /**
         * @brief
         *
         * @param tFieldNames
         * @param tFieldData
         */

        void
        create_basis_function_fields( Matrix< DDRMat > const & aProbeSpheres );

        //------------------------------------------------------------------------------

        /**
         * @brief Create a cell id fields object
         *
         */

        void
        create_cell_id_fields();

        //------------------------------------------------------------------------------
        // Function to write a field
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

        /**
         * @brief Get the coefficient indices of node object
         *
         * @param aNodeIndex
         * @param aDiscretizationMeshIndex
         * @return Matrix< IndexMat >
         */

        Matrix< IndexMat >
        get_coefficient_indices_of_node(
                uint aNodeIndex,
                uint aDiscretizationMeshIndex );

        //------------------------------------------------------------------------------

        /**
         * @brief Get the t matrix of node loc ind object
         *
         * @param aNodeIndex
         * @param aDiscretizationMeshIndex
         * @return const Matrix< DDRMat >&
         */

        const Matrix< DDRMat >&
        get_t_matrix_of_node_loc_ind(
                uint aNodeIndex,
                uint aDiscretizationMeshIndex );

        //------------------------------------------------------------------------------

        /**
         * @brief Get the coefficient owners of node object
         *
         * @param aNodeIndex
         * @param aBSplineMeshIndex
         * @return Matrix< IdMat >
         */

        Matrix< IdMat >
        get_coefficient_owners_of_node(
                uint aNodeIndex,
                uint aBSplineMeshIndex );

        //------------------------------------------------------------------------------

        /**
         * @brief Get the coefficient IDs of node object
         *
         * @param aNodeIndex
         * @param aDiscretizationIndex
         * @return Matrix< IdMat >
         */

        Matrix< IdMat > get_coefficient_IDs_of_node(
                uint aNodeIndex,
                uint aDiscretizationIndex );

        //------------------------------------------------------------------------------

        /**
         * @brief Get the entity owner object
         *
         * @param aEntityIndex
         * @param aEntityRank
         * @param aDiscretizationMeshIndex
         * @return uint
         */

        uint get_entity_owner(
                moris_index       aEntityIndex,
                mtk::EntityRank   aEntityRank,
                const moris_index aDiscretizationMeshIndex ) const;


        //------------------------------------------------------------------------------
        // Additional Field Functions
        //------------------------------------------------------------------------------

        /**
         * @brief Get the field names object
         *
         * @param aEntityRank
         * @return Vector< std::string >
         */

        Vector< std::string >
        get_field_names( mtk::EntityRank aEntityRank );

        //------------------------------------------------------------------------------

        /**
         * @brief Create a field
         * @param[in] aLabel Field label
         * @param[in] aEntityRank Field entity rank
         * @param[in]aBulkPhaseIndex Bulk phase field defined over
         * aBulkphaseIndex of MORIS_INDEX_MAX results in a field over all phases
         */

        moris::moris_index
        create_field( std::string  aLabel,
                mtk::EntityRank    aEntityRank,
                moris::moris_index aBulkPhaseIndex = MORIS_INDEX_MAX );

        //------------------------------------------------------------------------------

        /**
         * @brief Given a field name and rank, gets the field index(ordinal)
         * @return Field index
         */

        moris::moris_index
        get_field_index( std::string aLabel,
                mtk::EntityRank      aEntityRank );

        //------------------------------------------------------------------------------

        /**
         * @brief Add field data to created field
         * @param[in] aFieldIndex Field index (use fn:get_field_index(...) to get)
         * @param[in] aEntityRank Field entity rank
         * @param[in] aFieldData Field data
         */

        void
        add_field_data( moris::moris_index aFieldIndex,
                mtk::EntityRank            aEntityRank,
                Matrix< DDRMat > const &   aFieldData );

        //------------------------------------------------------------------------------

        /**
         * @brief Given a field index and field entity rank, get the field data
         * @param[in] aFieldIndex Field index (use fn:get_field_index(...) to get)
         * @param[in] aEntityRank Field entity rank
         * @return Field data
         */

        Matrix< DDRMat > const &
        get_field_data( moris::moris_index aFieldIndex,
                mtk::EntityRank            aEntityRank ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the entity rank field index object
         *
         * @param aEntityRank
         * @return moris_index
         */

        moris_index
        get_entity_rank_field_index( mtk::EntityRank aEntityRank );

        //------------------------------------------------------------------------------

        /**
         * @return  whether a field exists or not
         */
        bool
        field_exists( std::string aLabel,
                mtk::EntityRank   aEntityRank );

        //------------------------------------------------------------------------------

        /**
         * @brief
         *
         * @param aParameterList
         */

        void
        write_mesh( moris::Parameter_List* aParamList );

        //------------------------------------------------------------------------------

        /**
         * @brief Get the num blocks object
         *
         * @return moris_index
         */

        moris::uint
        get_num_blocks() const;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the block set index object
         *
         * @param aString
         * @return moris_index
         */

        moris_index
        get_block_set_index( std::string aString ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Get the set cells object
         *
         * @param aSetLabel
         * @return Vector< mtk::Cell const * >
         */

        Vector< mtk::Cell const * >
        get_set_cells( std::string aSetLabel ) const;

        //------------------------------------------------------------------------------

        void
        update_communication_table( Vector< moris_id > const & aNewCommunicationTable );
    };
}    // namespace moris::xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_ */
