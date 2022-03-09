/*
 * cl_XTK_Enriched_Interpolation_Mesh.hpp
 *
 *  Created on: Jul 10, 2019
 *      Author: doble
 */

#ifndef PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_
#define PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_

#include "cl_Cell.hpp"

#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_XTK_Interpolation_Vertex_Unzipped.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Vertex_Enrichment.hpp"
#include "cl_XTK_Enrichment.hpp"

using namespace moris;

namespace xtk
{

    class Model;
    class Enriched_Interpolation_Mesh : public mtk::Interpolation_Mesh
    {
        public:
            Enriched_Interpolation_Mesh( Model* aXTKModel );
            ~Enriched_Interpolation_Mesh();
            //------------------------------------------------------------------------------
            // MTK Mesh Core Impl (see base class mtk::Mesh for function descriptions)
            //------------------------------------------------------------------------------
            MeshType                                    get_mesh_type() const;
            moris::uint                                 get_spatial_dim() const;
            uint                                        get_num_entities( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
            uint                                        get_max_num_coeffs_on_proc( const uint aBSplineMeshIndex ) const;
            Matrix< IndexMat >                          get_entity_connected_to_entity_loc_inds( moris_index aEntityIndex, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const;
            Matrix< IndexMat >                          get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const;
            moris_id                                    get_glb_entity_id_from_entity_loc_index( moris_index aEntityIndex, enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
            moris_index                                 get_loc_entity_ind_from_entity_glb_id( moris_id aEntityId, enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
            std::unordered_map< moris_id, moris_index > get_vertex_glb_id_to_loc_vertex_ind_map() const;
            Cell< mtk::Vertex const* >                  get_all_vertices() const;
            Matrix< IdMat >                             get_entity_connected_to_entity_glob_ids( moris_id aEntityId, enum EntityRank aInputEntityRank, enum EntityRank aOutputEntityRank, const moris_index aIndex = 0 ) const;
            Matrix< DDRMat >                            get_node_coordinate( moris_index aNodeIndex ) const;
            Matrix< DDRMat >                            get_base_node_coordinate( moris_index aBaseNodeIndex ) const;
            mtk::Vertex&                                get_mtk_vertex( moris_index aVertexIndex );
            mtk::Vertex const&                          get_mtk_vertex( moris_index aVertexIndex ) const;
            mtk::Cell const&                            get_mtk_cell( moris_index aElementIndex ) const;
            mtk::Cell&                                  get_mtk_cell( moris_index aElementIndex );
            mtk::Cell&                                  get_writable_mtk_cell( moris_index aElementIndex );
            Matrix< IdMat >                             get_communication_table() const;
            uint                                        get_num_elements();
            Matrix< IndexMat >                          get_element_indices_in_block_set( uint aSetIndex );
            moris_id                                    get_max_entity_id( enum EntityRank aEntityRank, const moris_index aIndex = 0 ) const;
            void                                        get_adof_map( const uint aBSplineIndex, map< moris_id, moris_index >& aAdofMap ) const;
            enum CellTopology                           get_blockset_topology( const std::string& aSetName );
            enum CellShape                              get_IG_blockset_shape( const std::string& aSetName );
            enum CellShape                              get_IP_blockset_shape( const std::string& aSetName );
            uint                                        get_node_owner( moris_index aNodeIndex ) const;
            uint                                        get_element_owner( moris_index aElementIndex ) const;

            //------------------------------------------------------------------------------
            // end mesh core functions
            //------------------------------------------------------------------------------

            //------------------------------------------------------------------------------
            // multigrid accessor functions
            //------------------------------------------------------------------------------

            uint                    get_num_interpolations();
            uint                    get_max_level( const moris_index aInterpolationIndex );
            uint                    get_num_basis( const moris_index aInterpolationIndex );
            uint                    get_basis_level( const moris_index aInterpolationIndex, const moris_index aBasisIndex );
            uint                    get_num_coarse_basis_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex );
            uint                    get_coarse_basis_index_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex, const moris_index aCoarseParentIndex );
            moris::Matrix< DDSMat > get_fine_basis_inds_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex );
            moris::Matrix< DDRMat > get_fine_basis_weights_of_basis( const moris_index aInterpolationIndex, const moris_index aBasisIndex );

            //------------------------------------------------------------------------------
            // Accessor functions of XTK specific data structures
            //------------------------------------------------------------------------------
            /*!
             * Get the enriched interpolation coefficients associated with a background coefficient
             */
            Matrix< IndexMat > const&
            get_enriched_coefficients_at_background_coefficient(
                    moris_index const& aMeshIndex,
                    moris_index                                                         aBackgroundCoeffIndex ) const;
            //------------------------------------------------------------------------------
            /*!
             * get the enriched coefficients at all the background mesh coefficients
             */
            Cell< Matrix< IndexMat > > const&
            get_enriched_coefficients_to_background_coefficients( moris_index const& aMeshIndex ) const;
            //------------------------------------------------------------------------------
            /*!
             * get the local enriched coefficient to global map
             */
            Matrix< IndexMat > const&
            get_enriched_coefficient_local_to_global_map( moris_index const& aMeshIndex ) const;
            //------------------------------------------------------------------------------
            /*!
             * Return the vector of background coefficient local to global
             */
            Matrix< IndexMat >
            get_background_coefficient_local_to_global_map() const;
            //------------------------------------------------------------------------------
            uint
            get_num_background_coefficients( moris_index const& aMeshIndex ) const;

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
            Cell< Interpolation_Cell_Unzipped* > const&
            get_enriched_interpolation_cells() const;

            //------------------------------------------------------------------------------
            /*!
             * Get the number of interpolation (t-matrices) defined on this mesh
             */
            uint
            get_num_interpolation_types() const;

            /*
             * Get the interpolation index for a local index
             */
            moris_index
            get_interpolation_index( moris_index const& aLocalInterpIndex ) const;

            /*!
             * get basis owner
             */
            moris_index
            get_basis_owner( 
                    moris_index aBasisIndex,
                    moris_index aMeshIndex );

            /*!
             * get basis bulk phase
             */
            moris_index
            get_basis_bulk_phase( 
                    moris_index const& aBasisIndex,
                    moris_index const& aMeshIndex ) const;

        public:
            Cell< Interpolation_Cell_Unzipped* >&
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
                    Cell< Interpolation_Cell_Unzipped* >&         aOwnedInterpCells,
                    Cell< Cell< Interpolation_Cell_Unzipped* > >& aNotOwnedInterpCells,
                    Cell< uint >&                                 aProcRanks );

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
             * @return moris::Cell< Interpolation_Cell_Unzipped const* > 
             */
            moris::Cell< Interpolation_Cell_Unzipped const* >
            get_enriched_cells_from_base_cells( moris::Cell< moris::mtk::Cell const* > const& aBaseCells ) const;

            //------------------------------------------------------------------------------

            /*!
             *  Single cell version of above
             */
            moris::Cell< Interpolation_Cell_Unzipped const* >
            get_enriched_cells_from_base_cell( moris::mtk::Cell const* aBaseCells ) const;

            //------------------------------------------------------------------------------

            Interpolation_Vertex_Unzipped const&
            get_xtk_interp_vertex( moris::uint aVertexIndex ) const;

            //------------------------------------------------------------------------------

            moris_index
            get_enr_basis_index_from_enr_basis_id(
                    moris_index const& aMeshIndex,
                    moris_index const& aBasisId ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief Convert a entity indices to entity ids
             * 
             * @param aIndices 
             * @param aEntityRank 
             * @return Matrix< IdMat > 
             */
            Matrix< IdMat > convert_indices_to_ids( 
                    Matrix< IndexMat > const& aIndices,
                    enum EntityRank           aEntityRank ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief Convert a entity ids to entity indices
             * 
             * @param aIds 
             * @param aEntityRank 
             * @return Matrix< IndexMat > 
             */
            Matrix< IndexMat > convert_ids_to_indices( 
                    Matrix< IdMat > const& aIds,
                    enum EntityRank        aEntityRank ) const;

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
             * @brief convert enriched //?
             * 
             * @param aMeshIndex 
             * @param aEnrichedIndices 
             * @param aEnrichedIds 
             */
            void
            convert_enriched_basis_indices_to_ids( moris_index const& aMeshIndex,
                    Matrix< IndexMat > const& aEnrichedIndices,
                    Matrix< IdMat >&          aEnrichedIds ) const;

            void
            convert_enriched_basis_indices_to_ids( moris_index const& aMeshIndex,
                    Cell< Matrix< IndexMat > > const& aEnrichedIndices,
                    Cell< Matrix< IdMat > >&          aEnrichedIds ) const;

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
            void print_enriched_verts_interpolation(const moris_index & aMeshIndex, std::string aFileName);
            void print_vertex_maps() const;
            void print_enriched_cell_maps() const;
            void print_basis_to_enriched_basis() const;
            void print_vertex_interpolation() const;
            void print_basis_information() const;

            // verification functions
            bool
            verify_basis_interpolating_into_cluster(
                    mtk::Cluster const&     aCluster,
                    moris_index const&      aMeshIndex,
                    const mtk::Master_Slave aIsMaster = mtk::Master_Slave::MASTER );

            // friend class
            friend class Enrichment;
            friend class Enriched_Integration_Mesh;
            friend class Ghost_Stabilization;

        protected:
            // Model pointer
            Model* mXTKModel;

            // basis rank
            enum EntityRank mBasisRank = EntityRank::INVALID;

            // mesh index
            Matrix< IndexMat >                             mMeshIndices;
            std::unordered_map< moris_index, moris_index > mMeshIndexToLocMeshIndex;   // over allocated

            // enriched interpolation vertices
            moris::uint                            mNumVerts;
            Cell< Interpolation_Vertex_Unzipped* > mEnrichedInterpVerts;    // over allocated

            // enriched interpolation cells
            moris::uint                          mNumVertsPerInterpCell;
            Cell< Interpolation_Cell_Unzipped* > mEnrichedInterpCells;     // over allocated

            // for each outer cell (base interpolation vertex), indices of enriched vertices
            Cell< Cell< Cell< moris_index > > > mBaseInterpVertToVertEnrichmentIndex; 
            // input: Discretization mesh index (DMI), BG mesh vertex index || output: list of enriched vertex indices

            // vertex enrichments or t-matrix
            // outer cell - mesh index (i.e. linear or quadratic b-spline enrichment)
            // inner cell
            Cell< Cell< Vertex_Enrichment* > > mInterpVertEnrichment; // input: DMI || output: list of T-Matrices ("Vertex Enrichments")

            // vertex enrichment to parent vertex index (these are enriched interpolation vertex indices)
            Cell< Cell< moris_index > > mVertexEnrichmentParentVertexIndex; // input: DMI, Vertex Enrichment index || output: index parent node on BG mesh

            // Bulk phase of each interpolation vertex
            Matrix< IndexMat > mVertexBulkPhase; // input: index of unzipped vertex || output: bulk-phase index
            Matrix< IndexMat > mVertexMaxSubphase; // input: index of unzipped vertex || output: maximum subphase ID on vertex // TODO: why do we need this?

            // basis coefficient to enriched basis coefficient
            Cell< moris::Cell< moris::Matrix< moris::IndexMat > > > mCoeffToEnrichCoeffs; 
            // input: DMI, index of non-enriched BF coefficient || ouput: list of enriched BF indices associated with it

            // local to global enriched basis vector
            Cell< moris::Matrix< moris::IdMat > >               mEnrichCoeffLocToGlob; // input: DMI, enriched BF index || output: global ID of that enriched BF 
            Cell< std::unordered_map< moris_id, moris_index > > mGlobaltoLocalBasisMaps; // input: DMI || output: map ordered by global BF IDs with corresponding local BF index

            // basis ownership
            Cell< moris::Matrix< moris::IdMat > > mEnrichCoeffOwnership; // input: DMI, enriched BF index || output:

            // basis bulk phase
            Cell< moris::Matrix< moris::IdMat > > mEnrichCoeffBulkPhase; // input: DMI, enriched BF index || output: Proc ID owning enr. BF

            // Entity maps
            Cell< Matrix< IdMat > >                             mLocalToGlobalMaps; // TODO input: DMI || output:
            Cell< std::unordered_map< moris_id, moris_index > > mGlobaltoLobalMaps; // TODO input: DMI || output:

            // base interpolation cells to their enriched interpolation cells
            moris::Cell< moris::Cell< Interpolation_Cell_Unzipped* > > mBaseCelltoEnrichedCell; // TODO input: DMI || output:

            // a connecitivty pointer that all the enriched interpolation cells use
            std::shared_ptr< moris::mtk::Cell_Info > mCellInfo;

            // Not owned vertex list
            Cell< moris_index > mNotOwnedVerts; // TODO input:  || output:

            // not owned basis functions
            Cell< moris_index > mNotOwnedBasis; // TODO input:  || output:
            Cell< moris_index > mOwnedBasis; // TODO input:  || output:

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
            add_vertex_enrichment( moris_index const& aMeshIndex,
                    mtk::Vertex*                          aBaseInterpVertex,
                    Vertex_Enrichment&                    aVertexEnrichment,
                    bool&                                 aNewVertex );

            //------------------------------------------------------------------------------

            void
            merge_duplicate_interpolation_vertices();

            //------------------------------------------------------------------------------

            void
            collect_base_vertex_to_enriched_vertex_connectivity( moris::Cell< moris::Cell< Interpolation_Vertex_Unzipped* > >& aBaseVertexToEnrichedVertex );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the pointer to the vertex enrichment provided the vertex enrichment index.
             * 
             * @param aMeshIndex 
             * @param aVertexEnrichmentIndex 
             * @return Vertex_Enrichment* 
             */
            Vertex_Enrichment*
            get_vertex_enrichment( moris_index const& aMeshIndex,
                    moris_index const&                    aVertexEnrichmentIndex );

            //------------------------------------------------------------------------------

            /**
             * @brief Returns the vertex index corresponding to the vertex enrichment
             * 
             * @param aMeshIndex 
             * @param aVertexEnrichmentIndex 
             * @return moris_index 
             */
            moris_index
            get_vertex_related_to_vertex_enrichment( moris_index const& aMeshIndex,
                    moris_index                                             aVertexEnrichmentIndex ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief Converts a mesh index from an external mesh into the local mesh index
             * 
             * @param aMeshIndex 
             * @return moris_index 
             */
            moris_index
            get_local_mesh_index( moris_index const& aMeshIndex ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief Returns a cell of not owned vertex indices. These vertices do not have vertex interpolations and
             * need to be handled differently for ghost stabilization
             * 
             * @return Cell< moris_index > const& 
             */
            Cell< moris_index > const&
            get_not_owned_vertex_indices() const;


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
            basis_exists_on_partition( moris_index const& aMeshIndex,
                    moris_index const&                        aBasisId );

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
                    moris_index const& aMeshIndex,
                    moris_index const& aBasisIdToAdd,
                    moris_index const& aBasisOwner,
                    moris_index const& aBasisBulkPhase );

            void
            finalize_setup();

            bool
            verify_basis_support();

            // map setup
            void setup_local_to_global_maps();
            void setup_vertex_maps();
            void setup_vertex_to_bulk_phase();
            void setup_cell_maps();
            void setup_basis_maps();
            void setup_basis_ownership();
            void setup_basis_to_bulk_phase();
            void setup_mesh_index_map();

            // not owned vertex functions
            void setup_not_owned_vertices();

            void
            assign_ip_vertex_ids();

            void
            sort_ip_vertices_by_owned_and_not_owned(
                    Cell< uint >&                             aOwnedVertices,
                    Cell< Cell< uint > >&                     aNotOwnedVertices,
                    Cell< Cell< uint > >&                     aNotOwnedIPCells,
                    Cell< uint >&                             aProcRanks,
                    std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData );

            void
            assign_owned_ip_vertex_ids(
                    Cell< uint > const& aOwnedIpVerts,
                    moris::moris_id&    aNodeId );

            void
            setup_outward_ip_vertex_requests(
                    Cell< Cell< uint > > const&               aNotOwnedIpVerts,
                    Cell< Cell< uint > > const&               aNotOwnedIpCells,
                    Cell< uint > const&                       aProcRanks,
                    std::unordered_map< moris_id, moris_id >& aProcRankToIndexInData,
                    Cell< Matrix< IndexMat > >&               aOutwardBaseVertexIds,
                    Cell< Matrix< IndexMat > >&               aOutwardIpCellIds );

            void
            prepare_ip_vertex_id_answers(
                    Cell< Matrix< IndexMat > >& aReceivedBaseVertexIds,
                    Cell< Matrix< IndexMat > >& aReceivedIpCellIds,
                    Cell< Matrix< IndexMat > >& aVertexIdAnswer );

            void
            handle_received_ip_vertex_ids(
                    Cell< Cell< uint > > const&       aNotOwnedVertices,
                    Cell< Matrix< IndexMat > > const& aReceivedVertexIds );

            /**
             * @brief This function communicates select vertex t-matrices
             *
             * @param aVerticesToCommunicate
             */
            void
            communicate_select_vertex_interpolation(
                    moris::Cell< mtk::Vertex* >  & aVerticesToCommunicate);

            void
            prepare_t_matrix_request_answers(
                    moris_index const&                aMeshIndex,
                    Cell< Matrix< IndexMat > > const& aRequestedEnrIPVertexIds,
                    Cell< Matrix< DDRMat > >&         aTMatrixWeights,
                    Cell< Matrix< IndexMat > >&       aTMatrixIndices,
                    Cell< Matrix< IndexMat > >&       aBasisOwners,
                    Cell< Matrix< IndexMat > >&       aTMatrixOffsets );

            void
            add_vertex_interpolation_to_communication_data(
                    moris::uint&        aCount,
                    Vertex_Enrichment*  aInterpolation,
                    Matrix< DDRMat >&   aTMatrixWeights,
                    Matrix< IndexMat >& aTMatrixIndices,
                    Matrix< IndexMat >& aTMatrixOwners,
                    Matrix< IndexMat >& aTMatrixOffsets );

            void
            handle_received_interpolation_data(
                    moris_index const&                aMeshIndex,
                    Cell< Matrix< IndexMat > > const& aVertexIdsToProc,
                    Cell< Matrix< DDRMat > > const&   aRequestedTMatrixWeights,
                    Cell< Matrix< IndexMat > > const& aRequestedTMatrixIndices,
                    Cell< Matrix< IndexMat > > const& aRequestedBasisOwners,
                    Cell< Matrix< IndexMat > > const& aRequestedTMatrixOffsets );

            void
            extract_vertex_interpolation_from_communication_data(
                    moris::uint const&          aNumVerts,
                    Matrix< DDRMat > const&     aTMatrixWeights,
                    Matrix< IndexMat > const&   aTMatrixIndices,
                    Matrix< IndexMat > const&   aTMatrixOwners,
                    Matrix< IndexMat > const&   aTMatrixOffsets,
                    Cell< Matrix< DDRMat > >&   aExtractedTMatrixWeights,
                    Cell< Matrix< IndexMat > >& aExtractedTMatrixIndices,
                    Cell< Matrix< IndexMat > >& aExtractedBasisOwners );

            //------------------------------------------------------------------------------
            // Parallel functions
            //------------------------------------------------------------------------------
            moris_id allocate_entity_ids(
                    moris::size_t aNumReqs,
                    enum EntityRank aEntityRank,
                    bool            aStartFresh );
    };
}// namespace xtk

#endif /* PROJECTS_XTK_SRC_XTK_CL_XTK_ENRICHED_INTERPOLATION_MESH_HPP_ */
