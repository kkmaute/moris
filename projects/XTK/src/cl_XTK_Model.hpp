/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jul 2, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_MODEL_HPP_
#define SRC_XTK_CL_XTK_MODEL_HPP_

// Standard Include
#include <limits>
#include <unordered_map>
#include <mpi.h>
#include <ctime>

// XTKL: Mesh Includes
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Sets_Info.hpp"
#include "cl_MTK_Fields_Info.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Tools includes
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"
#include "cl_Interpolaton.hpp"

#include "cl_XTK_Output_Options.hpp"

// Linalg Includes
#include "cl_Matrix.hpp"
#include "typedefs.hpp"
#include "fn_print.hpp"

#include "cl_Communication_Tools.hpp"

// Topology
//TODO: MOVE THESE WITH CUTTING METHODS SOMEWHERE ELSE
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Quad_4_Topology.hpp"
#include "cl_XTK_Hexahedron_8_Topology.hpp"

#include "fn_tet_volume.hpp"

// general geometry engine class
#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_Param_List.hpp"


namespace xtk
{
    class Enrichment;
    class Enrichment_Parameters;
    class Ghost_Stabilization;
    class Enriched_Interpolation_Mesh;
    class Enriched_Integration_Mesh;
    class Ghost_Stabilization;
    class Multigrid;
}


using namespace moris;

namespace xtk
{
    class Model
    {
        public:
            // Public member functions/data
            bool mVerbose = false;

            // friend classes
            friend class Enrichment;
            friend class Enriched_Interpolation_Mesh;
            friend class Enriched_Integration_Mesh;
            friend class Ghost_Stabilization;
            friend class Multigrid;

            //--------------------------------------------------------------------------------
            // Initialization
            //--------------------------------------------------------------------------------
            Model(){};

            //--------------------------------------------------------------------------------
            /*!
             * @brief Constructor for standalone use of XTK
             *
             */
            Model(uint aModelDimension,
                    moris::mtk::Interpolation_Mesh* aMeshData,
                    moris::ge::Geometry_Engine* aGeometryEngine,
                    bool aLinkGeometryOnConstruction = true);

            //--------------------------------------------------------------------------------

            /*!
             * @brief Parameter list based XTK initialization
             */
            Model(moris::ParameterList const & aParameterList);

            // Indicates the background mesh and the geometry are the same thing
            //FIXME:still necessary?
            bool mSameMesh;

            //--------------------------------------------------------------------------------

            ~Model();

            //--------------------------------------------------------------------------------
            /*!
             * @brief Set the pointer to the geometry engine
             */
            void
            set_geometry_engine(moris::ge::Geometry_Engine* aGeometryEngine);

            //--------------------------------------------------------------------------------
            /*!
             * @brief Set the pointer to the background interpolation mesh
             */
            void
            set_mtk_background_mesh(moris::mtk::Interpolation_Mesh* aMesh);

            //--------------------------------------------------------------------------------

            /*!
             * @brief Set the pointer to the output performer
             */
            void
            set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

            //--------------------------------------------------------------------------------

            /*!
             * @brief Set the pointer to the input performer
             */
            void
            set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

            //--------------------------------------------------------------------------------

            /*!
             * @brief Initialize data using the interpolation mesh
             */
            void
            initialize( moris::mtk::Interpolation_Mesh* aMesh );


            //--------------------------------------------------------------------------------
            // Operations
            //--------------------------------------------------------------------------------
            /*!
             * @brief Workflow based perform call. Uses parameter list to decide operations to
             * perform
             */
            bool
            perform();

            //--------------------------------------------------------------------------------
            /*!
             * Decomposes a mesh to conform to a geometric interface
             * @param[in] aMethods Subdivision Methods
             */
            bool decompose(Cell<enum Subdivision_Method> aMethods);

            //--------------------------------------------------------------------------------
            /*!
            * Remove child meshes that have an intersection but all cells are in the same bulk phase
            */
            void
            set_cleanup_cut_mesh_flag( bool aCleanupMesh );

            //--------------------------------------------------------------------------------
            /*!
             * @return Timing data
             */
            Matrix<DDRMat>
            get_timing_data() const;

            //--------------------------------------------------------------------------------

            /*!
             * @return Timing data labels
             */
            Cell<std::string>
            get_timing_labels() const;

            //--------------------------------------------------------------------------------

            /*!
             * @brief Save timing data to hdf5 file
             */
            void
            save_timing_to_hdf5( const std::string & aFilePath ) const;

            //--------------------------------------------------------------------------------

            /*!
             * @brief Save model statistics to file
             */
            void
            save_model_statistics_to_file( const std::string & aFilePath );

            //--------------------------------------------------------------------------------

            /*!
             * @brief Cleanly print XTK timing data to the log
             */
            void
            print_timing_data() const;

            //--------------------------------------------------------------------------------
            /*!
            * @brief get the memory usage of XTK
            */
            moris::Memory_Map
            get_memory_usage();

            //--------------------------------------------------------------------------------
            //FIXME: REMOVE and related functions in child mesh
            /*!
             * Uses sub-phase information within a child mesh to construct one interpolation element for each sub-phase cluster
             */
            void
            unzip_child_mesh();

            // ----------------------------------------------------------------------------------

            // FIXME: REMOVE
            /*!
             * Unzipps the interface (must be called after decompose but prior
             * to enrichment)
             */
            void
            unzip_interface();

            // ----------------------------------------------------------------------------------
 
            /*!
             * Perform generalized heaviside  enrichment.
             */
            void
            perform_basis_enrichment(
                    enum EntityRank  const & aBasisRank,
                    moris_index      const & aMeshIndex = 0);

            // ----------------------------------------------------------------------------------
 
            /*!
             * @brief Perform basis enrichment
             * @param[in] aBasisRank Entity rank of basis functions
             * @param[in] aMeshIndex Mesh index
             */
            void
            perform_basis_enrichment(
                    enum EntityRank  const & aBasisRank,
                    Matrix<IndexMat> const & aMeshIndex);

            // ----------------------------------------------------------------------------------

            void perform_hanging_node_identification();

            // ----------------------------------------------------------------------------------
            /*!
             * @brief Probes and prints the information about a background cell
             */
            void
            probe_bg_cell(Matrix<IndexMat> const & tBGCellIds);

            // ----------------------------------------------------------------------------------

            /*!
             * @return Basis enrichment
             */
            Enrichment const &
            get_basis_enrichment();

            // ----------------------------------------------------------------------------------
            /*!
             * @param[in] aIndex Interpolation Index
             * @return Enriched interpolation mesh
             */
            Enriched_Interpolation_Mesh &
            get_enriched_interp_mesh(moris::moris_index aIndex = 0);

            // ----------------------------------------------------------------------------------
            /*!
             * @param[in] aIndex Interpolation Index
             * @return Enriched integration mesh
             */
            Enriched_Integration_Mesh &
            get_enriched_integ_mesh(moris::moris_index  aIndex = 0);


            // ----------------------------------------------------------------------------------
            /*!
             * @brief Constructs the face oriented ghost penalization
             */
            void
            construct_face_oriented_ghost_penalization_cells();

            // ----------------------------------------------------------------------------------
            /*!
             * @return Ghost stabilization
             */
            Ghost_Stabilization &
            get_ghost_stabilization(moris::moris_index  aIndex = 0);

            // ----------------------------------------------------------------------------------
            //FIXME: REMOVE
            /*!
             * Convert Tet4 elements to Tet10 Elements
             */
            void
            convert_mesh_tet4_to_tet10();

            // ----------------------------------------------------------------------------------

            /*!
             * @brief Construct multigrid information
             */
            void construct_multigrid();

            //--------------------------------------------------------------------------------
            // Member data access functions
            //--------------------------------------------------------------------------------
            /*!
             * @return Cut mesh (collection of child meshes)
             */
            Cut_Mesh &
            get_cut_mesh();

            // ----------------------------------------------------------------------------------
            /*!
             * @return Const cut mesh (collection of child meshes)
             */
            Cut_Mesh const &
            get_cut_mesh() const;

            // ----------------------------------------------------------------------------------
            /*!
             * @return Background mesh
             */
            Background_Mesh &
            get_background_mesh();

            // ----------------------------------------------------------------------------------
            /*!
             * @return Constant background mesh
             */
            Background_Mesh const &
            get_background_mesh() const;

            // ----------------------------------------------------------------------------------
            /*!
             * @return Geometry engine pointer
             */
            moris::ge::Geometry_Engine*
            get_geom_engine();

            // ----------------------------------------------------------------------------------
            // Outputting functions
            //-----------------------------------------------------------------------------------
            /*!
             * @brief Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
             * @param[in] aOutputOptions Mesh output options
             */
            moris::mtk::Integration_Mesh*
            get_output_mesh(Output_Options const & aOutputOptions = Output_Options());

            //--------------------------------------------------------------------------------
            // Cell Neighborhood creation and access
            //--------------------------------------------------------------------------------
            /*!
             * @brief constructs the cell neighborhood
             */
            void
            construct_neighborhood();

            //-----------------------------------------------------------------------------------
            /*!
             * @brief constructs subphase neighborhood
             */
            void
            construct_subphase_neighborhood();

            //-----------------------------------------------------------------------------------
            /*!
             * @brief constructs the simple cell neighborhood between intersected background cells
             */
            void
            construct_cut_mesh_simple_neighborhood();

            //-----------------------------------------------------------------------------------

            /*!
             * @brief construct neighborhood from cut to uncut cells
             */
            void
            construct_cut_mesh_to_uncut_mesh_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace);

            //-----------------------------------------------------------------------------------
            /*!
             * @briefHanging node considerations for cell neighborhood
             */
            void
            construct_complex_neighborhood();

            //-----------------------------------------------------------------------------------
            /*!
             * @briefConstruct neighborhood between uncut cells
             */
            void
            construct_uncut_neighborhood(moris::Cell<moris::Cell<moris_index>> & aCutToUncutFace);

            //-----------------------------------------------------------------------------------
            /*!
             *@brief Delete cell neighborhood data
             */
            void
            delete_neighborhood();

            //--------------------------------------------------------------------------------
            // Data access functions
            //--------------------------------------------------------------------------------

            /*!
             * get spatial dimension of model
             */
            moris::uint
            get_spatial_dim() const;

            //-----------------------------------------------------------------------------------

            /*!
             * @return the number of elements in the entire model
             * includes all child elements and all background elements (combination)
             */
            moris::uint
            get_num_elements_total();

            //-----------------------------------------------------------------------------------

            /*!
             * @return the number of elements in unzipped mesh (total - num child meshes)
             */
            moris::uint
            get_num_elements_unzipped();

            /*!
             * @returns element to element connectivity
             */
            moris::Cell<moris::Cell<moris::mtk::Cell*>> const &
            get_element_to_element(){ return mElementToElement; };

             /*!
             * @returns element to element my side ordinals
             */
            moris::Cell<moris::Cell<moris_index>> const &
            get_element_to_element_my_side_ords(){ return mElementToElementSideOrds; };

            /*!
             * @returns element to element neighbor side ordinals
             */
            moris::Cell<moris::Cell<moris_index>> const &
            get_element_to_element_neighbor_side_ords(){ return mElementToElementNeighborSideOrds; };

            //-----------------------------------------------------------------------------------

            moris_index
            get_cell_xtk_index(moris_id aCellId);

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase();

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase_my_side_ords();

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase_neighbor_side_ords();

            //-----------------------------------------------------------------------------------

            moris::Cell<moris::Cell<moris_index>>  const &
            get_subphase_to_subphase_transition_loc();

            //-----------------------------------------------------------------------------------

            bool
            subphase_is_in_child_mesh(moris_index aSubphaseIndex);

            //-----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat>
            get_element_to_subphase();

            //-----------------------------------------------------------------------------------

            moris_id
            get_subphase_id(moris_id aSubphaseIndex);

            //-----------------------------------------------------------------------------------

            moris_index
            get_subphase_index(moris_id aSubphaseId);

            //-----------------------------------------------------------------------------------

            std::shared_ptr< Multigrid > get_multigrid_ptr();

            //--------------------------------------------------------------------------------
            // Printing Functions
            //--------------------------------------------------------------------------------
            void
            print_cells();

            //------------------------------------------------------------------------------

            void
            print_neighborhood();

            //------------------------------------------------------------------------------

            void
            print_subphase_neighborhood();

            //------------------------------------------------------------------------------

            void
            print_interface_vertices();

            //------------------------------------------------------------------------------

        protected:
            moris::ParameterList               mParameterList;

            uint                               mModelDimension;
            Background_Mesh                    mBackgroundMesh;
            Cut_Mesh                           mCutMesh;

            moris::ge::Geometry_Engine*        mGeometryEngine;

            Enrichment*                        mEnrichment;
            Ghost_Stabilization*               mGhostStabilization;
            Cell<Enriched_Interpolation_Mesh*> mEnrichedInterpMesh;
            Cell<Enriched_Integration_Mesh*>   mEnrichedIntegMesh;

            std::shared_ptr< xtk::Multigrid >  mMultigrid;

        private:

            // XTK Model State Flags
            bool mDecomposed        = false; // Model has been decomposed
            bool mConvertedToTet10s = false; // Model has been converted from tet4's to tet10s
            bool mMeshDataFinalized = false;
            bool mEnriched          = false; // Model has been enriched
            bool mUnzipped          = false; // Model has been unzipped
            bool mGhost             = false; // Model has setup ghost stabilization

            // Flag to cleanup mesh at end of decomposition
            bool mTriangulateAll = false; // Triangulate all background cells
            bool mCleanupMesh    = false; // Cleanup the mesh 

            // cell map
            std::map< moris_id, moris_index> mCellGlbToLocalMap;

            // The midside nodes are stored here currently but this may change
            moris::Matrix< moris::IndexMat > mMidsideElementToNode;

            // element to element neighborhood
            moris::Cell<moris::Cell<moris::mtk::Cell*>> mElementToElement;
            moris::Cell<moris::Cell<moris_index>>       mElementToElementSideOrds;
            moris::Cell<moris::Cell<moris_index>>       mElementToElementNeighborSideOrds;
            moris::Cell<moris::Cell<moris_index>>       mSubphaseToSubPhase;
            moris::Cell<moris::Cell<moris_index>>       mSubphaseToSubPhaseMySideOrds;
            moris::Cell<moris::Cell<moris_index>>       mSubphaseToSubPhaseNeighborSideOrds;

            // in the case of a hierarchically refined mesh, there are transitions with hanging nodes
            // this data flags the transition from a large facet to a smaller facet. (this is trivial
            // for non hmr type meshes)
            moris::Cell<moris::Cell<moris_index>> mTransitionNeighborCellLocation;

            // local to global subphase map
            std::unordered_map<moris::moris_id,moris::moris_index> mGlobalToLocalSubphaseMap;

            std::shared_ptr< mtk::Mesh_Manager > mMTKInputPerformer  = nullptr;
            std::shared_ptr< mtk::Mesh_Manager > mMTKOutputPerformer = nullptr;

            bool mInitializeCalled = false;

            // timing information
            // the time
            Cell<real>        mTimingData;

            // the label of the time
            Cell<std::string> mTimingLabels;

            // category that the timing is in (i.e. decomp, enrich, ghost, overall)
            Cell<std::string> mTimingCategory;

            // Private Functions
        private:
            //------------------------------------------------------------------------------

            // Functions that take parameter inputs and setup XTK inputs
            bool
            has_parameter_list();


            //------------------------------------------------------------------------------

            /*!
             * Verifys provided parameter list
             */
            bool
            valid_parameters();

            //------------------------------------------------------------------------------

            /*
             * Using the parameter list, figure out the cell of subdivision methods
             */
            Cell<enum Subdivision_Method>
            get_subdivision_methods();

            //------------------------------------------------------------------------------
            // Internal Decomposition Functions
            //------------------------------------------------------------------------------

            /*!
             * formulates node requests in the geometry objects. Dependent on the type of decomposition
             * @param[in] aReqType- specifies which template mesh is going to be used later on
             */
            void decompose_internal(
                    enum Subdivision_Method    const & aSubdivisionMethod,
                    moris::uint                        aGeomIndex,
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    bool                       const & aFirstSubdivision   = true,
                    bool                       const & aSetIds = false);

            //------------------------------------------------------------------------------

            /*!
             * Regular subdivision for a 3D mesh
             */
            void
            decompose_internal_reg_sub_hex8(
                    moris::uint                        aGeomIndex,
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    bool                       const & aFirstSubdivision,
                    bool                       const & aSetIds );

            //------------------------------------------------------------------------------

            void
            decompose_internal_reg_sub_hex8_make_requests(
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    moris::Matrix< moris::IndexMat > & tNewPairBool,
                    Decomposition_Data               & tDecompData);

            //------------------------------------------------------------------------------

            /*!
             * Regular subdivision for a 2D mesh
             */
            void
            decompose_internal_reg_sub_quad4(
                    moris::uint                        aGeomIndex,
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    bool                       const & aFirstSubdivision,
                    bool                       const & aSetIds );

            //------------------------------------------------------------------------------

            void
            decompose_internal_reg_sub_quad4_make_requests(
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    moris::Matrix< moris::IndexMat > & tNewPairBool,
                    Decomposition_Data               & tDecompData);

            //------------------------------------------------------------------------------

            void
            decompose_internal_set_new_nodes_in_child_mesh_reg_sub(
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    moris::Matrix< moris::IndexMat > & tNewPairBool,
                    moris::real                        tNumParamCoords,
                    Decomposition_Data &               tDecompData);

            //------------------------------------------------------------------------------

            void
            decompose_internal_set_new_nodes_in_child_mesh_nh(
                    moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                    Decomposition_Data &               tDecompData);

            //------------------------------------------------------------------------------

            void
            create_new_node_association_with_geometry(Decomposition_Data & tDecompData);

            //------------------------------------------------------------------------------

            void
            catch_all_unhandled_interfaces( );


            //------------------------------------------------------------------------------
            bool
            check_for_degenerated_cells();

            //------------------------------------------------------------------------------

            bool
            check_for_all_cell_vertices_on_interface();

            //------------------------------------------------------------------------------

            /*
             * Parallel assignment of node request identifiers
             */
            void
            assign_node_requests_identifiers(
                    Decomposition_Data & aDecompData,
                    moris::moris_index   aMPITag);

            //------------------------------------------------------------------------------

            void
            sort_new_node_requests_by_owned_and_not_owned(
                    Decomposition_Data                    & tDecompData,
                    Cell<uint>                            & aOwnedRequests,
                    Cell<Cell<uint>>                      & aNotOwnedRequests,
                    Cell<uint>                            & aProcRanks,
                    std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData);

            //------------------------------------------------------------------------------

            void
            assign_owned_request_id(
                    Decomposition_Data & aDecompData,
                    Cell<uint> const   & aOwnedRequest,
                    moris::moris_id    & aNodeId);

            //------------------------------------------------------------------------------

            void assign_index(Decomposition_Data & aDecompData);

            //------------------------------------------------------------------------------

            /*!
             * Sets up the decomposition data request for node identifiers from owning proc
             *
             */
            void
            setup_outward_requests(Decomposition_Data              const & aDecompData,
                    Cell<Cell<uint>>                const & aNotOwnedRequests,
                    Cell<uint>                      const & aProcRanks,
                    std::unordered_map<moris_id,moris_id> & aProcRankToIndexInData,
                    Cell<Matrix<IndexMat>>                & aSentRequests);

             //------------------------------------------------------------------------------

             /*!
              * Verifies that the nodes in a decomposition have all been assigned node ids
              */
             bool
             verify_successful_node_assignment(Decomposition_Data & aDecompData);

             //------------------------------------------------------------------------------

        protected:
             void
             send_outward_requests(
                     moris_index            const & aMPITag,
                     Cell<uint>             const & aProcRanks,
                     Cell<Matrix<IndexMat>> & aOutwardRequests);

             //------------------------------------------------------------------------------

             void
             inward_receive_requests(
                     moris_index            const & aMPITag,
                     moris::uint                    aNumRows,
                     Cell<Matrix<IndexMat>> &       aReceivedData,
                     Cell<uint>             &       aProcRanksReceivedFrom);

             //------------------------------------------------------------------------------

             void
             prepare_request_answers(
                     Decomposition_Data           & aDecompData,
                     Cell<Matrix<IndexMat>> const & aReceiveData,
                     Cell<Matrix<IndexMat>>       & aReceivedRequestAnswers);

             //------------------------------------------------------------------------------

             void
             return_request_answers(
                     moris_index            const & aMPITag,
                     Cell<Matrix<IndexMat>> const & aRequestAnswers,
                     Cell<uint>             const & aProcRanks);

             //------------------------------------------------------------------------------

             void
             inward_receive_request_answers(
                     moris_index            const & aMPITag,
                     moris::uint            const & aNumRows,
                     Cell<uint>             const & aProcRanks,
                     Cell<Matrix<IndexMat>> &       aReceivedData);

             //------------------------------------------------------------------------------

             void
             handle_received_request_answers(
                     Decomposition_Data           & aDecompData,
                     Cell<Matrix<IndexMat>> const & aRequests,
                     Cell<Matrix<IndexMat>> const & aRequestAnswers,
                     moris::moris_id              & aNodeId);

             //------------------------------------------------------------------------------
             // moris real versions of above
             void
             send_outward_requests_reals(
                     moris_index const    & aMPITag,
                     Cell<uint>  const    & aProcRanks,
                     Cell<Matrix<DDRMat>> & aOutwardRequests);

             //------------------------------------------------------------------------------

             void
             inward_receive_requests_reals(
                     moris_index const &    aMPITag,
                     moris::uint            aNumRows,
                     Cell<Matrix<DDRMat>> & aReceivedData,
                     Cell<uint>           & aProcRanksReceivedFrom);

             //------------------------------------------------------------------------------

             void
             return_request_answers_reals(
                     moris_index const & aMPITag,
                     Cell<Matrix<DDRMat>> const & aRequestAnswers,
                     Cell<uint>              const & aProcRanks);

             //------------------------------------------------------------------------------

             void
             inward_receive_request_answers_reals(
                     moris_index          const & aMPITag,
                     moris::uint          const & aNumRows,
                     Cell<uint>           const & aProcRanks,
                     Cell<Matrix<DDRMat>>       & aReceivedData);

             //------------------------------------------------------------------------------

             void
             add_timing_data(
                     real        const & aTime,
                     std::string const & aLabel,
                     std::string const & aCategory);

        private:

             /*
              * Perform all tasks needed to finalize the decomposition process. Assignes cells
              */
             void
             finalize_decomp();

             void
             finalize_mesh_data();

             //------------------------------------------------------------------------------

            /*!
            * assign child element indices
            */
            void
            assign_child_element_indices(bool aUpdateAvailable );

             /*!
              * assign child element global ids
              */
             void
             assign_child_element_ids();

             //------------------------------------------------------------------------------

             void
             prepare_child_element_identifier_requests(
                     Cell<Cell<moris_id>>                  & aNotOwnedChildMeshesToProcs,
                     Cell<moris::Matrix<IdMat>>            & aOwnedParentCellId,
                     Cell<moris::Matrix<IdMat>>            & aNumOwnedCellIdsOffsets,
                     Cell<uint >                           & aProcRanks,
                     std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex);

             //------------------------------------------------------------------------------

             void
             prepare_child_cell_id_answers(Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
                     Cell<Matrix<IndexMat>> & aReceivedParentCellNumChildren,
                     Cell<Matrix<IndexMat>> & aChildCellIdOffset);

             //------------------------------------------------------------------------------

             void
             handle_received_child_cell_id_request_answers(
                     Cell<Cell<moris_index>> const & aChildMeshesInInNotOwned,
                     Cell<Matrix<IndexMat>>  const & aReceivedChildCellIdOffset);

             //------------------------------------------------------------------------------

             /*!
              * Add children elements to local to global map
              */
              void
              add_child_elements_to_local_to_global_map();

              //------------------------------------------------------------------------------

              /*!
               * Constructs the mtk cell interface for all child elements created during the
               * decomposition process
               */
              void
              create_child_element_mtk_cells();

              //------------------------------------------------------------------------------

              /*!
               * Add the vertex pointer to child meshes
               */
              void
              add_vertices_to_child_meshes();

              //------------------------------------------------------------------------------

              /*!
               * setup cell id to index map
               */
              void
              setup_cell_glb_to_local_map();

              //------------------------------------------------------------------------------

              /*
               *
               * Identifies local sub-phase clusters within a single child mesh.
               * The sub-phase data (result of floodfill) is stored as a member
               * variable in each child mesh as sub-phase bins
               */
              void
              identify_local_subphase_clusters_in_child_meshes();

              //------------------------------------------------------------------------------
              /*
               * Creates double side set data for interfaces internal to the child mesh
               */ 
              void
              construct_internal_double_sides_between_subphases();

              //------------------------------------------------------------------------------

              void
              prepare_subphase_identifier_requests(
                      Cell<Cell<moris_id>>       & aNotOwnedSubphasesToProcs,
                      Cell<Cell<moris_id>>       & aSubphaseCMIndices,
                      Cell<moris::Matrix<IdMat>> & aParentCellIds,
                      Cell<moris::Matrix<IdMat>> & aChildCellIds,
                      Cell<moris::Matrix<IdMat>> & aNumChildCellsInSubphase,
                      Cell<uint>                 & aProcRanks,
                      std::unordered_map<moris_id,moris_id> & aProcRankToDataIndex);

              //------------------------------------------------------------------------------

              void
              prepare_subphase_id_answers(
                      Cell<Matrix<IndexMat>> & aReceivedParentCellIds,
                      Cell<Matrix<IndexMat>> & aFirstChildCellIds,
                      Cell<Matrix<IndexMat>> & aReceivedNumChildCellsInSubphase,
                      Cell<Matrix<IndexMat>> & aSubphaseIds);

              void
              handle_received_subphase_id_request_answers( 
                      Cell<Cell<moris_index>>    const & aChildMeshesInNotOwned,
                      Cell<Cell<moris_index>>    const & aCMSubphaseIndices,
                      Cell<Matrix<IndexMat>>     const & aReceivedSubphaseIds);

              //------------------------------------------------------------------------------

              void
              assign_subphase_glob_ids();

              //------------------------------------------------------------------------------

              void
              setup_glob_to_loc_subphase_map();

              //------------------------------------------------------------------------------

              void
              unzip_child_mesh_internal();

              //------------------------------------------------------------------------------

              /**
               * Take the interface faces and create collapsed prisms
               */
              void
              unzip_interface_internal();

              //------------------------------------------------------------------------------

              void
              unzip_interface_internal_assign_node_identifiers(
                      moris::uint aNumNodes,
                      moris::Matrix<moris::IdMat> & aUnzippedNodeIndices,
                      moris::Matrix<moris::IdMat> & aUnzippedNodeIds);

              //------------------------------------------------------------------------------

              void
              unzip_interface_internal_modify_child_mesh(
                      moris::uint                         aGeometryIndex,
                      moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                      moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                      moris::Matrix<moris::IdMat> const & aUnzippedNodeIds);

              //------------------------------------------------------------------------------

              moris::Matrix< moris::IndexMat >
              unzip_interface_internal_assign_which_element_uses_unzipped_nodes(
                      moris::moris_index                       aGeometryIndex,
                      moris::Matrix< moris::IndexMat > const & aInterfaceElementPairs );

              //------------------------------------------------------------------------------

              moris::Cell<moris::Cell< moris::moris_index >>
              unzip_interface_internal_collect_child_mesh_to_interface_node(
                      moris::Matrix<moris::IdMat> const & aInterfaceNodeIndices,
                      moris::Matrix<moris::IdMat> const & aUnzippedNodeIndices,
                      moris::Matrix<moris::IdMat> const & aUnzippedNodeIds);

              void
              unzip_interface_construct_interface_elements(
                      moris::uint aGeometryIndex,
                      moris::Matrix< moris::IndexMat > const & aElementPairs,
                      moris::Matrix< moris::IndexMat > const & aSideOrdinalPairs);

              //------------------------------------------------------------------------------

              void
              unzip_interface_assign_element_identifiers();
              //------------------------------------------------------------------------------
              // Enrichment computation functions
              //------------------------------------------------------------------------------

              void
              perform_basis_enrichment_internal(
                      enum EntityRank   const & aBasisRank,
                      Matrix<IndexMat> const & aMeshIndex);

              //------------------------------------------------------------------------------

              /*!
               * Links the vertex enrichment to the mtk implementation of the vertex
               */
              void
              link_vertex_enrichment_to_vertex_interpolation();

              //------------------------------------------------------------------------------

              //------------------------------------------------------------------------------
              // internal ghost functions
              //------------------------------------------------------------------------------

              /*
               * Constructs child mesh groups (owned, owned shared, not owned shared)
               */
              void
              sort_children_meshes_into_groups();

              //------------------------------------------------------------------------------

              /*
               * For nodes that are created during the decomposition process, tell
               * the XTK mesh about where they live in child meshes.
               */
              void
              associate_nodes_created_during_decomp_to_child_meshes();

              //------------------------------------------------------------------------------

              /*
               * Set element phase index
               */
              void
              set_element_phases();

              //------------------------------------------------------------------------------

              /*
               * Tells the XTK mesh about where it's children live in the cut mesh
               */
              void
              set_downward_inheritance();

              //------------------------------------------------------------------------------
              /*
               * This algorithm sets up the active child mesh indices and registers new pairs in the downward inheritance
               */

              void
              run_first_cut_routine(
                      moris::uint                        aGeomIndex,
                      moris::Matrix< moris::IndexMat > & aActiveChildMeshIndices,
                      moris::Matrix< moris::IndexMat > & aNewPairBool);

             //------------------------------------------------------------------------------
              /*
               * Returns a flag that all intersected cells are on the same level if false the decomposition will return false
               */
              bool
              all_child_meshes_on_same_level();

              //------------------------------------------------------------------------------

              /*!
               * Constructs the output mesh using provided Output_Options
               */
              moris::mtk::Integration_Mesh*
              construct_output_mesh( Output_Options const & aOutputOptions );

              //------------------------------------------------------------------------------

              /*!
               * Setup interface single sided side sets
               */
              void
              setup_interface_single_side_sets(
                      Output_Options                   const  & aOutputOptions,
                      Cell<moris::Matrix<moris::IdMat>>       & aCellIdsAndSideOrds,
                      Cell<std::string>                       & aInterfaceSetNames);

              //------------------------------------------------------------------------------

              /*!
               * Takes the whole node local to global map and removes the nodes
               * which are not part of the phase being output
               */
              moris::Matrix<moris::IndexMat>
              get_node_map_restricted_to_output_phases(
                      Output_Options const &           aOutputOptions,
                      moris::Matrix<moris::IndexMat> & aOutputtedNodeInds);

              //------------------------------------------------------------------------------

              moris::Cell<std::string>
              assign_geometry_data_names();

              //------------------------------------------------------------------------------


              moris::Cell < enum moris::EntityRank >
              assign_geometry_data_field_ranks();

              //------------------------------------------------------------------------------

              /*!
               * Sets up background node sets for mesh output. Propogates the node set from
               * the background mesh to the output mesh
               */
              moris::Cell<moris::mtk::MtkNodeSetInfo>
              propogate_background_node_sets(
                      moris::Cell<moris::Matrix<IndexMat>>       & aNodeSetData,
                      Output_Options                       const & aOutputOptions);

              //------------------------------------------------------------------------------

              /*!
               * Sets up background side sets for mesh output. Propogates the side set from
               * the background mesh to the output mesh
               */
              moris::Cell<moris::mtk::MtkSideSetInfo>
              propogate_background_side_sets(
                      moris::Cell<moris::Matrix<IndexMat>>       & aSideSetData,
                      Output_Options                       const & aOutputOptions);

              //------------------------------------------------------------------------------

              /*!
               * Add a single side set from the background mesh
               */
              void
              propogate_background_side_set(
                      std::string             const &             aSideSetName,
                      moris::moris_index                          aNoChildIndex,
                      moris::moris_index                          aChildIndex,
                      moris::Cell<moris::Matrix<IndexMat>>      & aElementIdsAndSideOrd,
                      moris::Cell<moris::mtk::MtkSideSetInfo>   & aSideSetData,
                      Output_Options          const             & aOutputOptions,
                      bool                                        aOutputIndices);

              //------------------------------------------------------------------------------

              /*!
               * This function checks for the side sets which appear in a mesh that comes
               * from a SEACAS generated string. For example:
               * "generated:1x1x1|sideset:xXyYzZ"
               */
              moris::Cell<std::string>
              check_for_and_remove_internal_seacas_side_sets(moris::Cell<std::string> & aSideSetNames);

              //------------------------------------------------------------------------------

              /*!
               * Combine interface and non-interface blocks
               */
              Cell<moris::Matrix<moris::IdMat>>
              combine_interface_and_non_interface_blocks(
                      Cell<moris::Matrix<moris::IdMat>> & tChildElementsByPhase,
                      Cell<moris::Matrix<moris::IdMat>> & tNoChildElementsByPhase);

              //------------------------------------------------------------------------------

              /*
               * Returns the number of phases to output given the output options
               */
              uint
              get_num_phases_to_output(Output_Options const & aOutputOptions);

              //------------------------------------------------------------------------------

              /*!
               * Setup clustering data
               */
              void
              setup_cell_clusters_for_output(
                      moris::mtk::Cell_Cluster_Input       & aCellClusterInput,
                      Output_Options                 const & aOutputOptions,
                      moris::Cell<Matrix<IdMat>>           & tCellIds);

              //------------------------------------------------------------------------------

              void
              setup_interface_side_cluster(
                      std::string                            aInterfaceSideLabelBase,
                      moris::mtk::Side_Cluster_Input       & aCellClusterInput,
                      Output_Options                 const & aOutputOptions,
                      moris::Cell<Matrix<IdMat>>           & tCellIdsandSideOrds,
                      moris::Cell<Matrix<DDRMat>>          & tParametricCoordinates);

              //------------------------------------------------------------------------------

              bool
              output_node(
                      moris::moris_index     aNodeIndex,
                      Output_Options const & aOutputOptions);

              //------------------------------------------------------------------------------

              /*
               * Prints the method of decomposition, type of background mesh,
               */
              void
              print_decompsition_preamble(Cell<enum Subdivision_Method> aMethods);

              //------------------------------------------------------------------------------

              moris::size_t
              determine_element_phase_index(
                      moris::size_t aRowIndex,
                      moris::Matrix< moris::IndexMat > const & aElementToNodeIndex);

              //------------------------------------------------------------------------------

              void
              collect_subphases_attached_to_facet_on_cell(
                      moris::moris_index         aCellIndex,
                      moris::moris_index         aFacetIndex,
                      Cell<moris::moris_index> & aCellSubphaseIndices,
                      Cell<moris::moris_index> & aCellSubphaseBulkIndices);
                      


    };
}

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
