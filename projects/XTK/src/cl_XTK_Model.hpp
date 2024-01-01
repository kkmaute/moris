/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Model.hpp
 *
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

// XTKL: Container includes
#include "containers.hpp"

// XTKL: Tools includes
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Decomposition_Data.hpp"
#include "cl_Interpolaton.hpp"

#include "cl_XTK_Output_Options.hpp"

// Linalg Includes
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "fn_print.hpp"

#include "cl_Communication_Tools.hpp"

// Topology
// TODO: MOVE THESE WITH CUTTING METHODS SOMEWHERE ELSE
#include "cl_XTK_Topology.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Quad_4_Topology.hpp"
#include "cl_XTK_Hexahedron_8_Topology.hpp"

// general geometry engine class
#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_Param_List.hpp"

namespace xtk
{
    class Integration_Mesh_Generator;
    class Enrichment;
    class Enrichment_Parameters;
    class Ghost_Stabilization;
    class Enriched_Interpolation_Mesh;
    class Enriched_Integration_Mesh;
    class Ghost_Stabilization;
    class Multigrid;
    class Basis_Processor;
}    // namespace xtk

using namespace moris;

namespace xtk
{
    class Model
    {
        // friend classes
        friend class Integration_Mesh_Generator;
        friend class Cut_Integration_Mesh;
        friend class Enrichment;
        friend class Enriched_Interpolation_Mesh;
        friend class Enriched_Integration_Mesh;
        friend class Cut_Integration_Mesh;
        friend class Ghost_Stabilization;
        friend class Multigrid;
        friend class Basis_Processor;

        //--------------------------------------------------------------------------------

      protected:
        moris::ParameterList mParameterList;

        uint                                    mModelDimension;
        moris::mtk::Interpolation_Mesh*         mBackgroundMesh;
        Cut_Mesh                                mCutMesh;
        moris::ge::Geometry_Engine*             mGeometryEngine;
        std::shared_ptr< Cut_Integration_Mesh > mCutIntegrationMesh;
        Enrichment*                             mEnrichment;
        Ghost_Stabilization*                    mGhostStabilization;
        Vector< Enriched_Interpolation_Mesh* >    mEnrichedInterpMesh;
        Vector< Enriched_Integration_Mesh* >      mEnrichedIntegMesh;
        std::shared_ptr< xtk::Multigrid >       mMultigrid;
        Matrix< IndexMat >                      mBsplineMeshIndices;

        //--------------------------------------------------------------------------------

      private:
        // XTK Model State Flags
        bool mDecomposed        = false;    // Model has been decomposed
        bool mConvertedToTet10s = false;    // Model has been converted from tet4's to tet10s
        bool mMeshDataFinalized = false;
        bool mEnriched          = false;    // Model has been enriched
        bool mUnzipped          = false;    // Model has been unzipped
        bool mGhost             = false;    // Model has setup ghost stabilization

        // Flag to cleanup mesh at end of decomposition
        bool mTriangulateAll       = false;    // Triangulate all background cells
        bool mTriangulateAllInPost = false;    // Triangulate all background cells in post-processing of cut integration mesh (only for exo-mesh output)
        bool mOnlyGenerateXtkTemp  = false;    // Flag to kill the workflow right after generating the xtk_temp file. Good for phase assignment.
        bool mCleanupMesh          = false;    // Cleanup the mesh

        // polynomial order of IG elements, 1 is default, 2 or higher will invoke order elevation in decomposition
        uint mIgElementOrder = 1;    // Triangulate all background cells

        // diagnostic information
        bool        mDiagnostics    = false;
        std::string mDiagnosticPath = "";
        std::string mDiagnosticId   = "";

        // communication table
        moris::Matrix< IdMat > mCommunicationMap;

        // cell map
        std::map< moris_id, moris_index > mCellGlbToLocalMap;

        // element to element neighborhood
        Vector< Vector< moris::mtk::Cell* > > mElementToElement;
        Vector< Vector< moris_index > >       mElementToElementSideOrds;
        Vector< Vector< moris_index > >       mElementToElementNeighborSideOrds;
        Vector< Vector< moris_index > >       mSubphaseToSubPhase;
        Vector< Vector< moris_index > >       mSubphaseToSubPhaseMySideOrds;
        Vector< Vector< moris_index > >       mSubphaseToSubPhaseNeighborSideOrds;

        // in the case of a hierarchically refined mesh, there are transitions with hanging nodes
        // this data flags the transition from a large facet to a smaller facet. (this is trivial
        // for non hmr type meshes)
        Vector< Vector< moris_index > > mTransitionNeighborCellLocation;

        // local to global subphase map
        std::unordered_map< moris::moris_id, moris::moris_index > mGlobalToLocalSubphaseMap;

        std::shared_ptr< mtk::Mesh_Manager > mMTKInputPerformer  = nullptr;
        std::shared_ptr< mtk::Mesh_Manager > mMTKOutputPerformer = nullptr;

        bool mInitializeCalled = false;

        //--------------------------------------------------------------------------------

      public:
        // Public member functions/data
        bool        mVerbose      = false;
        moris::uint mVerboseLevel = 0;

        //--------------------------------------------------------------------------------
        // Initialization
        //--------------------------------------------------------------------------------
        Model(){};

        //--------------------------------------------------------------------------------
        /**
         * @brief Constructor for standalone use of XTK
         *
         */
        Model( uint                             aModelDimension,
                moris::mtk::Interpolation_Mesh* aMeshData,
                moris::ge::Geometry_Engine*     aGeometryEngine,
                bool                            aLinkGeometryOnConstruction = true );

        //--------------------------------------------------------------------------------

        /**
         * @brief Parameter list based XTK initialization
         */
        Model( moris::ParameterList const & aParameterList );

        //--------------------------------------------------------------------------------

        ~Model();

        //--------------------------------------------------------------------------------
        /*
         * @brief Set the pointer to the geometry engine
         */
        void
        set_geometry_engine( moris::ge::Geometry_Engine* aGeometryEngine );

        //--------------------------------------------------------------------------------
        /**
         * @brief Set the pointer to the background interpolation mesh
         */
        void
        set_mtk_background_mesh( moris::mtk::Interpolation_Mesh* aMesh );

        //--------------------------------------------------------------------------------

        /**
         * @brief Set the pointer to the output performer
         */
        void
        set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

        //--------------------------------------------------------------------------------

        /**
         * @brief Set the pointer to the input performer
         */
        void
        set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer );

        void
        set_cut_ig_mesh( std::shared_ptr< Cut_Integration_Mesh > aCutIgMesh );

        //--------------------------------------------------------------------------------

        /**
         * @brief Initialize data using the interpolation mesh
         */
        void
        initialize( moris::mtk::Interpolation_Mesh* aMesh );

        //--------------------------------------------------------------------------------
        // Operations
        //--------------------------------------------------------------------------------
        /**
         * @brief Workflow based perform call. Uses parameter list to decide operations to
         * perform
         */
        bool
        perform();

        bool
        perform_decomposition();

        void
        perform_enrichment();

        //--------------------------------------------------------------------------------
        /**
         * Decomposes a mesh to conform to a geometric interface
         * @param[in] aMethods Subdivision Methods
         */
        bool decompose( Vector< enum Subdivision_Method > aMethods );

        //--------------------------------------------------------------------------------
        /**
         * Remove child meshes that have an intersection but all cells are in the same bulk phase
         */
        void
        set_cleanup_cut_mesh_flag( bool aCleanupMesh );

        //--------------------------------------------------------------------------------
        /**
         * @brief get the memory usage of XTK
         */
        moris::Memory_Map
        get_memory_usage();

        // ----------------------------------------------------------------------------------

        /**
         * @brief Perform basis enrichment
         * @param[in] aBasisRank                 Entity rank of basis functions
         * @param[in] aMeshIndex                 Mesh index (default: 0)
         * @param[in] aSortBasisEnrichmentLevels Flag whether to sort basis enrichment levels
         *                                       (default: false)
         */
        void
        perform_basis_enrichment(
                mtk::EntityRank const & aBasisRank,
                moris_index const &     aBsplineMeshIndex          = 0,
                bool                    aSortBasisEnrichmentLevels = false,
                bool                    aUseSpgBasedEnrichment     = false );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Perform basis enrichment
         * @param[in] aBasisRank Entity rank of basis functions
         * @param[in] aMeshIndex List of mesh indices
         * @param[in] aSortBasisEnrichmentLevels Flag whether to sort basis enrichment levels
         *                                       (default: false)
         */
        void
        perform_basis_enrichment(
                mtk::EntityRank const &    aBasisRank,
                Matrix< IndexMat > const & aBsplineMeshIndices,
                bool                       aSortBasisEnrichmentLevels = false,
                bool                       aUseSpgBasedEnrichment     = false );

        // ----------------------------------------------------------------------------------

        void perform_hanging_node_identification();

        // ----------------------------------------------------------------------------------
        /**
         * @brief Probes and prints the information about a background cell
         */
        void
        probe_bg_cell( Matrix< IndexMat > const & tBGCellIds );

        // ----------------------------------------------------------------------------------

        Cut_Integration_Mesh*
        get_cut_integration_mesh();
        /**
         * @return Basis enrichment
         */
        Enrichment const &
        get_basis_enrichment();

        // ----------------------------------------------------------------------------------
        /**
         * @param[in] aIndex Interpolation Index
         * @return Enriched interpolation mesh
         */
        Enriched_Interpolation_Mesh&
        get_enriched_interp_mesh( moris::moris_index aIndex = 0 );

        // ----------------------------------------------------------------------------------
        /**
         * @param[in] aIndex Interpolation Index
         * @return Enriched integration mesh
         */
        Enriched_Integration_Mesh&
        get_enriched_integ_mesh( moris::moris_index aIndex = 0 );

        // ----------------------------------------------------------------------------------
        /**
         * @brief Constructs the face oriented ghost penalization
         */
        void
        construct_face_oriented_ghost_penalization_cells();

        void
        construct_face_oriented_ghost_penalization_cells_new();

        // ----------------------------------------------------------------------------------
        /**
         * @return Ghost stabilization
         */
        Ghost_Stabilization&
        get_ghost_stabilization( moris::moris_index aIndex = 0 );

        // ----------------------------------------------------------------------------------

        /**
         * @brief Construct multigrid information
         */
        void construct_multigrid();

        //--------------------------------------------------------------------------------
        // Member data access functions
        //--------------------------------------------------------------------------------
        /**
         * @return Cut mesh (collection of child meshes)
         */
        Cut_Mesh&
        get_cut_mesh();

        // ----------------------------------------------------------------------------------
        /**
         * @return Const cut mesh (collection of child meshes)
         */
        Cut_Mesh const &
        get_cut_mesh() const;

        // ----------------------------------------------------------------------------------
        /**
         * @return Background mesh
         */
        moris::mtk::Interpolation_Mesh&
        get_background_mesh();

        // ----------------------------------------------------------------------------------
        /**
         * @return Constant background mesh
         */
        moris::mtk::Interpolation_Mesh const &
        get_background_mesh() const;

        // ----------------------------------------------------------------------------------
        /**
         * @return Geometry engine pointer
         */
        moris::ge::Geometry_Engine*
        get_geom_engine();

        // ----------------------------------------------------------------------------------

        moris::Matrix< IndexMat >
        get_Bspline_mesh_indices() const;

        // ----------------------------------------------------------------------------------
        // Outputting functions
        //-----------------------------------------------------------------------------------
        /**
         * @brief Outputs the Mesh to a mesh data which can then be written to exodus files as desired.
         * @param[in] aOutputOptions Mesh output options
         */
        moris::mtk::Integration_Mesh*
        get_output_mesh( Output_Options const & aOutputOptions = Output_Options() );

        /**
         * get spatial dimension of model
         */
        moris::uint
        get_spatial_dim() const;

        //-----------------------------------------------------------------------------------

        /**
         * @returns element to element connectivity
         */
        Vector< Vector< moris::mtk::Cell* > > const &
        get_element_to_element()
        {
            return mElementToElement;
        };

        //-----------------------------------------------------------------------------------

        /**
         * @returns element to element my side ordinals
         */
        Vector< Vector< moris_index > > const &
        get_element_to_element_my_side_ords()
        {
            return mElementToElementSideOrds;
        };

        //-----------------------------------------------------------------------------------

        /**
         * @returns element to element neighbor side ordinals
         */
        Vector< Vector< moris_index > > const &
        get_element_to_element_neighbor_side_ords()
        {
            return mElementToElementNeighborSideOrds;
        };

        //-----------------------------------------------------------------------------------

        moris_index
        get_cell_xtk_index( moris_id aCellId );

        //-----------------------------------------------------------------------------------

        Vector< Vector< moris_index > > const &
        get_subphase_to_subphase();

        //-----------------------------------------------------------------------------------

        Vector< Vector< moris_index > > const &
        get_subphase_to_subphase_my_side_ords();

        //-----------------------------------------------------------------------------------

        Vector< Vector< moris_index > > const &
        get_subphase_to_subphase_neighbor_side_ords();

        //-----------------------------------------------------------------------------------
        moris::Matrix< moris::IndexMat >
        get_num_subphase_neighbors();

        //-----------------------------------------------------------------------------------

        Vector< Vector< moris_index > > const &
        get_subphase_to_subphase_transition_loc();

        //-----------------------------------------------------------------------------------

        bool
        subphase_is_in_child_mesh( moris_index aSubphaseIndex );

        //-----------------------------------------------------------------------------------

        moris::Matrix< moris::IndexMat >
        get_element_to_subphase();

        //-----------------------------------------------------------------------------------

        moris_id
        get_subphase_id( moris_id aSubphaseIndex );

        //-----------------------------------------------------------------------------------

        moris_index
        get_subphase_index( moris_id aSubphaseId );

        //-----------------------------------------------------------------------------------

        moris_id
        get_subphase_group_id(
                moris_id    aSubphaseGroupIndex,
                moris_index aBsplineMeshIndex );

        //-----------------------------------------------------------------------------------

        moris_index
        get_subphase_group_index(
                moris_id    aSubphaseGroupId,
                moris_index aBsplineMeshIndex );

        //-----------------------------------------------------------------------------------

        std::string
        get_global_T_matrix_output_file_name();

        std::string
        get_elemental_T_matrix_output_file_name();

        std::string
        get_MPC_output_file_name();

        //-----------------------------------------------------------------------------------

        bool
        only_generate_xtk_temp();

        bool
        kill_workflow_flag();

        //-----------------------------------------------------------------------------------

        std::shared_ptr< Multigrid > get_multigrid_ptr();

        //--------------------------------------------------------------------------------
        // Printing Functions
        //--------------------------------------------------------------------------------

        void
        print_cells();

        //------------------------------------------------------------------------------

        moris::ParameterList&
        get_parameter_list();

        //------------------------------------------------------------------------------

        void
        setup_diagnostics(
                bool                aDiagnostics,
                std::string const & aDiagnosticPath,
                std::string const & aDiagnosticLabel );

        //------------------------------------------------------------------------------

        std::string
        get_diagnostic_file_name( std::string const & aLabel ) const;

        //------------------------------------------------------------------------------

        /**
         * @brief Tells the integration mesh to either decompose all background cells or not
         *
         * @return true - decomposition mesh will triangulate/tessalate all cells
         * @return false
         */
        bool
        triangulate_all();

        //------------------------------------------------------------------------------

        /**
         * @brief Get the polynomial order of the TRI/TET elements
         *
         * @return uint - polynomial order of TRIs/TETs in the final integration mesh
         */
        uint
        ig_element_order();

        //------------------------------------------------------------------------------

        mtk::CellTopology
        get_parent_cell_topology() const
        {
            mtk::Cell& tCell = mBackgroundMesh->get_mtk_cell( 0 );
            return tCell.get_cell_info()->get_cell_topology();
        }

        //------------------------------------------------------------------------------

        mtk::Geometry_Type
        get_parent_cell_geometry() const
        {
            mtk::Cell& tCell = mBackgroundMesh->get_mtk_cell( 0 );
            return tCell.get_cell_info()->get_cell_geometry();
        }

        //------------------------------------------------------------------------------

        Matrix< IdMat >
        get_communication_table() const
        {
            return mCutIntegrationMesh->get_communication_table();
        }

        //------------------------------------------------------------------------------

        std::map< moris_id, moris_index >
        get_communication_map() const
        {
            return mCutIntegrationMesh->get_communication_map();
        }

        //------------------------------------------------------------------------------

        /**
         * @brief returns whether SPG based enrichment is used or not
         */
        bool
        uses_SPG_based_enrichment();

        //------------------------------------------------------------------------------

        /**
         * @brief function to activate deletion of the xtk and activate database in mtk
         *
         */
        bool
        delete_xtk_after_generation();

        // Private Functions

      private:
        //------------------------------------------------------------------------------

        // Functions that take parameter inputs and setup XTK inputs
        bool
        has_parameter_list();

        //------------------------------------------------------------------------------

        /**
         * Verifys provided parameter list
         */
        bool
        valid_parameters();

        //------------------------------------------------------------------------------

        /*
         * Using the parameter list, figure out the cell of subdivision methods
         */
        Vector< enum Subdivision_Method >
        get_subdivision_methods();

        //------------------------------------------------------------------------------

        void
        create_new_node_association_with_geometry( Decomposition_Data& tDecompData );

        //------------------------------------------------------------------------------

        /**
         * Verifies that the nodes in a decomposition have all been assigned node ids
         */
        bool
        verify_successful_node_assignment( Decomposition_Data& aDecompData );

        //------------------------------------------------------------------------------

      protected:
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
                moris::uint                       aGeomIndex,
                moris::Matrix< moris::IndexMat >& aActiveChildMeshIndices,
                moris::Matrix< moris::IndexMat >& aNewPairBool );

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
                Output_Options const &                 aOutputOptions,
                Vector< moris::Matrix< moris::IdMat > >& aCellIdsAndSideOrds,
                Vector< std::string >&                   aInterfaceSetNames );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

        Vector< std::string >
        assign_geometry_data_names();

        //------------------------------------------------------------------------------

        Vector< mtk::EntityRank >
        assign_geometry_data_field_ranks();

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

        /*!
         * Sets up background side sets for mesh output. Propogates the side set from
         * the background mesh to the output mesh
         */
        Vector< moris::mtk::MtkSideSetInfo >
        propogate_background_side_sets(
                Vector< moris::Matrix< IndexMat > >& aSideSetData,
                Output_Options const &                    aOutputOptions );

        //------------------------------------------------------------------------------

        /*!
         * Add a single side set from the background mesh
         */
        void
        propogate_background_side_set(
                std::string const &                        aSideSetName,
                moris::moris_index                         aNoChildIndex,
                moris::moris_index                         aChildIndex,
                Vector< moris::Matrix< IndexMat > >&  aElementIdsAndSideOrd,
                Vector< moris::mtk::MtkSideSetInfo >& aSideSetData,
                Output_Options const &                     aOutputOptions,
                bool                                       aOutputIndices );

        //------------------------------------------------------------------------------

        /*!
         * This function checks for the side sets which appear in a mesh that comes
         * from a SEACAS generated string. For example:
         * "generated:1x1x1|sideset:xXyYzZ"
         */
        Vector< std::string >
        check_for_and_remove_internal_seacas_side_sets( Vector< std::string >& aSideSetNames );

        //------------------------------------------------------------------------------

        /*!
         * Combine interface and non-interface blocks
         */
        Vector< moris::Matrix< moris::IdMat > >
        combine_interface_and_non_interface_blocks(
                Vector< moris::Matrix< moris::IdMat > >& tChildElementsByPhase,
                Vector< moris::Matrix< moris::IdMat > >& tNoChildElementsByPhase );

        //------------------------------------------------------------------------------

        /*
         * Returns the number of phases to output given the output options
         */
        uint
        get_num_phases_to_output( Output_Options const & aOutputOptions );

        //------------------------------------------------------------------------------

        /*!
         * Setup clustering data
         */
        void
        setup_cell_clusters_for_output(
                moris::mtk::Cell_Cluster_Input& aCellClusterInput,
                Output_Options const &          aOutputOptions,
                Vector< Matrix< IdMat > >& tCellIds );

        //------------------------------------------------------------------------------

        void
        setup_interface_side_cluster(
                std::string                      aInterfaceSideLabelBase,
                moris::mtk::Side_Cluster_Input&  aCellClusterInput,
                Output_Options const &           aOutputOptions,
                Vector< Matrix< IdMat > >&  tCellIdsandSideOrds,
                Vector< Matrix< DDRMat > >& tParametricCoordinates );

        //------------------------------------------------------------------------------

        bool
        output_node(
                moris::moris_index     aNodeIndex,
                Output_Options const & aOutputOptions );

        //------------------------------------------------------------------------------

        /*
         * Prints the method of decomposition, type of background mesh,
         */
        void
        print_decompsition_preamble( Vector< enum Subdivision_Method > aMethods );

        //------------------------------------------------------------------------------

        moris::size_t
        determine_element_phase_index(
                moris::size_t                            aRowIndex,
                moris::Matrix< moris::IndexMat > const & aElementToNodeIndex );

        //------------------------------------------------------------------------------

        void
        collect_subphases_attached_to_facet_on_cell(
                moris::moris_index          aCellIndex,
                moris::moris_index          aFacetOrdinal,
                Vector< moris::moris_index >& aCellSubphaseIndices,
                Vector< moris::moris_index >& aCellSubphaseBulkIndices,
                Vector< moris::moris_index >& aRepresentativeCellInd,
                Vector< moris::moris_index >& aRepresentativeCellSideOrdinal );

        //------------------------------------------------------------------------------

        /**
         * @brief perform unenrichment on an already enriched mesh by overrding the t-matrices id and index
         *
         * @param aUnenrichedBsplineMeshIndices
         */

        void
        perform_unenrichment( Matrix< IndexMat > const & aUnenrichedBsplineMeshIndices );

        //------------------------------------------------------------------------------

      private:
        //------------------------------------------------------------------------------
        // Enrichment computation functions
        //------------------------------------------------------------------------------

        void perform_basis_enrichment_internal(
                mtk::EntityRank const &    aBasisRank,
                Matrix< IndexMat > const & aBsplineMeshIndices,
                bool                       aSortBasisEnrichmentLevels,
                bool                       aUseSpgBasedEnrichment = false );
    };
}    // namespace xtk

#endif /* SRC_XTK_CL_XTK_MODEL_HPP_ */
