/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Model.cpp
 *
 */

#include "cl_XTK_Model.hpp"

#include <utility>
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_XTK_Mesh_Cleanup.hpp"
#include "cl_XTK_Diagnostics.hpp"
// #include "cl_XTK_Contact_Sandbox.hpp"
#include "cl_XTK_Multigrid.hpp"
#include "fn_all_true.hpp"
#include "fn_unique.hpp"
#include "op_equal_equal.hpp"
#include "fn_sort.hpp"
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "fn_equal_to.hpp"
#include "fn_generate_element_to_element.hpp"
#include "fn_create_faces_from_element_to_node.hpp"
#include "fn_create_edges_from_element_to_node.hpp"
#include "HDF5_Tools.hpp"
#include "cl_MTK_Visualization_STK.hpp"
#include "cl_MTK_Cell_Info_Factory.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_TOL_Memory_Map.hpp"
#include "cl_Tracer.hpp"
#include "fn_stringify_matrix.hpp"
#include "cl_XTK_Basis_Processor.hpp"

using namespace moris;

namespace moris::xtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------

    Model::~Model()
    {
        delete mEnrichment;
        mEnrichment = nullptr;

        delete mGhostStabilization;
        mGhostStabilization = nullptr;

        for ( auto tIt : mEnrichedInterpMesh )
        {
            delete tIt;
        }

        mEnrichedInterpMesh.clear();

        for ( auto tIt : mEnrichedIntegMesh )
        {
            delete tIt;
        }

        mEnrichedIntegMesh.clear();
    }

    // ----------------------------------------------------------------------------------

    /*
     * using the general geometry engine
     */
    Model::Model(
            uint                            aModelDimension,
            moris::mtk::Interpolation_Mesh *aMeshData,
            moris::gen::Geometry_Engine    *aGeometryEngine,
            bool                            aLinkGeometryOnConstruction )
            : mModelDimension( aModelDimension )
            , mBackgroundMesh( aMeshData )
            , mCutMesh( this, mModelDimension )
            , mGeometryEngine( aGeometryEngine )
            , mEnrichment( nullptr )
            , mGhostStabilization( nullptr )
            , mEnrichedInterpMesh( 0, nullptr )
            , mEnrichedIntegMesh( 0, nullptr )
            , mConvertedToTet10s( false )
    {
        // flag this as a non-parameter list based run
        mParameterList.insert( "has_parameter_list", false );

        // flag that this field is false, this is for unit testing where enrichment is used
        mParameterList.insert( "write_cell_enrichments_levels", false );
    }

    // ----------------------------------------------------------------------------------

    Model::Model( moris::Parameter_List const &aParameterList )
            : mParameterList( aParameterList )
            , mModelDimension( UINT_MAX )
            , mEnrichment( nullptr )
            , mGhostStabilization( nullptr )
            , mEnrichedInterpMesh( 0, nullptr )
            , mEnrichedIntegMesh( 0, nullptr )
            , mConvertedToTet10s( false )
    {
        // flag this as a parameter list based run
        mParameterList.insert( "has_parameter_list", true );

        this->setup_diagnostics(
                mParameterList.get< bool >( "diagnostics" ),
                mParameterList.get< std::string >( "diagnostics_path" ),
                mParameterList.get< std::string >( "diagnostics_id" ) );
    }

    // ----------------------------------------------------------------------------------

    void
    Model::set_geometry_engine( moris::gen::Geometry_Engine *aGeometryEngine )
    {
        mGeometryEngine = aGeometryEngine;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::set_mtk_background_mesh( moris::mtk::Interpolation_Mesh *aMesh )
    {
        this->initialize( aMesh );

        mInitializeCalled = true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::set_input_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
    {
        mMTKInputPerformer = std::move( aMTKPerformer );
    }
    // ----------------------------------------------------------------------------------

    void
    Model::set_cut_ig_mesh( std::shared_ptr< Cut_Integration_Mesh > aCutIgMesh )
    {
        mCutIntegrationMesh = std::move( aCutIgMesh );
        mDecomposed         = true;
    }
    // ----------------------------------------------------------------------------------

    void
    Model::set_output_performer( std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer )
    {
        mMTKOutputPerformer = std::move( aMTKPerformer );
    }

    // ----------------------------------------------------------------------------------

    void
    Model::initialize( moris::mtk::Interpolation_Mesh *aMesh )
    {
        mBackgroundMesh     = aMesh;
        mModelDimension     = aMesh->get_spatial_dim();
        mCutMesh            = Cut_Mesh( this, mModelDimension );
        mEnrichment         = nullptr;
        mGhostStabilization = nullptr;
        mEnrichedInterpMesh = Vector< Enriched_Interpolation_Mesh * >( 0, nullptr );
        mEnrichedIntegMesh  = Vector< Enriched_Integration_Mesh  *>( 0, nullptr );
        mConvertedToTet10s  = false;
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::perform()
    {
        bool tReturn = this->perform_decomposition();
        this->perform_enrichment();
        return tReturn;
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::perform_decomposition()
    {
        Tracer tTracer( "XTK", "Overall", "Run" );

        mVerbose      = mParameterList.get< bool >( "verbose" );
        mVerboseLevel = mParameterList.get< uint >( "verbose_level" );

        if ( !mInitializeCalled )
        {
            MORIS_ERROR( mMTKInputPerformer != nullptr, "xtk::Model::perform(), mMTKInputPerformer not set!" );

            // FIXME hardcodes to mesh pair index 0
            moris::mtk::Interpolation_Mesh *tMesh = mMTKInputPerformer->get_interpolation_mesh( 0 );

            this->initialize( tMesh );
        }

        MORIS_ASSERT( this->has_parameter_list(), "Perform can only be called on a parameter list based XTK" );
        MORIS_ERROR( this->valid_parameters(), "Invalid parameters detected in XTK." );

        if ( mParameterList.get< std::string >( "probe_bg_cells" ) != "" )
        {
            Matrix< IdMat > tBgCellIds;
            string_to_mat( mParameterList.get< std::string >( "probe_bg_cells" ), tBgCellIds );
            print( tBgCellIds, "tBgCellIds" );
            this->probe_bg_cell( tBgCellIds );
        }

        // perform decomposition if requested in parameter list
        if ( mParameterList.get< bool >( "decompose" ) )
        {
            // get parameters specifying the integration mesh from the parameter list
            mTriangulateAll       = mParameterList.get< bool >( "triangulate_all" );
            mTriangulateAllInPost = mParameterList.get< bool >( "triangulate_all_in_post" );
            mOnlyGenerateXtkTemp  = mParameterList.get< bool >( "only_generate_xtk_temp" );
            mIgElementOrder       = mParameterList.get< uint >( "ig_element_order" );

            // get and store indices of B-spline meshes wrt. which information needs to be constructed
            moris::string_to_mat( mParameterList.get< std::string >( "enrich_mesh_indices" ), mBsplineMeshIndices );

            // check the enriched B-spline mesh indices
            for ( uint iBspMesh = 0; iBspMesh < mBsplineMeshIndices.numel(); iBspMesh++ )
            {
                MORIS_ERROR( mBsplineMeshIndices( iBspMesh ) == (moris_index)iBspMesh,
                        "xtk::Model::perform_decomposition() - B-spline mesh indices marked for "
                        "enrichment should be in the same order as B-spline meshes associated with "
                        "the Lagrange mesh in HMR. They should also start with the first one." );
            }

            // if ( mParameterList.get< bool >( "cleanup_cut_mesh" ) )
            // {
            //     mCleanupMesh = true;
            // }

            // generate list of decomposition methods to be performed to generate the XTK mesh
            Vector< enum Subdivision_Method > tSubdivisionMethods = this->get_subdivision_methods();

            // perform the decomposition methods specified
            bool tSuccess = this->decompose( tSubdivisionMethods );

            // Return false if cells are on different refinement levels
            if ( !tSuccess )
            {
                return false;
            }
        }

        // perform mesh cleanup if requested through parameter list
        if ( mParameterList.get< bool >( "cleanup_cut_mesh" ) )
        {
            // set flag as member variable
            mCleanupMesh = true;

            // cleanup the mesh
            Mesh_Cleanup tMeshCleanup( this, &mParameterList );
            tMeshCleanup.perform();
        }

        // return that the decomposition was performed successfully
        return true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::perform_enrichment()
    {
        if ( mParameterList.get< bool >( "enrich" ) )
        {
            // get rank of the interpolation basis (B-spline or Lagrange Element)
            // determines which basis functions to perform enrichment on
            mtk::EntityRank tBasisRank = mtk::get_entity_rank_from_str( mParameterList.get< std::string >( "basis_rank" ) );

            // get flag whether basis enrichments need to be sorted
            bool tSortBasisEnrichmentLevels = mParameterList.get< bool >( "sort_basis_enrichment_levels" );

            // perform the enrichment
            this->perform_basis_enrichment( tBasisRank, mBsplineMeshIndices, tSortBasisEnrichmentLevels, this->uses_SPG_based_enrichment() );

            // if high to low double side sets need to be created
            if ( mParameterList.get< bool >( "high_to_low_dbl_side_sets" ) )
            {
                // log this operation
                Tracer tTracer( "XTK", "Create high to low double side sets" );

                // ----------------------------------------------------
                // collect all the double sided side sets to be created

                // get the number of bulk-phases that exist
                uint tNumBulkPhases = mGeometryEngine->get_num_bulk_phase();

                // initialize the lists of double sided side sets to be created
                // Note: the number of side set combinations is ( n^2 - n ) where n is the number of bulk phases - we're only looking at the lower triangular part, so divide by 2
                uint                  tNumDblSideSetsToCreate = tNumBulkPhases * ( tNumBulkPhases - 1 ) / 2;
                Vector< moris_index > tLeaderPhaseIndices( tNumDblSideSetsToCreate, MORIS_INDEX_MAX );
                Vector< moris_index > tFollowerPhaseIndices( tNumDblSideSetsToCreate, MORIS_INDEX_MAX );

                // count number of double side sets to be created
                uint tNumDblSideSetsCreated = 0;

                // loop over all possible combinations of bulk phases and collect the combinations to be created
                for ( uint iLeaderPhaseIndex = 0; iLeaderPhaseIndex < tNumBulkPhases; iLeaderPhaseIndex++ )
                {
                    for ( uint iFollowerPhaseIndex = 0; iFollowerPhaseIndex < tNumBulkPhases; iFollowerPhaseIndex++ )
                    {
                        if ( iLeaderPhaseIndex > iFollowerPhaseIndex )
                        {
                            tLeaderPhaseIndices( tNumDblSideSetsCreated )   = (moris_index)iLeaderPhaseIndex;
                            tFollowerPhaseIndices( tNumDblSideSetsCreated ) = (moris_index)iFollowerPhaseIndex;
                            tNumDblSideSetsCreated++;
                        }
                    }
                }

                // check that the number of double sided side sets to be created is correct
                MORIS_ASSERT(
                        tNumDblSideSetsCreated == tNumDblSideSetsToCreate,
                        "xtk::Model::perform_enrichment() - "
                        "Number of high-to-low double sided side sets created is incorrect." );

                // create the double sided side sets
                mEnrichedIntegMesh( 0 )->create_dbl_sided_interface_sets( tLeaderPhaseIndices, tFollowerPhaseIndices );
            }

            // get index of B-spline meshes indices that will be unenriched later
            Matrix< IndexMat > tUnenrichedBsplineMeshIndices;
            moris::string_to_mat( mParameterList.get< std::string >( "unenriched_mesh_indices" ), tUnenrichedBsplineMeshIndices );

            this->perform_unenrichment( tUnenrichedBsplineMeshIndices );

        } // end: enrichment

        if ( mParameterList.get< bool >( "identify_hanging_nodes" ) )
        {
            this->perform_hanging_node_identification();
        }

        if ( mParameterList.get< bool >( "ghost_stab" ) )
        {

            // construct ghost using the new procedure if specifically requested by user
            if ( this->uses_SPG_based_enrichment() )
            {
                this->construct_face_oriented_ghost_penalization_cells_new();
            }
            else    // otherwise, just use old way
            {
                this->construct_face_oriented_ghost_penalization_cells();
            }

            if ( mParameterList.get< bool >( "visualize_ghost" ) )
            {
                // log/trace the creation of Ghost mesh sets for visualization
                Tracer tTracer( "XTK", "GhostStabilization", "Visualize" );

                // visualize ghost using the new procedure if specifically requested by user
                if ( this->uses_SPG_based_enrichment() )
                {
                    // get the B-spline mesh indices
                    Matrix< IndexMat > tBsplineMeshIndices;
                    moris::string_to_mat( mParameterList.get< std::string >( "enrich_mesh_indices" ), tBsplineMeshIndices );

                    // visualize ghost mesh sets for all B-spline meshes and bulk-phases
                    for ( moris::moris_index iBM = 0; iBM < (moris_index)tBsplineMeshIndices.numel(); iBM++ )
                    {
                        for ( moris::moris_index iBP = 0; iBP < (moris_index)mGeometryEngine->get_num_bulk_phase(); iBP++ )
                        {
                            mGhostStabilization->visualize_ghost_on_mesh_new( iBM, iBP );
                        }
                    }
                }
                else    // otherwise, just use the old way
                {
                    // visualize ghost mesh sets for every bulk-phase
                    for ( moris::moris_index i = 0; i < (moris_index)mGeometryEngine->get_num_bulk_phase(); i++ )
                    {
                        mGhostStabilization->visualize_ghost_on_mesh( i );
                    }
                }
            }
        }

        std::string tUnionBlockStr = mParameterList.get< std::string >( "union_blocks" );
        if ( !tUnionBlockStr.empty() )
        {
            // get the blocks to unionize
            Vector< Vector< std::string > > tUnionBlockCells;
            moris::string_to_vector_of_vectors( tUnionBlockStr, tUnionBlockCells );

            // Row based
            Matrix< IndexMat > tUnionBlockColors      = string_to_mat< IndexMat >( mParameterList.get< std::string >( "union_block_colors" ) );
            std::string        tUnionNewBlockNamesStr = mParameterList.get< std::string >( "union_block_names" );

            Vector< Vector< std::string > > tNewBlockNames;
            moris::string_to_vector_of_vectors( tUnionNewBlockNamesStr, tNewBlockNames );

            MORIS_ERROR( tUnionBlockCells.size() == tNewBlockNames.size(), "Dimension Mismatch in number of union operations for block" );
            MORIS_ERROR( tUnionBlockCells.size() == tUnionBlockColors.n_rows(), "Dimension Mismatch in number of union operations for block" );

            for ( uint iUnion = 0; iUnion < tUnionBlockCells.size(); iUnion++ )
            {
                this->get_enriched_integ_mesh( 0 ).create_union_block( tUnionBlockCells( iUnion ), tNewBlockNames( iUnion )( 0 ), tUnionBlockColors.get_row( iUnion ) );
            }

            // communicate new list of block sets
            this->get_enriched_integ_mesh( 0 ).communicate_sets_of_type( mtk::SetType::BULK );
        }

        std::string tUnionSideSetStr = mParameterList.get< std::string >( "union_side_sets" );
        if ( !tUnionSideSetStr.empty() )
        {
            // get the blocks to unionize
            Vector< Vector< std::string > > tUnionSideSetCells;
            moris::string_to_vector_of_vectors( tUnionSideSetStr, tUnionSideSetCells );

            // Row based
            Matrix< IndexMat > tUnionSideSetColors      = string_to_mat< IndexMat >( mParameterList.get< std::string >( "union_side_set_colors" ) );
            std::string        tUnionNewSideSetNamesStr = mParameterList.get< std::string >( "union_side_set_names" );

            Vector< Vector< std::string > > tNewSideSetNames;
            moris::string_to_vector_of_vectors( tUnionNewSideSetNamesStr, tNewSideSetNames );

            MORIS_ERROR( tUnionSideSetCells.size() == tNewSideSetNames.size(), "Dimension Mismatch in number of union operations for side set" );
            MORIS_ERROR( tUnionSideSetCells.size() == tUnionSideSetColors.n_rows(), "Dimension Mismatch in number of union operations for side set" );

            for ( uint iUnion = 0; iUnion < tUnionSideSetCells.size(); iUnion++ )
            {
                this->get_enriched_integ_mesh( 0 ).create_union_side_set( tUnionSideSetCells( iUnion ), tNewSideSetNames( iUnion )( 0 ), tUnionSideSetColors.get_row( iUnion ) );
            }

            // communicate new list of side sets
            this->get_enriched_integ_mesh( 0 ).communicate_sets_of_type( mtk::SetType::SIDESET );
        }

        std::string tDeactivatedBlockStr = mParameterList.get< std::string >( "deactivate_all_but_blocks" );
        if ( !tDeactivatedBlockStr.empty() )
        {
            // get the blocks to unionize
            Vector< Vector< std::string > > tBlocksToKeepStr;
            moris::string_to_vector_of_vectors( tDeactivatedBlockStr, tBlocksToKeepStr );

            MORIS_ERROR( tBlocksToKeepStr.size() == 1, "deactivate_all_but_block issue: This operation can only be performed on time" );

            this->get_enriched_integ_mesh( 0 ).deactivate_all_blocks_except_selected( tBlocksToKeepStr( 0 ) );
        }

        std::string tDeactivatedSideSetStr = mParameterList.get< std::string >( "deactivate_all_but_side_sets" );
        if ( !tDeactivatedSideSetStr.empty() )
        {
            // get the blocks to unionize
            Vector< Vector< std::string > > tSideSetsToKeepStr;
            moris::string_to_vector_of_vectors( tDeactivatedSideSetStr, tSideSetsToKeepStr );

            MORIS_ERROR( tSideSetsToKeepStr.size() == 1, "deactivate_all_side_sets_except_selected issue: This operation can only be performed on time" );

            this->get_enriched_integ_mesh( 0 ).deactivate_all_side_sets_except_selected( tSideSetsToKeepStr( 0 ) );
        }

        if ( mParameterList.get< bool >( "multigrid" ) )
        {
            this->construct_multigrid();
        }

        if ( mEnriched )
        {
            // if SPG based enrichment is used, construct the cluster groups (i.e. B-spline mesh clusters) here
            if ( this->uses_SPG_based_enrichment() )
            {
                mEnrichedIntegMesh( 0 )->setup_cluster_groups();
            }

            // get the reference to the enriched IP mesh
            xtk::Enriched_Interpolation_Mesh &tEnrInterpMesh = this->get_enriched_interp_mesh();

            // get the reference to the enriched IG mesh
            xtk::Enriched_Integration_Mesh &tEnrIntegMesh = this->get_enriched_integ_mesh();

            // set a mesh name for the XTK-mesh in the MTK output performer
            std::string tXTKMeshName = "XTKMesh";

            // place the pair in mesh manager
            mMTKOutputPerformer->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh, false, tXTKMeshName );

            // basis agglomeration
            if ( mParameterList.get< bool >( "activate_basis_agglomeration" ) )
            {
                // create a static object
                Basis_Processor tBasisProcessor = Basis_Processor( this );

                // perform the basis agglomeration
                tBasisProcessor.perform_basis_extention();
            }

            else if ( mParameterList.get< bool >( "activate_cell_agglomeration" ) )
            {
                // create a static object
                Basis_Processor tBasisProcessor = Basis_Processor( this );

                // perform the basis agglomeration
                tBasisProcessor.perform_cell_agglomeration();
            }

            // write the enriched ip mesh
            if ( mParameterList.get< bool >( "exodus_output_XTK_ip_mesh" ) )
            {
                mEnrichedInterpMesh( 0 )->write_mesh( &mParameterList );
            }

            // output (to console) a list of all blocks & sets in the XIGA model
            if ( mParameterList.get< bool >( "print_enriched_ig_mesh" ) )
            {
                tEnrIntegMesh.print();
            }

            // Visualize XTK (write xtk_temp.exo)
            if ( mParameterList.get< bool >( "exodus_output_XTK_ig_mesh" ) )
            {
                // log the visualization step
                Tracer tTracer( "XTK", "Overall", "Visualize" );

                // compute and write cluster measures to exo output if requested
                if ( mParameterList.get< bool >( "write_cluster_measures_to_exo" ) )
                {
                    // pre-compute the cluster measures for writing to xtk_temp
                    mEnrichedIntegMesh( 0 )->visualize_cluster_measures();

                    // pre-compute cluster measures on B-spline meshes (only possible for SPG based enrichment)
                    if ( this->uses_SPG_based_enrichment() )
                    {
                        bool tWriteBsplineClusterInfo = mParameterList.get< bool >( "write_bspline_cluster_info" );
                        mEnrichedIntegMesh( 0 )->visualize_cluster_group_measures( tWriteBsplineClusterInfo );
                    }
                }

                // write the xtk_temp.exo
                tEnrIntegMesh.write_mesh( &mParameterList );
            }

            // print the memory usage of XTK
            if ( mParameterList.get< bool >( "print_memory" ) )
            {
                moris::Memory_Map tXtkMM = this->get_memory_usage();
                tXtkMM.par_print( "XTK Model" );
            }

            // print summary of mesh size to console
            MORIS_LOG_SPEC( "All_IG_verts", sum_all( tEnrIntegMesh.get_num_entities( mtk::EntityRank::NODE ) ) );
            MORIS_LOG_SPEC( "All_IG_cells", sum_all( tEnrIntegMesh.get_num_entities( mtk::EntityRank::ELEMENT ) ) );
            MORIS_LOG_SPEC( "All_IP_verts", sum_all( tEnrInterpMesh.get_num_entities( mtk::EntityRank::NODE ) ) );
            MORIS_LOG_SPEC( "All_IP_cells", sum_all( tEnrInterpMesh.get_num_entities( mtk::EntityRank::ELEMENT ) ) );
            MORIS_LOG_SPEC( "My_IG_verts", tEnrIntegMesh.get_num_entities( mtk::EntityRank::NODE ) );
            MORIS_LOG_SPEC( "My_IG_cells", tEnrIntegMesh.get_num_entities( mtk::EntityRank::ELEMENT ) );
            MORIS_LOG_SPEC( "My_IP_verts", tEnrInterpMesh.get_num_entities( mtk::EntityRank::NODE ) );
            MORIS_LOG_SPEC( "My_IP_cells", tEnrInterpMesh.get_num_entities( mtk::EntityRank::ELEMENT ) );
        }
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::has_parameter_list()
    {
        return mParameterList.get< bool >( "has_parameter_list" );
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::valid_parameters()
    {
        bool tDecompose = mParameterList.get< bool >( "decompose" );
        bool tEnrich    = mParameterList.get< bool >( "enrich" );
        bool tGhost     = mParameterList.get< bool >( "ghost_stab" );
        bool tMultigrid = mParameterList.get< bool >( "multigrid" );

        if ( tEnrich == true )
        {
            MORIS_ERROR( tDecompose, "To perform basis enrichment, decomposition is also required." );
        }

        if ( tGhost == true )
        {
            MORIS_ERROR( tDecompose && tEnrich, "To perform ghost stabilization, decomposition and enrichment are also required." );
        }

        if ( tMultigrid == true )
        {
            MORIS_ERROR( tDecompose && tEnrich, "To perform multigrid, decomposition and enrichment are also required." );
        }

        return true;
    }

    // ----------------------------------------------------------------------------------

    Vector< enum Subdivision_Method >
    Model::get_subdivision_methods()
    {
        MORIS_ASSERT( this->has_parameter_list(), "Perform can only be called on a parameter list based XTK" );

        Vector< enum Subdivision_Method > tSubdivisionMethods;

        uint               tSpatialDimension = this->get_spatial_dim();
        mtk::Geometry_Type tBGCellTopo       = this->get_parent_cell_geometry();
        std::string        tDecompStr        = mParameterList.get< std::string >( "decomposition_type" );
        moris::lint        tOctreeRefLevel   = std::stoi( mParameterList.get< std::string >( "octree_refinement_level" ) );

        if ( tDecompStr.compare( "octree_only" ) == 0 )
        {
            return { Subdivision_Method::NC_OCTREE };
        }

        if ( tOctreeRefLevel > -1 )
        {
            tSubdivisionMethods.push_back( Subdivision_Method::NC_OCTREE );
        }

        // determine if we are going conformal or not
        bool tConformal = true;
        if ( tDecompStr.compare( "conformal" ) == 0 )
        {
            tConformal = true;
        }
        else if ( tDecompStr.compare( "nonconformal" ) == 0 )
        {
            tConformal = false;
        }
        else
        {
            MORIS_ERROR( 0, "Invalid decomposition_type provided. Recognized Options: Conformal and Non-conformal" );
        }

        if ( tSpatialDimension == 2 )
        {
            if ( tBGCellTopo == mtk::Geometry_Type::QUAD && tConformal )
            {
                Vector< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
                tSubdivisionMethods.append( tMethods );
            }
            else if ( tBGCellTopo == mtk::Geometry_Type::QUAD && !tConformal )
            {
                Vector< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4 };
                tSubdivisionMethods.append( tMethods );
            }
        }
        else if ( tSpatialDimension == 3 )
        {
            if ( tBGCellTopo == mtk::Geometry_Type::HEX && tConformal )
            {
                Vector< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
                tSubdivisionMethods.append( tMethods );
            }
            else if ( tBGCellTopo == mtk::Geometry_Type::HEX && !tConformal )
            {
                Vector< enum Subdivision_Method > tMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8 };
                tSubdivisionMethods.append( tMethods );
            }
            else if ( tBGCellTopo == mtk::Geometry_Type::HEX && tConformal )
            {
                Vector< enum Subdivision_Method > tMethods = { Subdivision_Method::C_HIERARCHY_TET4 };
                tSubdivisionMethods.append( tMethods );
            }
        }
        else
        {
            MORIS_ASSERT( 0, "Invalid spatial dimension" );
        }

        return tSubdivisionMethods;
    }

    // ----------------------------------------------------------------------------------
    // Decomposition Source code
    // ----------------------------------------------------------------------------------

    bool
    Model::decompose( const Vector< enum Subdivision_Method > &aMethods )
    {
        // log/trace the mesh decomposition
        Tracer tTracer( "XTK", "Decomposition", "Decompose" );

        // initialize IMG object for performing decomposition
        Integration_Mesh_Generator tIntegrationGenerator( this, aMethods );

        // perform decomposition
        mCutIntegrationMesh = tIntegrationGenerator.perform();

        // check that only Lagrange elements on the lowest refinement-level are cut / got decomposed
        mDecomposed = mCutIntegrationMesh->mSameLevelChildMeshes;
        if ( !mDecomposed )
        {
            MORIS_LOG_INFO( "xtk::Model::decompose() - Decomposition Failed. Not all child meshes on the same level." );
        }

        // print diagnostic data if requested by user
        if ( mDiagnostics )
        {
            if ( interpolated_coordinate_check( mCutIntegrationMesh.get() ) )
            {
                MORIS_LOG_INFO( "Interpolated Coordinate Check: Pass" );
            }
            else
            {
                MORIS_LOG_INFO( "Interpolated Coordinate Check: Fail" );
            }
        }

        // return successful decomposition
        return mDecomposed;
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::verify_successful_node_assignment( Decomposition_Data &aDecompData )
    {
        uint tNumUnsuccessful = 0;
        for ( uint i = 0; i < aDecompData.tNewNodeId.size(); i++ )
        {
            if ( aDecompData.tNewNodeId( i ) == MORIS_INDEX_MAX )
            {
                tNumUnsuccessful++;
            }
        }

        if ( tNumUnsuccessful > 0 )
        {
            std::cout << "There were " << tNumUnsuccessful << " bad nodes of " << aDecompData.tNewNodeId.size() << " total nodes" << '\n';
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------------

    Vector< std::string >
    Model::check_for_and_remove_internal_seacas_side_sets( Vector< std::string > &aSideSetNames )
    {
        for ( std::vector< std::string >::iterator iSet = aSideSetNames.begin(); iSet != aSideSetNames.end(); ++iSet )
        {
            if ( iSet->compare( "surface_1_quad4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }

            else if ( iSet->compare( "surface_2_quad4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_3_quad4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_4_quad4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_5_quad4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_6_quad4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_hex8_quad_1" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_hex8_quad_2" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_hex8_quad4_1" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_hex8_quad4_2" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_hex8_quad4_3" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
            else if ( iSet->compare( "surface_hex8_quad4_4" ) == 0 )
            {
                aSideSetNames.data().erase( iSet-- );
            }
        }

        return aSideSetNames;
    }

    // ----------------------------------------------------------------------------------
    // Enrichment Source code
    // ----------------------------------------------------------------------------------

    void
    Model::perform_basis_enrichment(
            mtk::EntityRank const &aBasisRank,
            moris_index const     &tBsplineMeshIndex,
            bool                   aSortBasisEnrichmentLevels,
            bool                   aUseSpgBasedEnrichment )
    {
        // log/trace the whole enrichment process
        Tracer tTracer( "XTK", "Enrichment", "Enrich" );

        // check that the mesh has already been decomposed
        MORIS_ERROR( mDecomposed, "Model::perform_basis_enrichment() - Prior to computing basis enrichment, the decomposition process must be called." );

        // allocate some new enriched interpolation and integration meshes
        mEnrichedIntegMesh.resize( tBsplineMeshIndex + 1, nullptr );
        mEnrichedInterpMesh.resize( tBsplineMeshIndex + 1, nullptr );

        this->perform_basis_enrichment_internal( aBasisRank, { { tBsplineMeshIndex } }, aSortBasisEnrichmentLevels, aUseSpgBasedEnrichment );

        if ( this->mDiagnostics )
        {
            mEnrichment->write_diagnostics();
        }

        // Change the enrichment flag
        mEnriched = true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::perform_basis_enrichment(
            mtk::EntityRank const    &aBasisRank,
            Matrix< IndexMat > const &aBsplineMeshIndices,
            bool                      aSortBasisEnrichmentLevels,
            bool                      aUseSpgBasedEnrichment )
    {
        Tracer tTracer( "XTK", "Enrichment" );

        MORIS_ERROR( mDecomposed, "Prior to computing basis enrichment, the decomposition process must be called" );

        // allocate some new enriched interpolation and integration meshes
        mEnrichedIntegMesh.resize( aBsplineMeshIndices.numel() + 1, nullptr );
        mEnrichedInterpMesh.resize( aBsplineMeshIndices.numel() + 1, nullptr );

        this->perform_basis_enrichment_internal( aBasisRank, aBsplineMeshIndices, aSortBasisEnrichmentLevels, aUseSpgBasedEnrichment );

        if ( mDiagnostics )
        {
            mEnrichment->write_diagnostics();
        }

        // Change the enrichment flag
        mEnriched = true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::perform_hanging_node_identification()
    {
        Tracer tTracer( "XTK", "Identify hanging nodes" );
        MORIS_ERROR( 0, "NEEDS IMPLEMENTING" );
        // MORIS_ERROR( mDecomposed, "Mesh needs to be decomposed prior to identifying hanging nodes" );

        // MORIS_ERROR( mEnriched, "Mesh needs to be enriched prior to identifying hanging nodes" );

        // // iterate through child meshes
        // for ( uint iCM = 0; iCM < mCutMesh.get_num_child_meshes(); iCM++ )
        // {
        //     // active child mesh
        //     Child_Mesh &tChildMesh = mCutMesh.get_child_mesh( iCM );

        //     // get the neighbors
        //     Matrix< IndexMat > tElementNeighbors = mBackgroundMesh->get_elements_connected_to_element_and_face_ind_loc_inds( tChildMesh.get_parent_element_index() );

        //     Vector< moris_index > tTransitionFacets;

        //     // iterate through neighbor
        //     for ( uint iN = 0; iN < tElementNeighbors.n_cols(); iN++ )
        //     {
        //         moris_index tNeighborCellIndex = tElementNeighbors( 0, iN );

        //         moris_index tSharedFaceIndex = tElementNeighbors( 1, iN );

        //         if ( !mBackgroundMesh.entity_has_children( tNeighborCellIndex, mtk::EntityRank::ELEMENT ) )
        //         {
        //             tTransitionFacets.push_back( tSharedFaceIndex );
        //         }
        //     }

        //     tChildMesh.identify_hanging_nodes( tTransitionFacets );
        // }
    }

    void
    Model::probe_bg_cell( Matrix< IndexMat > const &tBGCellIds )
    {
        Tracer tTracer( "XTK", "BG Cell Probe" );

        for ( uint i = 0; i < tBGCellIds.numel(); i++ )
        {
            Tracer tTracer( "XTK", "BG Cell Probe", "Cell Id " + std::to_string( tBGCellIds( i ) ) );

            moris_index             tIndex      = mBackgroundMesh->get_loc_entity_ind_from_entity_glb_id( tBGCellIds( i ), mtk::EntityRank::ELEMENT );
            mtk::Cell              &tCell       = mBackgroundMesh->get_mtk_cell( tIndex );
            Matrix< IndexMat >      tVertexIds  = tCell.get_vertex_ids();
            Vector< mtk::Vertex * > tVertexPtrs = tCell.get_vertex_pointers();

            MORIS_LOG_SPEC( "Cell Id", tBGCellIds( i ) );
            MORIS_LOG_SPEC( "Cell Index", tIndex );
            MORIS_LOG_SPEC( "Cell Owner", tCell.get_owner() );
            MORIS_LOG_SPEC( "Vertex Ids", ios::stringify_log( tVertexIds ) );
        }
    }

    // ----------------------------------------------------------------------------------
    Cut_Integration_Mesh *
    Model::get_cut_integration_mesh()
    {
        MORIS_ASSERT( mDecomposed,
                "Cannot get cut integration mesh prior to the decomposition strategy " );

        return mCutIntegrationMesh.get();
    }

    Enrichment const &
    Model::get_basis_enrichment()
    {
        MORIS_ASSERT( mEnriched,
                "Cannot get basis enrichment from an XTK model which has not called perform_basis_enrichment " );

        return *mEnrichment;
    }

    // ----------------------------------------------------------------------------------

    Enriched_Interpolation_Mesh &
    Model::get_enriched_interp_mesh( moris::moris_index aIndex )
    {
        MORIS_ASSERT( mEnriched,
                "Cannot get enriched interpolation mesh from an XTK model which has not called perform_basis_enrichment " );

        return *( mEnrichedInterpMesh( aIndex ) );
    }

    // ----------------------------------------------------------------------------------

    Enriched_Integration_Mesh &
    Model::get_enriched_integ_mesh( moris::moris_index aIndex )
    {
        MORIS_ASSERT( mEnriched,
                "Cannot get enriched integration mesh from an XTK model which has not called perform_basis_enrichment " );

        return *( mEnrichedIntegMesh( aIndex ) );
    }

    // ----------------------------------------------------------------------------------

    void
    Model::perform_basis_enrichment_internal(
            mtk::EntityRank const    &aBasisRank,
            Matrix< IndexMat > const &aBsplineMeshIndices,
            bool                      aSortBasisEnrichmentLevels,
            bool                      aUseSpgBasedEnrichment )
    {
        // initialize enrichment (ptr because of circular dependency)
        mEnrichment = new Enrichment(
                Enrichment_Method::USE_INTERPOLATION_CELL_BASIS,
                aBasisRank,
                aBsplineMeshIndices,
                mGeometryEngine->get_num_phases(),
                this,
                mBackgroundMesh,
                aSortBasisEnrichmentLevels,
                aUseSpgBasedEnrichment );

        // Set verbose flag to match XTK.
        mEnrichment->mVerbose = mVerbose;

        // perform the enrichment
        if ( aUseSpgBasedEnrichment )
        {
            // FIXME: this needs to go once the SPG based enrichment is validated
            mEnrichment->perform_enrichment_new();
        }
        else
        {
            mEnrichment->perform_enrichment();
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_face_oriented_ghost_penalization_cells()
    {
        Tracer tTracer( "XTK", "GhostStabilization", "Stabilize" );

        MORIS_ERROR( mDecomposed, "Mesh needs to be decomposed prior to calling ghost penalization" );

        MORIS_ERROR( !mGhost, "Ghost penalization has already been called" );

        mGhostStabilization = new Ghost_Stabilization( this );

        mGhostStabilization->setup_ghost_stabilization();

        mGhost = true;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_face_oriented_ghost_penalization_cells_new()
    {
        Tracer tTracer( "XTK", "Ghost Stabilization", "Construct Ghost Facets (new approach)" );

        MORIS_ERROR( mDecomposed, "Mesh needs to be decomposed prior to calling ghost penalization" );

        MORIS_ERROR( !mGhost, "Ghost penalization has already been called" );

        mGhostStabilization = new Ghost_Stabilization( this );

        mGhostStabilization->setup_ghost_stabilization_new();

        mGhost = true;
    }

    // ----------------------------------------------------------------------------------

    Ghost_Stabilization &
    Model::get_ghost_stabilization( moris::moris_index aIndex )
    {
        MORIS_ERROR( mGhost, "Ghost has not been constructed on this model." );
        return *mGhostStabilization;
    }

    // ----------------------------------------------------------------------------------

    void
    Model::construct_multigrid()
    {
        Tracer tTracer( "XTK", "Multigrid", "Run" );

        mMultigrid = std::make_shared< xtk::Multigrid >( this );

        mMultigrid->build_enriched_coeff_to_background_coeff_map();

        mMultigrid->create_fine_to_coarse_relationship();

        mMultigrid->create_coarse_to_fine_relationship();

        mMultigrid->create_coarse_to_fine_weights();

        std::string tName = "Enriched_bspline_1.exo";
        mMultigrid->build_basis_exodus_information( tName );
    }

    // ----------------------------------------------------------------------------------

    Cut_Mesh &
    Model::get_cut_mesh()
    {
        return mCutMesh;
    }

    // ----------------------------------------------------------------------------------

    Cut_Mesh const &
    Model::get_cut_mesh() const
    {
        return mCutMesh;
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Interpolation_Mesh &
    Model::get_background_mesh()
    {
        return *mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------------

    moris::mtk::Interpolation_Mesh const &
    Model::get_background_mesh() const
    {
        return *mBackgroundMesh;
    }

    // ----------------------------------------------------------------------------------

    moris::gen::Geometry_Engine *
    Model::get_geom_engine()
    {
        return mGeometryEngine;
    }

    // ----------------------------------------------------------------------------------

    Matrix< IndexMat >
    Model::get_Bspline_mesh_indices() const
    {
        return mBsplineMeshIndices;
    }

    // ----------------------------------------------------------------------------------

    bool
    Model::subphase_is_in_child_mesh( moris_index aSubphaseIndex )
    {
        MORIS_ERROR( 0, "TO REMOVE" );
        return false;
    }

    // ----------------------------------------------------------------------------------

    moris::Parameter_List &
    Model::get_parameter_list()
    {
        return mParameterList;
    }

    //------------------------------------------------------------------------------
    void
    Model::setup_diagnostics(
            bool               aDiagnostics,
            std::string const &aDiagnosticPath,
            std::string const &aDiagnosticLabel )
    {
        mDiagnostics = aDiagnostics;

        if ( mDiagnostics )
        {
            mDiagnosticPath = aDiagnosticPath;
            mDiagnosticId   = aDiagnosticLabel;

            MORIS_ERROR( !mDiagnosticPath.empty(), "If diagnostics are turned on, a diagnostics path must be specified" );
            if ( mDiagnosticId.empty() )
            {
                mDiagnosticId = "no_spec";
            }
        }
    }
    std::string
    Model::get_diagnostic_file_name( std::string const &aLabel ) const
    {
        MORIS_ASSERT( mDiagnostics, "Only callable with diagnostics on" );
        return mDiagnosticPath + "/id_" + mDiagnosticId + "_ps_" + std::to_string( moris::par_size() ) + "_pr_" + std::to_string( moris::par_rank() ) + "_" + aLabel + ".csv";
    }

    bool
    Model::triangulate_all()
    {
        return mTriangulateAll;
    }

    uint
    Model::ig_element_order()
    {
        return mIgElementOrder;
    }

    //------------------------------------------------------------------------------

    moris::mtk::Integration_Mesh *
    Model::get_output_mesh( Output_Options const &aOutputOptions )

    {
        MORIS_ERROR( 0, "Deprecated. (Removal in progress)" );

        return nullptr;
    }

    //------------------------------------------------------------------------------

    moris::uint
    Model::get_spatial_dim() const
    {
        return mModelDimension;
    }

    //------------------------------------------------------------------------------

    moris_index
    Model::get_cell_xtk_index( moris_id aCellId )
    {
        auto tIter = mCellGlbToLocalMap.find( aCellId );
        MORIS_ASSERT( tIter != mCellGlbToLocalMap.end(), "Id not in map" );
        return tIter->second;
    }

    //------------------------------------------------------------------------------

    Vector< Vector< moris_index > > const &
    Model::get_subphase_to_subphase()
    {
        MORIS_ERROR( 0, "Deprecated." );
        return mSubphaseToSubPhase;
    }

    //------------------------------------------------------------------------------

    Vector< Vector< moris_index > > const &
    Model::get_subphase_to_subphase_my_side_ords()
    {
        MORIS_ERROR( 0, "Deprecated." );
        return mSubphaseToSubPhaseMySideOrds;
    }

    //------------------------------------------------------------------------------

    Vector< Vector< moris_index > > const &
    Model::get_subphase_to_subphase_transition_loc()
    {
        MORIS_ERROR( 0, "Deprecated." );
        return mTransitionNeighborCellLocation;
    }

    //------------------------------------------------------------------------------

    Vector< Vector< moris_index > > const &
    Model::get_subphase_to_subphase_neighbor_side_ords()
    {
        MORIS_ERROR( 0, "Deprecated." );
        return mSubphaseToSubPhaseNeighborSideOrds;
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< Multigrid >
    Model::get_multigrid_ptr()
    {
        return mMultigrid;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Model::get_num_subphase_neighbors()
    {
        Vector< Vector< moris_index > > const &tSubPhaseToSubphase = this->get_subphase_to_subphase();
        Matrix< IndexMat >                     tSubphaseNumNeighbors( 1, tSubPhaseToSubphase.size() );
        for ( size_t iSP = 0; iSP < tSubPhaseToSubphase.size(); iSP++ )
        {
            tSubphaseNumNeighbors( iSP ) = tSubPhaseToSubphase( iSP ).size();
        }
        return tSubphaseNumNeighbors;
    }

    //------------------------------------------------------------------------------

    moris_id
    Model::get_subphase_id( moris_id aSubphaseIndex )
    {

        return mCutIntegrationMesh->get_subphase_id( aSubphaseIndex );
    }

    // -----------------------------------------------------------------------------

    moris_index
    Model::get_subphase_index( moris_id aSubphaseId )
    {
        return mCutIntegrationMesh->get_subphase_index( aSubphaseId );
    }

    //------------------------------------------------------------------------------

    moris_id
    Model::get_subphase_group_id(
            moris_id    aSubphaseGroupIndex,
            moris_index aBsplineMeshIndex )
    {
        return mCutIntegrationMesh->get_subphase_group_id( aSubphaseGroupIndex, aBsplineMeshIndex );
    }

    // -----------------------------------------------------------------------------

    moris_index
    Model::get_subphase_group_index(
            moris_id    aSubphaseGroupId,
            moris_index aBsplineMeshIndex )
    {
        return mCutIntegrationMesh->get_subphase_group_index( aSubphaseGroupId, aBsplineMeshIndex );
    }

    // -----------------------------------------------------------------------------

    std::string
    Model::get_global_T_matrix_output_file_name()
    {
        // get value from parameter list
        return mParameterList.get< std::string >( "global_T_matrix_output_file" );
    }

    std::string
    Model::get_elemental_T_matrix_output_file_name()
    {
        // get value from parameter list
        return mParameterList.get< std::string >( "elemental_T_matrix_output_file" );
    }

    std::string
    Model::get_MPC_output_file_name()
    {
        // get value from parameter list
        return mParameterList.get< std::string >( "MPC_output_file" );
    }

    // -----------------------------------------------------------------------------

    bool
    Model::only_generate_xtk_temp()
    {
        // indicate to kill workflow after XTK mesh output if requested by user
        if ( mOnlyGenerateXtkTemp )
        {
            return true;
        }
        else    // otherwise don't kill the workflow
        {
            return false;
        }
    }

    // -----------------------------------------------------------------------------

    bool
    Model::kill_workflow_flag()
    {
        // indicate to kill workflow if T-matrix output of full triangulation in post-processing of the cut IG mesh has been requested
        if ( this->get_global_T_matrix_output_file_name() != "" || this->get_elemental_T_matrix_output_file_name() != "" || mTriangulateAllInPost )
        {
            return true;
        }
        else    // otherwise don't kill the workflow
        {
            return false;
        }
    }

    //------------------------------------------------------------------------------

    moris::Memory_Map
    Model::get_memory_usage()
    {
        // memory map for model
        moris::Memory_Map tXTKModelMM;

        // member data that have memory maps
        moris::Memory_Map tCutMeshMM;
        moris::Memory_Map tBGMeshMM;
        moris::Memory_Map tEnrichmentMM;
        moris::Memory_Map tGhostMM;
        moris::Memory_Map tIgMeshMM;
        moris::Memory_Map tIpMeshMM;

        if ( mDecomposed )
        {
            tCutMeshMM = mCutMesh.get_memory_usage();
        }

        if ( mEnriched )
        {
            tEnrichmentMM = mEnrichment->get_memory_usage();
            tIgMeshMM     = this->get_enriched_integ_mesh().get_memory_usage();
            tIpMeshMM     = this->get_enriched_interp_mesh().get_memory_usage();
        }

        if ( mGhost )
        {
            tGhostMM = mGhostStabilization->get_memory_usage();
        }

        // make the sum of the cut mesh memory map the cut mesh memory
        tXTKModelMM.mMemoryMapData[ "Cut Mesh" ]                            = tCutMeshMM.sum();
        tXTKModelMM.mMemoryMapData[ "Enrichment" ]                          = tEnrichmentMM.sum();
        tXTKModelMM.mMemoryMapData[ "Enriched Ig Mesh" ]                    = tIgMeshMM.sum();
        tXTKModelMM.mMemoryMapData[ "Enriched Ip Mesh" ]                    = tIpMeshMM.sum();
        tXTKModelMM.mMemoryMapData[ "Ghost" ]                               = tGhostMM.sum();
        tXTKModelMM.mMemoryMapData[ "Background Mesh" ]                     = tBGMeshMM.sum();
        tXTKModelMM.mMemoryMapData[ "mElementToElement ptrs" ]              = moris::internal_capacity( mElementToElement );
        tXTKModelMM.mMemoryMapData[ "mElementToElement ptrs" ]              = moris::internal_capacity( mElementToElement );
        tXTKModelMM.mMemoryMapData[ "mSubphaseToSubPhase" ]                 = moris::internal_capacity( mSubphaseToSubPhase );
        tXTKModelMM.mMemoryMapData[ "mSubphaseToSubPhaseMySideOrds" ]       = moris::internal_capacity( mSubphaseToSubPhaseMySideOrds );
        tXTKModelMM.mMemoryMapData[ "mSubphaseToSubPhaseNeighborSideOrds" ] = moris::internal_capacity( mSubphaseToSubPhaseNeighborSideOrds );

        tIgMeshMM.par_print( "Ig Mesh" );
        return tXTKModelMM;
    }

    //------------------------------------------------------------------------------
    void
    Model::perform_unenrichment( Matrix< IndexMat > const &aUnenrichedBsplineMeshIndices )
    {
        // if there is any elements in the matrix
        if ( aUnenrichedBsplineMeshIndices.numel() )
        {
            Tracer tTracer( "XTK", "No-Type", "Unenrichment" );

            // set the mesh indices
            mEnrichedInterpMesh( 0 )->set_unenriched_mesh_indices( aUnenrichedBsplineMeshIndices );

            // override id and index of t-matrices
            mEnrichedInterpMesh( 0 )->override_vertex_enrichment_id_index();

            // override the required maps
            mEnrichedInterpMesh( 0 )->override_maps();
        }
    }

    //------------------------------------------------------------------------------

    bool
    Model::uses_SPG_based_enrichment()
    {
        if ( mParameterList.get< bool >( "has_parameter_list" ) )
        {
            if ( mParameterList.get< bool >( "use_SPG_based_enrichment" ) )
            {
                return true;
            }
        }

        return false;
    }

    //------------------------------------------------------------------------------

    bool
    Model::delete_xtk_after_generation()
    {
        // get the ouput from parameter
        return mParameterList.get< bool >( "delete_xtk_after_generation" );
    }

}    // namespace moris::xtk
