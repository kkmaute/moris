/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_HMR_XTK.cpp
 *
 */

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_HMR.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Checker.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_WRK_GEN_Performer.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "cl_WRK_perform_remeshing.hpp"

#include "fn_norm.hpp"

#include "cl_WRK_DataBase_Performer.hpp"
#include "cl_MIG.hpp"

#include "cl_WRK_Reinitialize_Performer.hpp"

namespace moris::wrk
{

    //--------------------------------------------------------------------------------------------------------------

    Workflow_HMR_XTK::Workflow_HMR_XTK( wrk::Performer_Manager* aPerformerManager )
            : Workflow( aPerformerManager )
    {
        // log & trace this function
        Tracer tTracer( "WRK", "HMR-XTK Workflow", "Create" );

        // Performer set for this workflow
        mPerformerManager->mHMRPerformer.resize( 1 );
        mPerformerManager->mGENPerformer.resize( 1 );
        mPerformerManager->mXTKPerformer.resize( 1 );
        mPerformerManager->mMTKPerformer.resize( 2 );
        mPerformerManager->mMDLPerformer.resize( 1 );
        mPerformerManager->mRemeshingMiniPerformer.resize( 1 );
        mPerformerManager->mDataBasePerformer.resize( 1 );

        // load parameter lists
        Module_Parameter_Lists tHMRParameterList   = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::HMR );
        Module_Parameter_Lists tGENParameterList   = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::GEN );
        Module_Parameter_Lists tMORISParameterList = aPerformerManager->mLibrary->get_parameters_for_module( Module_Type::MORISGENERAL );

        // create HMR performer
        mPerformerManager->mHMRPerformer( 0 ) = std::make_shared< hmr::HMR >( tHMRParameterList, mPerformerManager->mLibrary );

        // create MTK performer - will be used for HMR mesh
        mPerformerManager->mMTKPerformer( 0 ) = std::make_shared< mtk::Mesh_Manager >();

        // Create GE performer
        mPerformerManager->mGENPerformer( 0 ) = std::make_shared< gen::Geometry_Engine >(
                tGENParameterList,
                mPerformerManager->mLibrary );

        // create MTK performer - will be used for XTK mesh
        mPerformerManager->mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();

        // create MDL performer
        mPerformerManager->mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mPerformerManager->mLibrary, 0 );

        // Set MTK performer to HMR
        mPerformerManager->mHMRPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 0 ) );

        // Set performer to MDL
        mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

        if ( tMORISParameterList.size() != 0 )
        {
            if ( !tMORISParameterList( 0 ).empty() )
            {
                // create re-meshing performer - will be used for HMR mesh
                mPerformerManager->mRemeshingMiniPerformer( 0 ) =
                        std::make_shared< wrk::Remeshing_Mini_Performer >( tMORISParameterList( 0 )( 0 ), mPerformerManager->mLibrary );
            }

            if ( !tMORISParameterList( 2 ).empty() )
            {
                // allocate size
                mPerformerManager->mReinitializePerformer.resize( 1 );

                // construct the parameter list
                mPerformerManager->mReinitializePerformer( 0 ) = std::make_shared< wrk::Reinitialize_Performer >( mPerformerManager->mLibrary );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Workflow_HMR_XTK::initialize(
            Vector< real >&  aADVs,
            Vector< real >&  aLowerBounds,
            Vector< real >&  aUpperBounds,
            Matrix< IdMat >& aIjklIDs )
    {
        // Trace & log this function
        Tracer tTracer( "WRK", "HMR-XTK Workflow", "Initialize" );

        mInitializeOptimizationRestart = false;

        mIter = 0;

        Vector< std::shared_ptr< mtk::Field > > tFieldsIn;
        Vector< std::shared_ptr< mtk::Field > > tFieldsOut;

        // perform in first optimization iteration only
        if ( tIsFirstOptSolve )
        {
            // Step 1: HMR refinement -------------------------------------------------------------------

            // mPerformerManager->mHMRPerformer( 0 )->reset_HMR();

            // check whether HMR should be NOT initialized from restart file
            if ( not mPerformerManager->mHMRPerformer( 0 )->get_restarted_from_file() )
            {
                // uniform initial refinement
                mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement();

                // HMR refined by GE
                Refinement_Mini_Performer tRefinementPerformer;

                // get the GEN interface performer
                std::shared_ptr< gen::Geometry_Engine > tGeometryEngine = mPerformerManager->mGENPerformer( 0 );

                std::shared_ptr< Performer > tGenInterfacePerformer = std::make_shared< wrk::Gen_Performer >( tGeometryEngine );

                // perform initial refinement for GEN and user defined refinement functions
                tRefinementPerformer.perform_refinement_old( mPerformerManager->mHMRPerformer( 0 ), { tGenInterfacePerformer } );
            }

            // HMR finalize
            mPerformerManager->mHMRPerformer( 0 )->perform();

            // switch flag for first solve in optimization process
            tIsFirstOptSolve = false;
        }
        // after first solve in optimization process
        else
        {
            // get access to the performers
            std::shared_ptr< gen::Geometry_Engine >         tGeometryEngine     = mPerformerManager->mGENPerformer( 0 );
            std::shared_ptr< mdl::Model >                   tModelPerformer     = mPerformerManager->mMDLPerformer( 0 );
            std::shared_ptr< Remeshing_Mini_Performer >     tRemeshingPerformer = mPerformerManager->mRemeshingMiniPerformer( 0 );
            Vector< std::shared_ptr< hmr::HMR > >&          tHmrPerformers      = mPerformerManager->mHMRPerformer;
            Vector< std::shared_ptr< mtk::Mesh_Manager > >& tMtkPerformers      = mPerformerManager->mMTKPerformer;

            // get refinement fields from GEN and MDL performers
            tFieldsIn.append( tGeometryEngine->get_mtk_fields() );
            tFieldsIn.append( tModelPerformer->get_mtk_fields() );

            // check remeshing mini-performer has been built
            MORIS_ERROR(
                    tRemeshingPerformer,
                    "Workflow_HMR_XTK::initialize() - remeshing performer has not been built." );

            // refine meshes
            tRemeshingPerformer->perform_remeshing(
                    tFieldsIn,
                    tHmrPerformers,
                    tMtkPerformers,
                    tFieldsOut );

            // get access to the library
            std::shared_ptr< Library_IO > tLibrary = mPerformerManager->mLibrary;

            // re-initialize GEN
            Module_Parameter_Lists tGENParameterList = tLibrary->get_parameters_for_module( Module_Type::GEN );
            tGeometryEngine                       = std::make_shared< gen::Geometry_Engine >( tGENParameterList, tLibrary );
        }

        // Step 2: Initialize Level set field in GEN -----------------------------------------------
        {
            // retrieve the mesh pair
            const mtk::Mesh_Pair& tMeshPair = mPerformerManager->mMTKPerformer( 0 )->get_mesh_pair( 0 );

            // initialize GEN
            std::shared_ptr< gen::Geometry_Engine > tGeometryEngine = mPerformerManager->mGENPerformer( 0 );
            tGeometryEngine->distribute_advs( tMeshPair, tFieldsOut );

            // Get ADVs
            aADVs        = tGeometryEngine->get_advs();
            aLowerBounds = tGeometryEngine->get_lower_bounds();
            aUpperBounds = tGeometryEngine->get_upper_bounds();
            aIjklIDs     = tGeometryEngine->get_IjklIDs();
        }

    }    // end function: Workflow_HMR_XTK::initialize()

    //--------------------------------------------------------------------------------------------------------------

    Vector< real >
    Workflow_HMR_XTK::perform( Vector< real >& aNewADVs )
    {
        // get optimization iteration
        sint tOptIter = gLogger.get_iteration( "OPT", "Manager", "Perform" );

        // compute total optimization iteration accounting for iterations performed //
        // in previous instances of optimization algorithm
        tOptIter = tOptIter + mIter;

        // return vector of design criteria with NANs causing the optimization algorithm to restart
        if ( mIter >= mReinitializeIterInterval or ( uint ) tOptIter == mReinitializeIterFirst )
        {
            mInitializeOptimizationRestart = true;

            Vector< real > tVector( mNumCriteria, std::numeric_limits< real >::quiet_NaN() );

            return tVector;
        }

        // Stage *: Re-initialization of the adv field
        if ( mPerformerManager->mReinitializePerformer.size() > 0 )
        {
            // decide if the re-initialization would be required
            sint tReinitFreq = mPerformerManager->mReinitializePerformer( 0 )->get_reinitialization_frequency();

            if ( tOptIter > 0 and tOptIter % tReinitFreq == 0 )
            {
                // Set new advs in GE
                Tracer tTracer( "GEN", "Levelset", "Re-InitializeADVs" );

                mPerformerManager->mGENPerformer( 0 )->distribute_advs(
                        mPerformerManager->mMTKPerformer( 0 )->get_mesh_pair( 0 ),
                        mPerformerManager->mReinitializePerformer( 0 )->get_mtk_fields() );

                // get advs from GE and overwrite them
                aNewADVs = mPerformerManager->mGENPerformer( 0 )->get_advs();
            }
            else
            {
                mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );
            }
        }
        else
        {
            // Set new advs in GE
            mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );
        }

        // Stage 1: HMR refinement
        if ( mPerformerManager->mRemeshingMiniPerformer( 0 ) )
        {
            sint tReinitFreq = mPerformerManager->mRemeshingMiniPerformer( 0 )->get_reinitialization_frequency();

            if ( tOptIter > 0 and tOptIter % tReinitFreq == 0 )
            {

                // allocate auxiliary arrays
                Vector< real >  tADVs;
                Vector< real >  tLowerBounds;
                Vector< real >  tUpperBounds;
                Matrix< IdMat > tIjklIDs;

                // initialize HMR and GEN
                this->initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );

                // Set new advs in GE
                mPerformerManager->mGENPerformer( 0 )->set_advs( aNewADVs );
            }
        }

        // Stage 2: XTK -----------------------------------------------------------------------------

        // Read parameter list from shared object
        Module_Parameter_Lists tXTKParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::XTK );

        // Create XTK
        xtk::Model* tXTKPerformer = new xtk::Model( tXTKParameterList( 0 )( 0 ) );

        std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer = std::make_shared< mtk::Mesh_Manager >();

        // Set performers
        tXTKPerformer->set_geometry_engine( mPerformerManager->mGENPerformer( 0 ).get() );
        tXTKPerformer->set_input_performer( mPerformerManager->mMTKPerformer( 0 ) );
        tXTKPerformer->set_output_performer( tMTKPerformer );

        // Compute level set data in GEN
        // FIXME: HMR stores mesh with aura on 0
        mPerformerManager->mGENPerformer( 0 )->reset_mesh_information(
                mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

        // Output GEN fields, if requested
        mPerformerManager->mGENPerformer( 0 )->output_fields(
                mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ) );

        // mtk::Mesh_Checker tMeshCheckerHMR(
        //         0,
        //         mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh(0),
        //         mPerformerManager->mMTKPerformer( 0 )->get_integration_mesh(0));
        // tMeshCheckerHMR.perform();
        // tMeshCheckerHMR.print_diagnostics();

        bool tDeleteXTK = tXTKPerformer->delete_xtk_after_generation();
        // XTK perform - decompose - enrich - ghost - multigrid
        bool tFlag = tXTKPerformer->perform_decomposition();

        if ( not tFlag )
        {
            mInitializeOptimizationRestart = true;

            MORIS_ERROR( mNumCriteria != MORIS_UINT_MAX,
                    "Workflow_HMR_XTK::perform() problem with mNumCriteria. "
                    "This can happen if the xtk interface interfaces different refinement level in the first optimization iteration" );

            Vector< real > tVector( mNumCriteria, std::numeric_limits< real >::quiet_NaN() );

            if ( tDeleteXTK )
            {
                // delete the xtk
                delete tXTKPerformer;
            }

            return tVector;
        }

        // store whether the new ghost has been used
        bool tUseNewGhostSets = tXTKPerformer->uses_SPG_based_enrichment();

        // XTK perform - enrich - ghost - multigrid
        tXTKPerformer->perform_enrichment();

        if ( tDeleteXTK )
        {
            // construct the data base with the mtk performer from xtk
            mPerformerManager->mDataBasePerformer( 0 ) = std::make_shared< DataBase_Performer >( tMTKPerformer );

            // create the mtk performer that will hold the data base mesh pair and set it
            std::shared_ptr< mtk::Mesh_Manager > tMTKDataBasePerformer = std::make_shared< mtk::Mesh_Manager >();
            mPerformerManager->mDataBasePerformer( 0 )->set_output_performer( tMTKDataBasePerformer );

            // turn off the mesh check if no FEM-model will be constructed on the mesh
            mPerformerManager->mDataBasePerformer( 0 )->set_mesh_check( !tXTKPerformer->kill_workflow_flag() );

            // perform the mtk data base
            mPerformerManager->mDataBasePerformer( 0 )->perform();

            // set the mtk performer
            mPerformerManager->mMTKPerformer( 1 ) = tMTKDataBasePerformer;
        }

        // stop workflow if only pre-processing output is requested
        if ( tXTKPerformer->only_generate_xtk_temp() )
        {
            MORIS_LOG( "------------------------------------------------------------------------------" );
            MORIS_LOG( "Only output of the foreground mesh requested. Stopping workflow after XTK/MTK." );
            MORIS_LOG( "------------------------------------------------------------------------------" );

            // delete data that needs to be deleted explicitly to prevent memory leaks
            delete tXTKPerformer;
            mPerformerManager->mDataBasePerformer( 0 )->free_memory();

            // return empty vector for design criteria
            Vector< real > tVector( 1, std::numeric_limits< real >::quiet_NaN() );
            return tVector;
        }

        // output T-matrices and/or MPCs if requested
        this->output_T_matrices( tMTKPerformer, tXTKPerformer );

        // stop workflow if T-Matrices have been outputted
        if ( tXTKPerformer->kill_workflow_flag() )
        {
            MORIS_LOG( "----------------------------------------------------------------------------------------------------" );
            MORIS_LOG( "T-Matrix output or triangulation of all elements in post requested. Stopping workflow after XTK/MTK." );
            MORIS_LOG( "----------------------------------------------------------------------------------------------------" );

            // delete data that needs to be deleted explicitly to prevent memory leaks
            delete tXTKPerformer;
            mPerformerManager->mDataBasePerformer( 0 )->free_memory();

            // return empty vector for design criteria
            Vector< real > tVector( 1, std::numeric_limits< real >::quiet_NaN() );
            return tVector;
        }

        if ( tDeleteXTK )
        {
            // delete the xtk-performer
            delete tXTKPerformer;
        }
        else
        {
            // IMPORTANT!!! do not overwrite previous XTK and MTK performer before we know if this XTK performer triggers a restart.
            // otherwise the fem::field meshes are deleted and cannot be used anymore.
            mPerformerManager->mXTKPerformer( 0 ) = std::shared_ptr< xtk::Model >( tXTKPerformer );
            mPerformerManager->mMTKPerformer( 1 ) = tMTKPerformer;
        }

        // mtk::Mesh_Checker tMeshCheckerXTK(
        //         0,
        //         mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_interpolation_mesh(),
        //         mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh() );
        // tMeshCheckerXTK.perform();
        // tMeshCheckerXTK.print_diagnostics();

        // mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh()->save_MPC_to_hdf5();
        // mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh()->save_IG_global_T_matrix_to_file();

        // get the MIG parameter list
        Module_Parameter_Lists tMIGParameterList = mPerformerManager->mLibrary->get_parameters_for_module( Module_Type::MIG );

        // check if there are MIG parameters specified
        if ( tMIGParameterList.size() > 0 )
        {
            moris::mig::MIG tMIGPerformer = moris::mig::MIG(
                    mPerformerManager->mMTKPerformer( 1 ),
                    tMIGParameterList( 0 )( 0 ),
                    mPerformerManager->mGENPerformer( 0 ).get() );

            tMIGPerformer.perform();
        }

        if ( tDeleteXTK )
        {
            // free the memory and delete the unused data
            mPerformerManager->mDataBasePerformer( 0 )->free_memory();
        }

        mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

        // Assign PDVs
        mPerformerManager->mGENPerformer( 0 )->create_pdvs( mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ) );

        // Stage 3: MDL perform ---------------------------------------------------------------------

        mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface(
                mPerformerManager->mGENPerformer( 0 )->get_design_variable_interface() );

        // store whether the new ghost sets are being used
        mPerformerManager->mMDLPerformer( 0 )->set_use_new_ghost_mesh_sets( tUseNewGhostSets );

        mPerformerManager->mMDLPerformer( 0 )->initialize();

        mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

        // Build MDL components and solve
        mPerformerManager->mMDLPerformer( 0 )->perform();

        // perform mapping at this stage between solution field and adv field as some data will be deleted
        if ( mPerformerManager->mReinitializePerformer.size() > 0 )
        {
            // decide if we need to perform mapping
            if ( tOptIter % ( mPerformerManager->mReinitializePerformer( 0 )->get_reinitialization_frequency() ) ==    //
                    mPerformerManager->mReinitializePerformer( 0 )->get_reinitialization_frequency() - 1 )
            {
                // perform mapping of the solution to the adv
                mPerformerManager->mReinitializePerformer( 0 )->perform( mPerformerManager->mHMRPerformer,
                        mPerformerManager->mGENPerformer,
                        mPerformerManager->mMTKPerformer,
                        mPerformerManager->mMDLPerformer );
            }
        }

        // evaluate IQIs
        Vector< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

        // get number of design criteria
        mNumCriteria = tVal.size();

        // Communicate IQIs
        for ( uint iIQIIndex = 0; iIQIIndex < mNumCriteria; iIQIIndex++ )
        {
            tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
        }

        // build vector of design criteria
        Vector< real > tVector( mNumCriteria, 0.0 );
        for ( uint iCriteriaIndex = 0; iCriteriaIndex < mNumCriteria; iCriteriaIndex++ )
        {
            tVector( iCriteriaIndex ) = tVal( iCriteriaIndex )( 0 );
        }

        // increment optimization iteration counter
        mIter++;

        // return vector of design criteria
        return tVector;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Workflow_HMR_XTK::compute_dcriteria_dadv()
    {
        mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

        mPerformerManager->mMDLPerformer( 0 )->perform( 1 );

        Matrix< DDRMat > tDCriteriaDAdv = mPerformerManager->mGENPerformer( 0 )->get_dcriteria_dadv();

        if ( par_rank() == 0 )
        {
            MORIS_LOG_INFO( "--------------------------------------------------------------------------------" );
            MORIS_LOG_INFO( "Gradients of design criteria wrt ADVs:" );

            for ( uint i = 0; i < tDCriteriaDAdv.n_rows(); ++i )
            {
                Matrix< DDRMat > tdIQIdAdv = tDCriteriaDAdv.get_row( i );

                auto tItrMin = std::min_element( tdIQIdAdv.data(), tdIQIdAdv.data() + tdIQIdAdv.numel() );
                auto tIndMin = std::distance( tdIQIdAdv.data(), tItrMin );

                auto tItrMax = std::max_element( tdIQIdAdv.data(), tdIQIdAdv.data() + tdIQIdAdv.numel() );
                auto tIndMax = std::distance( tdIQIdAdv.data(), tItrMax );

                MORIS_LOG_INFO( "Criteria(%i): norm = %e   min = %e  (index = %li)   max = %e  (index = %li)",
                        i,
                        norm( tdIQIdAdv ),
                        tdIQIdAdv.min(),
                        tIndMin,
                        tdIQIdAdv.max(),
                        tIndMax );
            }

            MORIS_LOG_INFO( "--------------------------------------------------------------------------------" );
        }

        mPerformerManager->mMDLPerformer( 0 )->free_memory();

        return tDCriteriaDAdv;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Workflow_HMR_XTK::output_T_matrices(
            const std::shared_ptr< mtk::Mesh_Manager >& aMTKPerformer,
            xtk::Model* const &                         aXTKPerformer )
    {
        // Output T-matrices if requested
        std::string tTmatrixFileName = aXTKPerformer->get_global_T_matrix_output_file_name();
        if ( tTmatrixFileName != "" )
        {
            mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_IG_global_T_matrix_to_file( tTmatrixFileName );

            // return flag stopping the workflow after the T-Matrix output
            return true;
        }

        // Output T-matrices if requested
        std::string tElementalTmatrixFileName = aXTKPerformer->get_elemental_T_matrix_output_file_name();
        if ( tElementalTmatrixFileName != "" )
        {
            uint tNumBsplineMeshes = mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_interpolation_mesh()->get_num_interpolations();
            mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_elemental_T_matrices_to_file( tElementalTmatrixFileName, tNumBsplineMeshes );

            // return flag stopping the workflow after the T-Matrix output
            return true;
        }

        // Output MPCs if requested
        std::string tMpcFileName = aXTKPerformer->get_MPC_output_file_name();
        if ( tMpcFileName != "" )
        {
            mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_MPC_to_hdf5( tMpcFileName );
        }

        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::wrk
