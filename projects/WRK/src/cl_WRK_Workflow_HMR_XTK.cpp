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
namespace moris
{
    namespace wrk
    {

        //------------------------------------------------------------------------------

        // Parameter function
        typedef void ( *Parameter_Function )( moris::Cell< moris::Cell< moris::ParameterList > >& aParameterList );

        //--------------------------------------------------------------------------------------------------------------

        Workflow_HMR_XTK::Workflow_HMR_XTK( wrk::Performer_Manager* aPerformerManager )
                : Workflow( aPerformerManager )
        {
            // Performer set for this workflow
            mPerformerManager->mHMRPerformer.resize( 1 );
            mPerformerManager->mGENPerformer.resize( 1 );
            // mPerformerManager->mXTKPerformer.resize( 1 );
            mPerformerManager->mMTKPerformer.resize( 2 );
            mPerformerManager->mMDLPerformer.resize( 1 );
            mPerformerManager->mRemeshingMiniPerformer.resize( 1 );

            // load the HMR parameter list
            std::string tHMRString = "HMRParameterList";

            Parameter_Function tHMRParameterListFunc =    //
                    mPerformerManager->mLibrary->load_function< Parameter_Function >( tHMRString );

            moris::Cell< moris::Cell< ParameterList > > tHMRParameterList;
            tHMRParameterListFunc( tHMRParameterList );

            // load the GEN parameter list
            std::string tGENString = "GENParameterList";

            Parameter_Function tGENParameterListFunc =    //
                    mPerformerManager->mLibrary->load_function< Parameter_Function >( tGENString );

            moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
            tGENParameterListFunc( tGENParameterList );

            // load the MORIS General parameter list
            std::string tMORISString = "MORISGENERALParameterList";

            Parameter_Function tMORISParameterListFunc =    //
                    mPerformerManager->mLibrary->load_function< Parameter_Function >( tMORISString );

            moris::Cell< moris::Cell< ParameterList > > tMORISParameterList;
            tMORISParameterListFunc( tMORISParameterList );

            // create HMR performer
            mPerformerManager->mHMRPerformer( 0 ) =    //
                    std::make_shared< hmr::HMR >( tHMRParameterList( 0 )( 0 ), mPerformerManager->mLibrary );

            // create MTK performer - will be used for HMR mesh
            mPerformerManager->mMTKPerformer( 0 ) = std::make_shared< mtk::Mesh_Manager >();

            // Create GE performer
            mPerformerManager->mGENPerformer( 0 ) = std::make_shared< ge::Geometry_Engine >(
                    tGENParameterList,
                    mPerformerManager->mLibrary );

            // create MTK performer - will be used for XTK mesh
            mPerformerManager->mMTKPerformer( 1 ) = std::make_shared< mtk::Mesh_Manager >();

            // create MDL performer
            mPerformerManager->mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mPerformerManager->mLibrary, 0 );

            // Set MTK performer to HMR
            mPerformerManager->mHMRPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 0 ) );

            // Set HMR performer to MTK performer
            mPerformerManager->mMTKPerformer( 0 )->set_performer( mPerformerManager->mHMRPerformer( 0 ) );

            // Set performer to MDL
            mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

            if ( !tMORISParameterList.empty() )
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
                Matrix< DDRMat >& aADVs,
                Matrix< DDRMat >& aLowerBounds,
                Matrix< DDRMat >& aUpperBounds,
                Matrix< IdMat  >& aIjklIDs )
        {
            mInitializeOptimizationRestart = false;

            mIter = 0;

            moris::Cell< std::shared_ptr< mtk::Field > > tFieldsIn;
            moris::Cell< std::shared_ptr< mtk::Field > > tFieldsOut;

            // perform in first optimization iteration only
            if ( tIsFirstOptSolve )
            {
                // Stage 1: HMR refinement -------------------------------------------------------------------

                // Trace HMR
                Tracer tTracer( "HMR", "HMRmesh", "Create" );

                // mPerformerManager->mHMRPerformer( 0 )->reset_HMR();

                // check whether HMR should be NOT initialized from restart file
                if ( not mPerformerManager->mHMRPerformer( 0 )->get_restarted_from_file() )
                {
                    // uniform initial refinement
                    mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement();

                    // HMR refined by GE
                    Refinement_Mini_Performer tRefinementPerfomer;

                    // GEN interface performer
                    std::shared_ptr< Performer > tGenPerformer =
                            std::make_shared< wrk::Gen_Performer >( mPerformerManager->mGENPerformer( 0 ) );

                    tRefinementPerfomer.perform_refinement_old( mPerformerManager->mHMRPerformer( 0 ), { tGenPerformer } );
                }

                // HMR finalize
                mPerformerManager->mHMRPerformer( 0 )->perform();

                // switch flag for first solve in optimization process
                tIsFirstOptSolve = false;
            }
            // after first solve in optimization process
            else
            {
                // get refinement fields from GEN and MDL performers
                tFieldsIn.append( mPerformerManager->mGENPerformer( 0 )->get_mtk_fields() );
                tFieldsIn.append( mPerformerManager->mMDLPerformer( 0 )->get_mtk_fields() );

                // check remeshing mini-performer has been built
                MORIS_ERROR( mPerformerManager->mRemeshingMiniPerformer( 0 ),
                        " Workflow_HMR_XTK::initialize - remeshing performer has not been built." );

                // refine meshes
                mPerformerManager->mRemeshingMiniPerformer( 0 )->perform_remeshing(
                        tFieldsIn,
                        mPerformerManager->mHMRPerformer,
                        mPerformerManager->mMTKPerformer,
                        tFieldsOut );

                // create new GE performer
                std::string tGENString = "GENParameterList";

                Parameter_Function tGENParameterListFunc =    //
                        mPerformerManager->mLibrary->load_function< Parameter_Function >( tGENString );

                moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
                tGENParameterListFunc( tGENParameterList );

                mPerformerManager->mGENPerformer( 0 ) =
                        std::make_shared< ge::Geometry_Engine >( tGENParameterList, mPerformerManager->mLibrary );
            }

            // Stage 2: Initialize Level set field in GEN -----------------------------------------------
            {
                // Trace GEN
                Tracer tTracer( "GEN", "Levelset", "InitializeADVs" );

                mPerformerManager->mGENPerformer( 0 )->distribute_advs(
                        mPerformerManager->mMTKPerformer( 0 )->get_mesh_pair( 0 ),
                        tFieldsOut );

                // Get ADVs
                aADVs        = mPerformerManager->mGENPerformer( 0 )->get_advs();
                aLowerBounds = mPerformerManager->mGENPerformer( 0 )->get_lower_bounds();
                aUpperBounds = mPerformerManager->mGENPerformer( 0 )->get_upper_bounds();
                aIjklIDs     = mPerformerManager->mGENPerformer( 0 )->get_IjklIDs();
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Workflow_HMR_XTK::perform( Matrix< DDRMat >& aNewADVs )
        {
            // get optimization iteration
            sint tOptIter = gLogger.get_iteration(
                    "OptimizationManager",
                    LOGGER_ARBITRARY_DESCRIPTOR,
                    LOGGER_ARBITRARY_DESCRIPTOR );

            // compute total optimization iteration accounting for iterations performed //
            // in previous instances of optimization algorithm
            tOptIter = tOptIter + mIter;

            // return vector of design criteria with NANs causing the optimization algorithm to restart
            if ( mIter >= mReinitializeIterIntervall or ( uint ) tOptIter == mReinitializeIterFirst )
            {
                mInitializeOptimizationRestart = true;

                moris::Matrix< DDRMat > tMat( mNumCriterias, 1, std::numeric_limits< real >::quiet_NaN() );

                return tMat;
            }

            // Stage *: Reinitialization of the adv field
            if ( mPerformerManager->mReinitializePerformer.size() > 0 )
            {
                // decide if the reinitialization would be required
                if ( tOptIter > 0 and tOptIter % ( mPerformerManager->mReinitializePerformer( 0 )->get_reinitialization_frequency() ) == 0 )
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
            if ( false && tOptIter > 0 )
            {
                Matrix< DDRMat > tLowerBounds = mPerformerManager->mGENPerformer( 0 )->get_lower_bounds();
                Matrix< DDRMat > tUpperBounds = mPerformerManager->mGENPerformer( 0 )->get_upper_bounds();

                this->initialize( aNewADVs, tLowerBounds, tUpperBounds );
            }

            // Stage 2: XTK -----------------------------------------------------------------------------
            // Read parameter list from shared object
            Parameter_Function tXTKParameterListFunc =    //
                    mPerformerManager->mLibrary->load_function< Parameter_Function >( "XTKParameterList" );

            moris::Cell< moris::Cell< ParameterList > > tXTKParameterList;
            tXTKParameterListFunc( tXTKParameterList );

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

            //             mtk::Mesh_Checker tMeshCheckerHMR(
            //                     0,
            //                     mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh(0),
            //                     mPerformerManager->mMTKPerformer( 0 )->get_integration_mesh(0));
            //             tMeshCheckerHMR.perform();
            //             tMeshCheckerHMR.print_diagnostics();

            // XTK perform - decompose - enrich - ghost - multigrid
            bool tFlag = tXTKPerformer->perform_decomposition();

            if ( not tFlag )
            {
                mInitializeOptimizationRestart = true;

                MORIS_ERROR( mNumCriterias != MORIS_UINT_MAX,
                        "Workflow_HMR_XTK::perform() problem with mNumCriterias. "
                        "This can happen if the xtk interface interfaces different refinement level in the first optimization iteration" );

                moris::Matrix< DDRMat > tMat( mNumCriterias, 1, std::numeric_limits< real >::quiet_NaN() );

                // delete the xtk
                delete tXTKPerformer;

                return tMat;
            }

            // XTK perform - enrich - ghost - multigrid
            tXTKPerformer->perform_enrichment();

            // constrcut the data base with the mtk performer from xtk
            DataBase_Performer tDataBasePerformer = DataBase_Performer( tMTKPerformer );

            // create the mtk performer that will hold the data base mesh pair and set it
            std::shared_ptr< mtk::Mesh_Manager > tMTKDataBasePerformer = std::make_shared< mtk::Mesh_Manager >();
            tDataBasePerformer.set_output_performer( tMTKDataBasePerformer );

            // turn off the mesh check if no FEM-model will be constructed on the mesh
            tDataBasePerformer.set_mesh_check( !tXTKPerformer->kill_workflow_flag() );

            // perform the mtk data base
            tDataBasePerformer.perform();

            // set the mtk performer
            mPerformerManager->mMTKPerformer( 1 ) = tMTKDataBasePerformer;

            // output T-matrices and/or MPCs if requested
            this->output_T_matrices( tMTKPerformer, tXTKPerformer );

            // stop workflow if T-Matrices have been outputted
            if ( tXTKPerformer->kill_workflow_flag() )
            {
                MORIS_LOG( "----------------------------------------------------------------------------------------------------" );
                MORIS_LOG( "T-Matrix output or triangulation of all elements in post requested. Stopping workflow after XTK/MTK." );
                MORIS_LOG( "----------------------------------------------------------------------------------------------------" );
                moris::Matrix< DDRMat > tMat( 1, 1, std::numeric_limits< real >::quiet_NaN() );
                return tMat;
            }

            // delete the xtk-performer
            delete tXTKPerformer;

            // IMPORTANT!!! do not overwrite previous XTK  and MTK performer before we know if this XTK performer triggers a restart.
            // otherwise the fem::field meshes are deleted and cannot be used anymore.
            // mPerformerManager->mXTKPerformer( 0 ) = std::shared_ptr<xtk::Model> (tXTKPerformer);
            // mPerformerManager->mMTKPerformer( 1 ) = tMTKPerformer;

            //            mtk::Mesh_Checker tMeshCheckerXTK(
            //                    0,
            //                    mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_interpolation_mesh(),
            //                    mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh());
            //            tMeshCheckerXTK.perform();
            //            tMeshCheckerXTK.print_diagnostics();

            // mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh()->save_MPC_to_hdf5();
            // mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair(0).get_integration_mesh()->save_IG_node_TMatrices_to_file();


            Parameter_Function tMIGParameterListFunc =
                    mPerformerManager->mLibrary->load_function< Parameter_Function >( "MIGParameterList", false );

            if ( tMIGParameterListFunc )
            {
                moris::Cell< moris::Cell< ParameterList > > tMIGParameterList;
                tMIGParameterListFunc( tMIGParameterList );

                moris::mig::MIG tMIGPerformer = moris::mig::MIG(
                        mPerformerManager->mMTKPerformer( 1 ),
                        tMIGParameterList( 0 )( 0 ),
                        mPerformerManager->mGENPerformer( 0 ).get() );

                tMIGPerformer.perform();
            }

            // free the memory and delete the unused data
            tDataBasePerformer.free_memory();

            mPerformerManager->mMDLPerformer( 0 )->set_performer( mPerformerManager->mMTKPerformer( 1 ) );

            // Assign PDVs
            mPerformerManager->mGENPerformer( 0 )->create_pdvs( mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ) );

            // Stage 3: MDL perform ---------------------------------------------------------------------

            mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface(
                    mPerformerManager->mGENPerformer( 0 )->get_design_variable_interface() );

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
            moris::Cell< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

            // get number of design criteria
            mNumCriterias = tVal.size();

            // Communicate IQIs
            for ( uint iIQIIndex = 0; iIQIIndex < mNumCriterias; iIQIIndex++ )
            {
                tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
            }

            // build vector of design criteria
            moris::Matrix< DDRMat > tMat( mNumCriterias, 1, 0.0 );

            for ( uint Ik = 0; Ik < mNumCriterias; Ik++ )
            {
                tMat( Ik ) = tVal( Ik )( 0 );
            }

            // increment optimization iteration counter
            mIter++;

            // return vector of design criteria
            return tMat;
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
                    Matrix< DDRMat > tDIQIDAdv = tDCriteriaDAdv.get_row( i );

                    auto tItrMin = std::min_element( tDIQIDAdv.data(), tDIQIDAdv.data() + tDIQIDAdv.numel() );
                    auto tIndMin = std::distance( tDIQIDAdv.data(), tItrMin );

                    auto tItrMax = std::max_element( tDIQIDAdv.data(), tDIQIDAdv.data() + tDIQIDAdv.numel() );
                    auto tIndMax = std::distance( tDIQIDAdv.data(), tItrMax );

                    MORIS_LOG_INFO( "Criteria(%i): norm = %e   min = %e  (index = %i)   max = %e  (index = %i)",
                            i,
                            norm( tDIQIDAdv ),
                            tDIQIDAdv.min(),
                            tIndMin,
                            tDIQIDAdv.max(),
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
                const std::shared_ptr< mtk::Mesh_Manager > aMTKPerformer,
                xtk::Model* const &                        aXTKPerformer )
        {
            // Output T-matrices if requested
            std::string tTmatrixFileName = aXTKPerformer->get_T_matrix_output_file_name();
            if ( tTmatrixFileName != "" )
            {
                mPerformerManager->mMTKPerformer( 1 )->get_mesh_pair( 0 ).get_integration_mesh()->save_IG_node_TMatrices_to_file( tTmatrixFileName );

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

    }    // namespace wrk
} /* namespace moris */
