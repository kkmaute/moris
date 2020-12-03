#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "fn_WRK_perform_refinement.hpp"

#include "cl_HMR.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MDL_Model.hpp"

#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "cl_Stopwatch.hpp"

#include "fn_norm.hpp"

namespace moris
{
    namespace wrk
    {
        //--------------------------------------------------------------------------------------------------------------

        Workflow::Workflow( wrk::Performer_Manager * aPerformerManager )
        : mPerformerManager( aPerformerManager )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Workflow::initialize(
                Matrix<DDRMat>& aADVs,
                Matrix<DDRMat>& aLowerBounds,
                Matrix<DDRMat>& aUpperBounds)
        {
            // Stage 1: HMR refinement -------------------------------------------------------------------

            // start timer
            tic tTimer;

            // uniform initial refinement
            mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement();

            // HMR refined by GE
            perform_refinement(mPerformerManager->mHMRPerformer( 0 ), {mPerformerManager->mGENPerformer( 0 )});

            // HMR finalize
            mPerformerManager->mHMRPerformer( 0 )->perform();

            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
            moris::real tElapsedTimeMax = max_all( tElapsedTime );

            if ( par_rank() == 0 )
            {
                MORIS_LOG_INFO( "HMR: Total time for refinement and mesh creation is %5.3f seconds.",
                        ( double ) tElapsedTimeMax / 1000);
            }

            // Stage 2: Initialize Level set field in GEN -----------------------------------------------
            mPerformerManager->mGENPerformer( 0 )->compute_level_set_data(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh(0) );

            // Get ADVs
            aADVs        = mPerformerManager->mGENPerformer( 0 )->get_advs();
            aLowerBounds = mPerformerManager->mGENPerformer( 0 )->get_lower_bounds();
            aUpperBounds = mPerformerManager->mGENPerformer( 0 )->get_upper_bounds();
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow::perform(Matrix<DDRMat> aNewADVs)
        {
            // Set new advs in GE
            mPerformerManager->mGENPerformer( 0 )->set_advs(aNewADVs);

            // Stage 1: HMR refinement

            // Stage 2: XTK -----------------------------------------------------------------------------
            mPerformerManager->create_xtk();

            // Compute level set data in GEN
            mPerformerManager->mGENPerformer( 0 )->compute_level_set_data(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ));

            // Output GEN fields, if requested
            mPerformerManager->mGENPerformer( 0 )->output_fields(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh( 0 ));

            // XTK perform - decompose - enrich - ghost - multigrid
            mPerformerManager->mXTKPerformer( 0 )->perform();

            // Assign PDVs
            mPerformerManager->mGENPerformer( 0 )->create_pdvs( mPerformerManager->mMTKPerformer( 1 ) );

            // Stage 3: MDL perform ---------------------------------------------------------------------

            mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface(
                                mPerformerManager->mGENPerformer( 0 )->get_design_variable_interface() );

            mPerformerManager->mMDLPerformer( 0 )->initialize();

            mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

            // Build MDL components and solve
            mPerformerManager->mMDLPerformer( 0 )->perform();

            moris::Cell< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

            // Communicate IQIs
            for( uint iIQIIndex = 0; iIQIIndex < tVal.size(); iIQIIndex++ )
            {
                tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
            }

            moris::Matrix< DDRMat > tMat( tVal.size(), 1, 0.0 );

            for( uint Ik = 0; Ik < tVal.size(); Ik ++ )
            {
                tMat( Ik ) = tVal( Ik )( 0 );
            }

            return tMat;
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow::compute_dcriteria_dadv()
        {
            mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

            mPerformerManager->mMDLPerformer( 0 )->perform( 1 );

            Matrix<DDRMat> tDCriteriaDAdv = mPerformerManager->mGENPerformer( 0 )->get_dcriteria_dadv();

            if (par_rank() == 0)
            {
                MORIS_LOG_INFO ( "--------------------------------------------------------------------------------");
                MORIS_LOG_INFO ( "Gradients of design criteria wrt ADVs:");

                for (uint i=0;i<tDCriteriaDAdv.n_rows();++i)
                {
                    Matrix<DDRMat> tDIQIDAdv = tDCriteriaDAdv.get_row(i);

                    auto tItrMin = std::min_element(tDIQIDAdv.data(),tDIQIDAdv.data()+tDIQIDAdv.numel());
                    auto tIndMin = std::distance(tDIQIDAdv.data(),tItrMin);

                    auto tItrMax = std::max_element(tDIQIDAdv.data(),tDIQIDAdv.data()+tDIQIDAdv.numel());
                    auto tIndMax = std::distance(tDIQIDAdv.data(),tItrMax);

                    MORIS_LOG_INFO ( "Criteria(%i): norm = %e   min = %e  (index = %i)   max = %e  (index = %i)",
                                     i, norm(tDIQIDAdv),tDIQIDAdv.min(),tIndMin,tDIQIDAdv.max(),tIndMax);
                }

                MORIS_LOG_INFO ( "--------------------------------------------------------------------------------");
            }

            return tDCriteriaDAdv;
        }

        //--------------------------------------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */
