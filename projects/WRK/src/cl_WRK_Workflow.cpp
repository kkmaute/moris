#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "fn_WRK_perform_refinement.hpp"

#include "cl_HMR.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_MDL_Model.hpp"

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

            // uniform initial refinement
            mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement( 0 );

            // HMR refined by GE
            perform_refinement(mPerformerManager->mHMRPerformer( 0 ), {mPerformerManager->mGENPerformer( 0 )});

            // HMR finalize
            mPerformerManager->mHMRPerformer( 0 )->perform();

            // Stage 2: Initialize Level set field in GEN -----------------------------------------------
            mPerformerManager->mGENPerformer( 0 )->compute_level_set_data(
                    mPerformerManager->mMTKPerformer( 0 )->get_interpolation_mesh(0));

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

            // XTK perform - decompose - enrich - ghost - multigrid
            mPerformerManager->mXTKPerformer( 0 )->perform();

            // Output GEN fields, if requested
            mPerformerManager->mGENPerformer( 0 )->output_fields(
                    mPerformerManager->mMTKPerformer( 1 )->get_integration_mesh( 0 ));

            // Assign PDVs
            mPerformerManager->mGENPerformer( 0 )->create_pdvs( mPerformerManager->mMTKPerformer( 1 ) );

            // Stage 3: MDL perform ---------------------------------------------------------------------

            mPerformerManager->mMDLPerformer( 0 )->initialize();

            mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface(
                    mPerformerManager->mGENPerformer( 0 )->get_design_variable_interface() );

            mPerformerManager->mGENPerformer( 0 )->communicate_requested_IQIs();

            // Build MDL components and solve
            mPerformerManager->mMDLPerformer( 0 )->perform();

            moris::Cell< moris::Matrix< DDRMat > > tVal = mPerformerManager->mMDLPerformer( 0 )->get_IQI_values();

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

            return mPerformerManager->mGENPerformer( 0 )->get_dcriteria_dadv();
        }

        //--------------------------------------------------------------------------------------------------------------
    } /* namespace mdl */
} /* namespace moris */
