
#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"

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

        Workflow::Workflow( wrk::Performer_Manager * aPerformerManager ) : mPerformerManager( aPerformerManager )
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        void Workflow::initialize(Matrix<DDRMat>& aADVs, Matrix<DDRMat>& aLowerBounds, Matrix<DDRMat>& aUpperBounds)
        {
            aADVs = mPerformerManager->mGENPerformer( 0 )->get_advs();
            aLowerBounds = mPerformerManager->mGENPerformer( 0 )->get_lower_bounds();
            aUpperBounds = mPerformerManager->mGENPerformer( 0 )->get_upper_bounds();

        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow::perform(Matrix<DDRMat> aNewADVs)
        {
            //---------------------------------------------------------------------------------------
            //                               Stage 1: HMR refinement
            //---------------------------------------------------------------------------------------

            // Set new advs in GE
            mPerformerManager->mGENPerformer( 0 )->set_advs(aNewADVs);

            // uniform initial refinement
            mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement( 0 );

            // HMR refined by GE
            mPerformerManager->mGENPerformer( 0 )->perform();

            // HMR finalize
            mPerformerManager->mHMRPerformer( 0 )->perform();

            //---------------------------------------------------------------------------------------
            //                               Stage 2: XTK
            //---------------------------------------------------------------------------------------

            // Register Mesh to Ge
            mPerformerManager->mGENPerformer( 0 )->register_mesh( mPerformerManager->mMTKPerformer( 0 ).get() );

            // XTK perform - decompose - enrich - ghost - multigrid
            mPerformerManager->mXTKPerformer( 0 )->perform();

            //---------------------------------------------------------------------------------------
            //                               Stage 3: MDL perform
            //---------------------------------------------------------------------------------------

            // Build MDL components and solve
            mPerformerManager->mMDLPerformer( 0 )->perform();

            return Matrix<DDRMat>(1, 1 , 0.0);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow::get_criteria(Matrix<DDRMat> aNewADVs)
        {
            mPerformerManager->mMDLPerformer( 0 )->set_design_variable_interface( mPerformerManager->mGENPerformer( 0 )
                                                            ->get_design_variable_interface() );

            // FIXME set requeted IQIs
//            mPerformerManager->mMDLPerformer( 0 )->

            mPerformerManager->mMDLPerformer( 0 )->perform_sensitivity_analysis();

            return Matrix<DDRMat>(1, 1 , 0.0);
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix<DDRMat> Workflow::get_dcriteria_dadv()
        {
            return Matrix<DDRMat>(1, 1 , 0.0);
        }

        //--------------------------------------------------------------------------------------------------------------

    } /* namespace mdl */
} /* namespace moris */
