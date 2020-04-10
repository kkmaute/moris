
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
//------------------------------------------------------------------------------
Workflow::Workflow( wrk::Performer_Manager * aPerformerManager ) : mPerformerManager( aPerformerManager )
{

}

void Workflow::initialize()
{
}

void Workflow::perform()
{
    //---------------------------------------------------------------------------------------
    //                               Stage 1: HMR refinement
    //---------------------------------------------------------------------------------------

    // uniform initial refinement
    mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement( 0 );

    // perform refinement
    mPerformerManager->mGENPerformer( 0 )->perform();

    // perform finalize HMR
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
}



    } /* namespace mdl */
} /* namespace moris */
