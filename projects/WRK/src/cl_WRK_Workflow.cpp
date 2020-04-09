

#include "cl_Stopwatch.hpp" //CHR/src

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_HMR.hpp"

#include "cl_MDL_Model.hpp"

#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"


#include "fn_Exec_load_user_library.hpp"

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

    // uniform initial refinrment
    mPerformerManager->mHMRPerformer( 0 )->perform_initial_refinement( 0 );

    // perform refinement
    mPerformerManager->mGENPerformer( 0 )->perform();

    // perform finalize HMR
    mPerformerManager->mHMRPerformer( 0 )->perform();

    //-----------------------------------------------------------------------------z----------
    //                               Stage 2: XTK
    //---------------------------------------------------------------------------------------

    // Register Mesh to Ge
    mPerformerManager->mGENPerformer( 0 )->register_mesh( mPerformerManager->mMTKPerformer( 0 ).get() );

//    mPerformerManager->mXTKPerformer( 0 )->mVerbose = true;       //FIXME replace with logger

    mPerformerManager->mXTKPerformer( 0 )->perform();

    // get meshes
    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = mPerformerManager->mXTKPerformer( 0 )->get_enriched_interp_mesh();
    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = mPerformerManager->mXTKPerformer( 0 )->get_enriched_integ_mesh();

    // place the pair in mesh manager
    mPerformerManager->mMTKPerformer( 1 )->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

    tEnrIntegMesh.print();

    if(true)
    {
        tEnrIntegMesh.deactivate_empty_sets();
        // Write mesh
        moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
        writer.write_mesh("", "./xtk_exo/xtk_temp.exo");

        // Write the fields
        writer.set_time(0.0);
        writer.close_file();

//            moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh();
//            tIntegMesh1->create_output_mesh(tEnrIgMeshFileName);
//            delete tIntegMesh1;
    }



    //---------------------------------------------------------------------------------------
    //                               Stage 3: MDL perform
    //---------------------------------------------------------------------------------------
    // create model
    mPerformerManager->mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mPerformerManager->mLibrary,
                                                                            mPerformerManager->mMTKPerformer( 1 ).get(),
                                                                            0 );

    mPerformerManager->mMDLPerformer( 0 )->solve();
}



    } /* namespace mdl */
} /* namespace moris */
