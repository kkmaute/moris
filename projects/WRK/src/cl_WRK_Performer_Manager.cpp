
#include "cl_Stopwatch.hpp" //CHR/src

// fixme: temporary
#include "cl_Map.hpp"
#include "fn_unique.hpp"
#include "fn_sum.hpp" // for check
#include "fn_iscol.hpp"
#include "fn_trans.hpp"
#include "op_equal_equal.hpp"

#include "MTK_Tools.hpp"
#include "cl_MTK_Enums.hpp"              //MTK/src
#include "cl_MTK_Mesh_Manager.hpp"       //MTK/src

#include "cl_WRK_Performer_Manager.hpp"


#include "fn_Exec_load_user_library.hpp"

namespace moris
{
    namespace wrk
    {
//------------------------------------------------------------------------------
    Performer_Manager::Performer_Manager( const std::string & aInputFilePath ) : mInputFile( aInputFilePath )
    {
        mLibrary = std::make_shared< Library_IO >( mInputFile );

        // load the HMR parameter list
        std::string tHMRString = "HMRParameterList";
        MORIS_PARAMETER_FUNCTION tHMRParameterListFunc = mLibrary->load_parameter_file( tHMRString );
        moris::Cell< moris::Cell< ParameterList > > tHMRParameterList;
        tHMRParameterListFunc( tHMRParameterList );

        std::string tGENString = "GENParameterList";
        MORIS_PARAMETER_FUNCTION tGENParameterListFunc = mLibrary->load_parameter_file( tGENString );
        moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
        tGENParameterListFunc( tGENParameterList );

        uint tLagrangeMeshIndex = 0;

        hmr::HMR tHMR( tParameters( 0 )( 0 ) );

        // initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        for( uint k=0; k<2; ++k )
        {
            GEN_Geometry_Engine tGE( tGENParameterList(0)(0) );
            tGE.register_mesh( tMesh );

            tGE.initialize( mLibrary );

            moris::Cell< moris::Matrix< DDRMat > > tValues;

            tGE.get_field_values_for_all_geometries( tValues );


            for( uint Ik=0; Ik < tValues.size(); ++Ik )
            {
                tHMR.based_on_field_put_elements_on_queue( tValues( Ik ), tLagrangeMeshIndex );
            }

            tHMR.perform_refinement_based_on_working_pattern( 0, false );
        }

        tHMR.finalize();

//        tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = tHMR.create_integration_mesh( 1, 0, *tInterpolationMesh );

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair( tInterpolationMesh.get(), tIntegrationMesh.get() );
        size_t tSpatialDimension = tInterpolationMesh->get_spatial_dim();

        //-------------------------------------------------------------------------------

        // Create GE with parameter list
        GEN_Geometry_Engine tGE( tGENParameterList(0)(0) );

        // Register Mesh to Ge
        tGE.register_mesh( tMeshManager );

        // Initialize call - build Geometries
        tGE.initialize( mLibrary );

        // Create phase table for geometry engine
        moris::ge::GEN_Phase_Table tPhaseTable ( 1, Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tSpatialDimension);

        xtk::Model tXTKModel(tSpatialDimension, tInterpolationMesh.get(), tGeometryEngine);

        tXTKModel.mVerbose = false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);

        // get meshes
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0 );

        tModel->solve();
    }



    } /* namespace mdl */
} /* namespace moris */
