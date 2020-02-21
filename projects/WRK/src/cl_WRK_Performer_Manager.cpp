
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

#include "cl_HMR.hpp"

#include "cl_MDL_Model.hpp"

#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_WRK_Performer_Manager.hpp"


#include "fn_Exec_load_user_library.hpp"

namespace moris
{
namespace wrk
{
//------------------------------------------------------------------------------
Performer_Manager::Performer_Manager( std::shared_ptr< Library_IO > aLibrary ) :mLibrary(aLibrary)
{

}

void Performer_Manager::initialize()
{
    mHMRPerformer.resize( 1 );
    mGENPerformer.resize( 1 );
    mXTKPerformer.resize( 1 );
    mMTKPerformer.resize( 2 );
    mMDLPerformer.resize( 1 );

    // load the HMR parameter list
    std::string tHMRString = "HMRParameterList";
    MORIS_PARAMETER_FUNCTION tHMRParameterListFunc = mLibrary->load_parameter_file( tHMRString );
    moris::Cell< moris::Cell< ParameterList > > tHMRParameterList;
    tHMRParameterListFunc( tHMRParameterList );

    std::string tGENString = "GENParameterList";
    MORIS_PARAMETER_FUNCTION tGENParameterListFunc = mLibrary->load_parameter_file( tGENString );
    moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
    tGENParameterListFunc( tGENParameterList );

    mHMRPerformer( 0 ) = std::make_shared< hmr::HMR >( tHMRParameterList( 0 )( 0 ) );

    // Create GE with parameter list
    mGENPerformer( 0 ) = std::make_shared< ge::GEN_Geometry_Engine >( tGENParameterList(0)(0) );
}

void Performer_Manager::perform_refinement()
{
    std::string tGENString = "GENParameterList";
    MORIS_PARAMETER_FUNCTION tGENParameterListFunc = mLibrary->load_parameter_file( tGENString );
    moris::Cell< moris::Cell< ParameterList > > tGENParameterList;
    tGENParameterListFunc( tGENParameterList );

    std::string tHMRString = "HMRParameterList";
    MORIS_PARAMETER_FUNCTION tHMRParameterListFunc = mLibrary->load_parameter_file( tHMRString );
    moris::Cell< moris::Cell< ParameterList > > tHMRParameterList;
    tHMRParameterListFunc( tHMRParameterList );

    // initial refinement
    mHMRPerformer( 0 )->perform_initial_refinement( 0 );

    uint tLagrangeMeshIndex = 0;     // Fixme get from parameter list - output mesh

    std::shared_ptr< moris::hmr::Mesh > tMesh = mHMRPerformer( 0 )->create_mesh( tLagrangeMeshIndex );

    for( sint k=0; k< tHMRParameterList( 0 )( 0 ).get< moris::sint >( "adaptive_refinement_level"); ++k )
    {
        ge::GEN_Geometry_Engine tGE( tGENParameterList(0)(0) );
        tGE.register_mesh( tMesh );

        tGE.initialize( mLibrary );

        moris::Cell< moris::Matrix< DDRMat > > tValues;

        tGE.get_field_values_for_all_geometries( tValues );


        for( uint Ik=0; Ik < tValues.size(); ++Ik )
        {
            mHMRPerformer( 0 )->based_on_field_put_elements_on_queue( tValues( Ik ), tLagrangeMeshIndex );
        }

        mHMRPerformer( 0 )->perform_refinement_based_on_working_pattern( 0, false );
    }
}

void Performer_Manager::perform()
{
    this->perform_refinement();

    uint tLagrangeMeshIndex = 0;     // Fixme get from parameter list - output mesh

    mHMRPerformer( 0 )->finalize();

    mHMRPerformer( 0 )->save_to_exodus( 0, "./hmr_exo/Box.e" );

    std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = mHMRPerformer( 0 )->create_interpolation_mesh( tLagrangeMeshIndex );
    std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = mHMRPerformer( 0 )->create_integration_mesh( 1, 0, *tInterpolationMesh );

    mMTKPerformer( 0 ) =std::make_shared< mtk::Mesh_Manager >();
    mMTKPerformer( 0 )->register_mesh_pair( tInterpolationMesh.get(), tIntegrationMesh.get() );
    size_t tSpatialDimension = tInterpolationMesh->get_spatial_dim();

    // Register Mesh to Ge
    mGENPerformer( 0 )->register_mesh( mMTKPerformer( 0 ).get() );

    // Initialize call - build Geometries
    mGENPerformer( 0 )->initialize( mLibrary );

    mXTKPerformer( 0 ) = std::make_shared< xtk::Model >( tSpatialDimension, tInterpolationMesh.get(), *(mGENPerformer( 0 ).get()) );

    mXTKPerformer( 0 )->mVerbose = true;

    Cell<enum Subdivision_Method> tDecompositionMethods;

    if( tSpatialDimension == 2 )
    {
        //Specify decomposition Method and Cut Mesh ---------------------------------------
        tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
    }
    else if( tSpatialDimension == 3 )
    {
        tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    }

    mXTKPerformer( 0 )->decompose(tDecompositionMethods);

    mXTKPerformer( 0 )->perform_basis_enrichment( EntityRank::BSPLINE,0 );
    mXTKPerformer( 0 )->construct_face_oriented_ghost_penalization_cells();

    // get meshes
    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = mXTKPerformer( 0 )->get_enriched_interp_mesh();
    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = mXTKPerformer( 0 )->get_enriched_integ_mesh();

    // place the pair in mesh manager
    mMTKPerformer( 1 ) =std::make_shared< mtk::Mesh_Manager >();
    mMTKPerformer( 1 )->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

    // create model
    mMDLPerformer( 0 ) = std::make_shared< mdl::Model >( mLibrary,
                                                         mMTKPerformer( 1 ).get(),
                                                         0 );

    mMDLPerformer( 0 )->solve();
}



    } /* namespace mdl */
} /* namespace moris */
