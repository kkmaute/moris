

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

#include "cl_Library_IO.hpp"
#include "cl_Communication_Tools.hpp"

#include "typedefs.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Field_Analytic.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "st_MTK_Mesh_Pair.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Database.hpp" //HMR/src

#include "cl_WRK_perform_refinement.hpp"
#include "cl_WRK_perform_remeshing.hpp"

#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_MORIS_GENERAL_Parameters.hpp"


using namespace moris;

moris::real tCircle_1( const Matrix< DDRMat > & aCoords )
{
    return norm( aCoords ) - 1.5237 ;
}

moris::real tCircle_2( const Matrix< DDRMat > & aCoords )
{
    return norm( aCoords ) - 0.8713 ;
}

TEST_CASE("WRK L2 test","[WRK_L2_test]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Geometry";

         //----- HRM parameter list --------

         ParameterList tParameters = prm::create_hmr_parameter_list();

         tParameters.set( "number_of_elements_per_dimension", "4,   4" );
         tParameters.set( "domain_dimensions",                "4.0,   4.0" );
         tParameters.set( "domain_offset",                    "-2.0,  -2.0"  );
         tParameters.set( "domain_sidesets",                  "1,2,3,4");
         tParameters.set( "lagrange_output_meshes",           "0");

         tParameters.set( "lagrange_orders",  "1" );
         tParameters.set( "lagrange_pattern",  "0" )  ;
         tParameters.set( "bspline_orders",   "1" );
         tParameters.set( "bspline_pattern",   "0" )  ;

         tParameters.set( "lagrange_to_bspline", "0") ;

         tParameters.set( "truncate_bsplines",  1 );
         tParameters.set( "refinement_buffer",  0 );
         tParameters.set( "staircase_buffer",   1 );
         tParameters.set( "initial_refinement", "1" );
         tParameters.set( "initial_refinement_pattern", "0" );

         tParameters.set( "use_number_aura", 1);

         tParameters.set( "use_multigrid",  0 );
         tParameters.set( "severity_level", 0 );

         ParameterList tRefinementParameters;
         prm::create_refinement_parameterlist( tRefinementParameters );
         tRefinementParameters.set( "field_names" , "Circle_1" );
         tRefinementParameters.set( "levels_of_refinement" , "1" );
         tRefinementParameters.set( "refinement_pattern" , "0" );

         ParameterList tRemeshingParameters;
         prm::create_remeshing_parameterlist( tRemeshingParameters );
         tRemeshingParameters.set( "mode" , "ab_initio" );
         tRemeshingParameters.set( "remeshing_levels_of_refinement" , "1" );
         tRemeshingParameters.set( "remeshing_refinement_pattern" , "0" );

         //---------------------------------------------------------------------------------------
         //                               Stage 1: HMR refinement
         //---------------------------------------------------------------------------------------

         std::shared_ptr< hmr::HMR > tHMRPerformer = std::make_shared< hmr::HMR >( tParameters );

         // uniform initial refinement
         tHMRPerformer->perform_initial_refinement();

         std::shared_ptr< hmr::Database > tHMRDatabase = tHMRPerformer->get_database();

         hmr::Interpolation_Mesh_HMR * tInterpolationMesh = new hmr::Interpolation_Mesh_HMR(
                 tHMRDatabase,
                 1,
                 0,
                 1,
                 0);

         mtk::Mesh_Pair tMeshPair;
         tMeshPair.mInterpolationMesh = tInterpolationMesh;
         tMeshPair.mIsOwned   = true;

         moris::Cell< mtk::Field * > tFields( 1, nullptr );
         tFields( 0 ) = new mtk::Field_Analytic( &tMeshPair );
         tFields( 0 )->set_label( "Circle_1" );
         reinterpret_cast< mtk::Field_Analytic* >(tFields( 0 ))->evaluate_scalar_function( tCircle_1 );

         //tFields( 0 )->save_field_to_exodus( "Remeshing_Field1.exo");

         wrk::Refinement_Mini_Performer tRefinementPerformer( tRefinementParameters );
         tRefinementPerformer.perform_refinement( tFields, tHMRPerformer );

         //delete( tInterpolationMesh );

         std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer_HMR = std::make_shared< mtk::Mesh_Manager >();
         tHMRPerformer->set_performer( tMTKPerformer_HMR );

         hmr::Interpolation_Mesh_HMR * tInterpolationMeshNew = new hmr::Interpolation_Mesh_HMR(
                 tHMRDatabase,
                 1,
                 0,
                 1,
                 0);

         delete tMeshPair.mInterpolationMesh;
         tMeshPair.mInterpolationMesh = tInterpolationMeshNew;
         reinterpret_cast< mtk::Field_Analytic* >(tFields( 0 ))->evaluate_scalar_function( tCircle_1 );
         //tFields( 0 )->save_field_to_exodus( "Remeshing_Field2.exo");
         reinterpret_cast< mtk::Field_Analytic* >(tFields( 0 ))->evaluate_scalar_function( tCircle_2 );
         //tFields( 0 )->save_field_to_exodus( "Remeshing_Field3.exo");


         //print( tFields( 0 )->get_nodal_values(), "Val1");

         tHMRDatabase->get_background_mesh()->update_database();
         tHMRDatabase->update_bspline_meshes();
         tHMRDatabase->update_lagrange_meshes();

         // HMR finalize
         tHMRPerformer->perform();

         Cell< std::shared_ptr< hmr::HMR > > tHMRPerformers( 1, tHMRPerformer );

         wrk::Remeshing_Mini_Performer tRemeshingPerformer( tRemeshingParameters );
         tRemeshingPerformer.perform_remeshing( tFields( 0 ), tHMRPerformers );

         std::shared_ptr< hmr::Database > tHMRDatabaseNew = tHMRPerformers( 0 )->get_database();

         hmr::Interpolation_Mesh_HMR * tInterpolationMeshNewMesh = new hmr::Interpolation_Mesh_HMR(
                 tHMRDatabaseNew,
                 1,
                 0,
                 1,
                 0);

         delete tMeshPair.mInterpolationMesh;
         tMeshPair.mInterpolationMesh = tInterpolationMeshNewMesh;
         reinterpret_cast< mtk::Field_Analytic* >(tFields( 0 ))->evaluate_scalar_function( tCircle_2 );
         //tFields( 0 )->save_field_to_exodus( "Remeshing_Field4.exo");


         delete tFields( 0 );



    }
}