

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "fn_Exec_load_user_library.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_XTK_Model.hpp"
#include "typedefs.hpp"
#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_PRM_HMR_Parameters.hpp"
#include "cl_PRM_XTK_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

#include "cl_GEN_Geometry_Field_HMR.hpp"

using namespace moris;

Matrix< DDRMat > tSimpleFunc( moris::Cell< Matrix< DDRMat > >         & aCoeff )
{
    // f(x,y) = 2*x + 3*(y*y)

    uint tSize = aCoeff.size();

    Matrix< DDRMat > tAllVals(tSize,1);

    for(uint i=0; i<tSize; i++)
    {
        tAllVals(i) = 2*aCoeff(i)(0) + 3*(aCoeff(i)(1)*aCoeff(i)(1));
    }

    return tAllVals;
}

TEST_CASE("WRK L2 test","[WRK_L2_test]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Geometry";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         //----- HRM parameter list --------

         ParameterList tParameters = prm::create_hmr_parameter_list();
         tParameters.set( "number_of_elements_per_dimension", std::string("4,4"));
         tParameters.set( "domain_dimensions",                std::string("4,4"));
         tParameters.set( "domain_offset",                    std::string("-2.0,-2.0") );
         tParameters.set( "domain_sidesets",                  std::string("1,2,3,4"));
         tParameters.set( "lagrange_output_meshes",           std::string("0") );
         tParameters.set( "lagrange_orders",                  std::string("2"));
         tParameters.set( "lagrange_pattern",                 std::string("0") );
         tParameters.set( "bspline_orders",                   std::string("1"));
         tParameters.set( "bspline_pattern",                 std::string("0"));
         tParameters.set( "lagrange_to_bspline",              std::string("0,-1") );
         tParameters.set( "truncate_bsplines", 1 );
         tParameters.set( "refinement_buffer", 2 );
         tParameters.set( "staircase_buffer", 2 );
         tParameters.set( "initial_refinement", 0 );

         tParameters.set( "use_multigrid", 0 );
         tParameters.set( "severity_level", 2 );
         tParameters.set( "use_number_aura", 1);

         //----- GE parameter list --------

         moris::Cell<moris::Cell<ParameterList>> tParameterLists(2);
         tParameterLists(0).resize(1);
         tParameterLists(1).resize(1);
         tParameterLists(0)(0) = moris::prm::create_gen_parameter_list();
         tParameterLists(0)(0).set("initial_advs", "0.017, 0.9753, 0.0");
         tParameterLists(0)(0).set("HMR_refinements", 2);
         tParameterLists(1)(0) = moris::prm::create_geometry_parameter_list();
         tParameterLists(1)(0).set("type", "circle");
         tParameterLists(1)(0).set("geometry_variable_indices", "0, 2");
         tParameterLists(1)(0).set("adv_indices", "0, 1");
         tParameterLists(1)(0).set("constant_parameters", "1.0");

         //----- XTK parameter list --------

         ParameterList tParametersXTK = prm::create_xtk_parameter_list();
         tParametersXTK.set( "decompose",                 true );
         tParametersXTK.set( "decomposition_type",        std::string("conformal") );
         tParametersXTK.set( "enrich",                    true );
         tParametersXTK.set( "basis_rank",                std::string("bspline") );
         tParametersXTK.set( "enrich_mesh_indices",       std::string("0,1") );
         tParametersXTK.set( "ghost_stab",                false );
         tParametersXTK.set( "multigrid",                 false );
         tParametersXTK.set( "verbose",                   true );
         tParametersXTK.set( "print_enriched_ig_mesh",    true );
         tParametersXTK.set( "exodus_output_XTK_ig_mesh", true );
         tParametersXTK.set( "high_to_low_dbl_side_sets", false );

         //---------------------------------------------------------------------------------------
         //                               Stage 1: HMR refinement
         //---------------------------------------------------------------------------------------

         std::shared_ptr< hmr::HMR > tHMRPerformer = std::make_shared< hmr::HMR >( tParameters );

         std::shared_ptr< ge::Geometry_Engine > tGENPerformer = std::make_shared< ge::Geometry_Engine >( tParameterLists );

         tGENPerformer->register_mesh( tHMRPerformer );

         // uniform initial refinement
         tHMRPerformer->perform_initial_refinement( 0 );

         // HMR refined by GE
         tGENPerformer->perform();

         std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer_HMR = std::make_shared< mtk::Mesh_Manager >();
         tHMRPerformer->set_performer( tMTKPerformer_HMR );

         // HMR finalize
         tHMRPerformer->perform();

         //---------------------------------------------------------------------------------------
         //                               Stage 2: XTK
         //---------------------------------------------------------------------------------------

         // Register Mesh to Ge
         tGENPerformer->register_mesh( tMTKPerformer_HMR.get() );

         // create MTK performer - will be used for XTK mesh
         std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer_XTK = std::make_shared< mtk::Mesh_Manager >();

         // create XTK performer
         std::shared_ptr< xtk::Model > tXTKPerformer = std::make_shared< xtk::Model >( tParametersXTK );
         tXTKPerformer->set_geometry_engine( tGENPerformer.get() );
         tXTKPerformer->set_input_performer( tMTKPerformer_HMR );
         tXTKPerformer->set_output_performer( tMTKPerformer_XTK );

         tXTKPerformer->perform();

         //---------------------------------------------------------------------------------------
         //                               Stage 3: build field
         //---------------------------------------------------------------------------------------

//         std::shared_ptr< ge::GEN_Property > tFieldProperty = std::make_shared< ge::GEN_Property >();
//         tFieldProperty->set_val_function( tSimpleFunc );
//
//         ge::GEN_Field tField( tFieldProperty );
//         tField.initialize( tMTKPerformer_XTK.get() );

         //---------------------------------------------------------------------------------------
         //                               Stage 3: MDL
         //---------------------------------------------------------------------------------------

//         std::shared_ptr< mdl::Model > tMDLPerformer = std::make_shared< mdl::Model >( nullptr, 0 );
    }
}
