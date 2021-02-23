

#include "catch.hpp"

#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "op_times.hpp"
#include "op_equal_equal.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow.hpp"
#include "cl_Library_IO.hpp"
#include "cl_Communication_Tools.hpp"

#include "typedefs.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "fn_PRM_HMR_Parameters.hpp"

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

         //----- HRM parameter list --------

         ParameterList tParameters = prm::create_hmr_parameter_list();
         tParameters.set( "number_of_elements_per_dimension", "4,4");
         tParameters.set( "domain_dimensions",                "4,4");
         tParameters.set( "domain_offset",                    "-2.0,-2.0" );
         tParameters.set( "domain_sidesets",                  "1,2,3,4");
         tParameters.set( "lagrange_output_meshes",           "0" );
         tParameters.set( "lagrange_orders",                  "2");
         tParameters.set( "lagrange_pattern",                 "0" );
         tParameters.set( "bspline_orders",                   "1");
         tParameters.set( "bspline_pattern",                 "0");
         tParameters.set( "lagrange_to_bspline",              "0,-1" );
         tParameters.set( "truncate_bsplines", 1 );
         tParameters.set( "refinement_buffer", 2 );
         tParameters.set( "staircase_buffer", 2 );
         tParameters.set( "initial_refinement", "0" );
         tParameters.set( "initial_refinement_pattern", "0" );

         tParameters.set( "use_multigrid", 0 );
         tParameters.set( "severity_level", 2 );
         tParameters.set( "use_number_aura", 1);

         //---------------------------------------------------------------------------------------
         //                               Stage 1: HMR refinement
         //---------------------------------------------------------------------------------------

         std::shared_ptr< hmr::HMR > tHMRPerformer = std::make_shared< hmr::HMR >( tParameters );

         // uniform initial refinement
         tHMRPerformer->perform_initial_refinement();

         std::shared_ptr< mtk::Mesh_Manager > tMTKPerformer_HMR = std::make_shared< mtk::Mesh_Manager >();
         tHMRPerformer->set_performer( tMTKPerformer_HMR );

         // HMR finalize
         tHMRPerformer->perform();

    }
}
