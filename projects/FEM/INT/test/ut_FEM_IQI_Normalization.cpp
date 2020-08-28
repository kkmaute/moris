#include "catch.hpp"

#include "cl_FEM_Model.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_PRM_HMR_Parameters.hpp"

namespace moris
{
    namespace fem
    {
        TEST_CASE("IQI Scaling", "[fem], [IQI]")
        {
            // Create meshes
            mtk::Mesh_Manager tMeshManager;

            // FEM parameter lists
            Cell<Cell<ParameterList>> tParameterList;
            tParameterList.resize( 6 );

            // create parameter list for constitutive model 1
//            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
//            tParameterList( 1 )( 0 ).set( "constitutive_name", "CMStrucLinIso1");
//            tParameterList( 1 )( 0 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
//            tParameterList( 1 )( 0 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
//            tParameterList( 1 )( 0 ).set( "properties",        "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio");
//
//            // create parameter list for IWG 1
//            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
//            tParameterList( 3 )( 0 ).set( "IWG_name",                   "IWGBulkU_1");
//            tParameterList( 3 )( 0 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
//            tParameterList( 3 )( 0 ).set( "dof_residual",               "UX,UY");
//            tParameterList( 3 )( 0 ).set( "master_dof_dependencies",    "UX,UY");
//            tParameterList( 3 )( 0 ).set( "master_constitutive_models", "CMStrucLinIso1,ElastLinIso");
//            tParameterList( 3 )( 0 ).set( "mesh_set_names",             "HMR_dummy");

            // create parameter list for IQI 4
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( 0 ).set( "IQI_name",                   "Volume");
            tParameterList( 4 )( 0 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::VOLUME ) );
            //tParameterList( 4 )( 0 ).set( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UX ) );
            //tParameterList( 4 )( 0 ).set( "master_dof_dependencies",    "UX,UY");
            //tParameterList( 4 )( 0 ).set( "vectorial_field_index",      0 );
            tParameterList( 4 )( 0 ).set( "mesh_set_names",             "HMR_dummy");

            // create computation  parameter list
            tParameterList( 5 ).resize( 1 );
            tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();

            // HMR
            ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "2, 2");
            tParameters.set( "domain_dimensions", "2, 2");
            tParameters.set( "domain_offset", "-1.0, -1.0");
            tParameters.set( "domain_sidesets", "1,2,3,4");
            tParameters.set( "lagrange_output_meshes", "0");

            tParameters.set( "lagrange_orders", "1");
            tParameters.set( "lagrange_pattern", "0");
            tParameters.set( "bspline_orders", "2");
            tParameters.set( "bspline_pattern", "0");

            tParameters.set( "lagrange_to_bspline", "0");

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );

            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement( 0 );
            tHMR.finalize();

            // Register mesh pair
            mtk::Interpolation_Mesh* tInterpolationMesh = tHMR.create_interpolation_mesh(0);
            tMeshManager.register_mesh_pair(tInterpolationMesh,
                    create_integration_mesh_from_interpolation_mesh(MeshType::HMR, tInterpolationMesh));

            // Create FEM model
            FEM_Model tModel(&tMeshManager, 0, tParameterList, nullptr);

            // Set requested IQI name
            tModel.set_requested_IQI_names({"Volume"});

            // Test IQI
            tModel.compute_IQIs();
            moris::print(tModel.get_IQI_values(), "get IQI values");
        }
    }
}
