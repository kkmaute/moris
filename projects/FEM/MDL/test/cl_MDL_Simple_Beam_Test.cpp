
#include "catch.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src
#include "fn_inv.hpp" // to allow for matrix inverse call

#include <stdlib.h> // to include random number generator

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"               //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"               //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"            //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_PRM_HMR_Parameters.hpp" // paramter list for hmr mesh

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

namespace moris
{
namespace mdl
{


void tConstantValueFunction( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                           moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                           moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tSelectValueFunction( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                           moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                           moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                      0.0,                     0.0 },
                  {                    0.0,    aParameters( 0 )( 1 ),                     0.0 },
                  {                    0.0,                      0.0,    aParameters( 0 )( 2 )}};
}




////--------------------------------------------------------------------------------------------------------------------//
//TEST_CASE( "LinElasticBeam", "[moris],[mdl],[LinElasticBeam_s]" )
//{
//    uint tSizeIndicator = 1;
//    function_TEST_LinElastic_Beam(tSizeIndicator);
//}
////--------------------------------------------------------------------------------------------------------------------//
//TEST_CASE( "LinElasticBeam", "[moris],[mdl],[LinElasticBeam_m]" )
//{
//    uint tSizeIndicator = 2;
//    function_TEST_LinElastic_Beam(tSizeIndicator);
//}
////--------------------------------------------------------------------------------------------------------------------//
//TEST_CASE( "LinElasticBeam", "[moris],[mdl],[LinElasticBeam_l]" )
//{
//    uint tSizeIndicator = 3;
//    function_TEST_LinElastic_Beam(tSizeIndicator);
//}


void function_TEST_LinElastic_Beam(uint aSizeIndicator)
{
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%% CREATE MESH FOR PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    // Pretiming
    std::clock_t tTimeStampBeginSetup = std::clock();

    // create a list to write mesh parameters into (standard constructor used, list starts with default settings)
    moris::ParameterList tParameters = prm::create_hmr_parameter_list();

    // set parameters to create tensor-grid linear Lagrange mesh
    uint tLagrangeMeshIndex = 0;

    switch (aSizeIndicator) {

    case 1:
        tParameters.set( "number_of_elements_per_dimension", std::string("200,  5,  5") );
        break;

    case 2:
        tParameters.set( "number_of_elements_per_dimension", std::string("400, 10, 10") );
        break;

    case 3:
        tParameters.set( "number_of_elements_per_dimension", std::string("600, 15, 15") );
        break;

    default:
        MORIS_ASSERT( false, "Size of problem not indicated correctly.");
        break;
    }

//    tParameters.set( "number_of_elements_per_dimension", "400,  10,  10" );

    tParameters.set( "domain_dimensions", std::string("4.0, 0.1, 0.1") );
    tParameters.set( "domain_offset",     std::string("0.0, 0.0, 0.0") );

    tParameters.set( "domain_sidesets", std::string("1, 2, 3, 4, 5, 6") );

    tParameters.set( "truncate_bsplines",  1  );
    tParameters.set(   "lagrange_orders", std::string("1") );
    tParameters.set(  "lagrange_pattern", std::string("0") );
    tParameters.set(    "bspline_orders", std::string("1") );
    tParameters.set(   "bspline_pattern", std::string("0") );

    tParameters.set( "lagrange_output_meshes", std::string("0") );
    tParameters.set(  "lagrange_input_meshes", std::string("0") );

    tParameters.set( "lagrange_to_bspline", std::string("0") );

    tParameters.set( "use_multigrid", 0 );

    tParameters.set( "refinement_buffer", 1 );
    tParameters.set( "staircase_buffer", 1 );

    tParameters.insert( "initial_refinement", 0 );

    // initialize HMR mesh
    moris::hmr::HMR tHMR( tParameters );

    // NO outputs
    // get HMR Mesh Database for potential modifications to it
    auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >

    tDatabase->update_bspline_meshes();
    tDatabase->update_lagrange_meshes();

    // NO outputs
    // create initial Mesh from the information fed into the HMR database (parameters are part of the database)
    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

    // finalize HMR mesh to be used
    tHMR.finalize();

    // NO outputs
//    std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR >  tInterpolationMesh  = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//
//    // NO outputs
//    std::shared_ptr< moris::hmr::Integration_Mesh_HMR >    tIntegrationMesh    = tHMR.create_integration_mesh( 1, 0,*tInterpolationMesh );

    hmr::Interpolation_Mesh_HMR *      tInterpolationMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
    moris::hmr::Integration_Mesh_HMR * tIntegrationMesh = tHMR.create_integration_mesh( 1, 0,*tInterpolationMesh );

    // NO outputs
    // create Mesh Manager to be passed on to the Model/FEM
    moris::mtk::Mesh_Manager tMeshManager;
    // NO outputs
    tMeshManager.register_mesh_pair( tInterpolationMesh, tIntegrationMesh); // the mesh mangager ist not able to deal with HMR


    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    // --------------------------- PROPERTIES --------------------------- //

    // create Bulk Properties
    std::shared_ptr< fem::Property > tPropertyEModBulk = std::make_shared< fem::Property >();
    tPropertyEModBulk->set_parameters( { {{ 1e11 }} } );  // Young's Moudulus
    tPropertyEModBulk->set_val_function( tConstantValueFunction ); //sets the rule by which the Coefficient is computed

    std::shared_ptr< fem::Property > tPropertyNuBulk = std::make_shared< fem::Property >();
    tPropertyNuBulk->set_parameters( { {{ 0.3 }} } );        // Poisson's ratio
    tPropertyNuBulk->set_val_function( tConstantValueFunction );


    // create the Dirichlet Boundary Condition (Properties)
    std::shared_ptr< fem::Property > tPropertyDirichletValue = std::make_shared< fem::Property >();   // declare x,y,z displacements to be zero
    tPropertyDirichletValue->set_parameters( { {{ 0.0 }, {0.0}, {0.0}} } );
    tPropertyDirichletValue->set_val_function( tConstantValueFunction );

    std::shared_ptr< fem::Property > tPropertyDirichletSelect = std::make_shared< fem::Property >();  // enforce displacement in x, y, z
    tPropertyDirichletSelect->set_parameters( { {{ 1.0 }, {1.0}, {1.0}} } );
    tPropertyDirichletSelect->set_val_function( tSelectValueFunction );


    // create the Neumann Boundary Condition (Properties)
    std::shared_ptr< fem::Property > tPropertyNeumann = std::make_shared< fem::Property >();
    tPropertyNeumann->set_parameters( {{{ 0.0 } , { 0.0 }, {-4e5}}}  );
    tPropertyNeumann->set_val_function( tConstantValueFunction );



    // --------------------------- CONSTITUTIVE MODEL --------------------------- //

    // create Factory to construct constitutive models
    fem::CM_Factory tCMFactory;

    // define CONSTITUTIVE MODEL for the Bulk as 3D linear elastic
    std::shared_ptr< fem::Constitutive_Model > tConstitutiveModelMaterial = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tConstitutiveModelMaterial->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY , MSI::Dof_Type::UZ }} );
    tConstitutiveModelMaterial->set_property( tPropertyEModBulk, "YoungsModulus" );
    tConstitutiveModelMaterial->set_property( tPropertyNuBulk, "PoissonRatio" );
    tConstitutiveModelMaterial->set_space_dim( 3 );



    // ------------------------ STABILIZATION PARAMETERS ------------------------- //

    // create Factory to construct SPs
    fem::SP_Factory tSPFactory;

    // define SP needed for Nitsche formulation for the dirichlet boundary
    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheBC = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitscheBC->set_parameters( {{{ 100.0 }}} );
    tSPDirichletNitscheBC->set_property( tPropertyEModBulk, "Material", mtk::Master_Slave::MASTER );



    // ---------------------------------- IWGs ----------------------------------- //

    // create Factory to construct IWGs
    fem::IWG_Factory tIWGFactory;

    // Bulk material
    std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
    tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
    tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
    tIWGBulk->set_constitutive_model( tConstitutiveModelMaterial, "ElastLinIso", mtk::Master_Slave::MASTER );

    // Dirichlet Boundary
    std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
    tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
    tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
    tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitscheBC, "DirichletNitsche" );
    tIWGDirichlet->set_constitutive_model( tConstitutiveModelMaterial, "ElastLinIso", mtk::Master_Slave::MASTER );
    tIWGDirichlet->set_property( tPropertyDirichletValue, "Dirichlet", mtk::Master_Slave::MASTER );
    tIWGDirichlet->set_property( tPropertyDirichletSelect, "Select", mtk::Master_Slave::MASTER );

    // Neumann Boundary
    std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
    tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
    tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
    tIWGNeumann->set_property( tPropertyNeumann, "Neumann", mtk::Master_Slave::MASTER );


    // ------------------------- Create Set Info Table ---------------------------- //


    fem::Set_User_Info tFemSetBulk;
//    tFemSetBulk.set_mesh_index( tIntegrationMesh->get_set_index("HMR_dummy") );
    tFemSetBulk.set_mesh_index( 0 );
    tFemSetBulk.set_IWGs( { tIWGBulk } );

    fem::Set_User_Info tFemSetDirichlet;
    tFemSetDirichlet.set_mesh_index( tIntegrationMesh->get_side_set_index("SideSet_4") );
    tFemSetDirichlet.set_IWGs( { tIWGDirichlet } );

    fem::Set_User_Info tFemSetNeumann;
    tFemSetNeumann.set_mesh_index( tIntegrationMesh->get_side_set_index("SideSet_2") );
    tFemSetNeumann.set_IWGs( { tIWGNeumann } );

    // create a cell of set info
    moris::Cell< fem::Set_User_Info > tSetInfo( 3 );
    tSetInfo( 0 ) = tFemSetBulk;
    tSetInfo( 1 ) = tFemSetDirichlet;
    tSetInfo( 2 ) = tFemSetNeumann;

    // create model
    mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                           0,
                                           tSetInfo,
                                           0,
                                           false);


    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP SOLVER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

    moris::Cell< enum MSI::Dof_Type > tDofTypes( 3 );
    tDofTypes( 0 ) = MSI::Dof_Type::UX;
    tDofTypes( 1 ) = MSI::Dof_Type::UY;
    tDofTypes( 2 ) = MSI::Dof_Type::UZ;

    // Output set up time
    moris::real tTimeForSetUp = (moris::real) ( clock() - tTimeStampBeginSetup ) / CLOCKS_PER_SEC;
    MORIS_LOG_INFO( "Time to complete set up: %5.3f seconds.", ( double ) tTimeForSetUp );


     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 1: create linear solver and algorithm
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     dla::Solver_Factory  tSolverFactory;
     std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolverFactory.create_solver( sol::SolverType::AMESOS_IMPL );

     //if needed, specify parameters for algorithm

     dla::Linear_Solver tLinearSolver;
     tLinearSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 2: create nonlinear solver and algorithm
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     NLA::Nonlinear_Solver_Factory tNonlinearFactory;
     std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm =
             tNonlinearFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

     // specify parameters
     tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
     tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
     tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
     tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

     tNonlinearSolverAlgorithm->set_linear_solver( &tLinearSolver );

     NLA::Nonlinear_Solver tNonlinearSolver;
     tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 3: create time Solver and algorithm
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     tsa::Time_Solver_Factory tTimeSolverFactory;
     std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm =
             tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

     tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

     // create time solver
     tsa::Time_Solver tTimeSolver;

     tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

     sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
     tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
     tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

     tNonlinearSolver.set_dof_type_list( tDofTypes );
     tTimeSolver.set_dof_type_list( tDofTypes );

     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // STEP 4: Solve and check
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     tTimeSolver.solve();

     moris::Matrix< DDRMat > tFinalSolution;

     tTimeSolver.get_full_solution( tFinalSolution );

} //End main_function


//--------------------------------------------------------------------------------------------------------------------//
TEST_CASE( "LinElasticBeam_s", "[moris],[mdl],[LinElasticBeam_s]" )
{
    uint p_size = moris::par_size();

    if( p_size == 1 ) // specify it is a serial test only
    {
        uint tSizeIndicator = 1;
        function_TEST_LinElastic_Beam(tSizeIndicator);
    }
}
//--------------------------------------------------------------------------------------------------------------------//
//TEST_CASE( "LinElasticBeam_m", "[moris],[mdl],[LinElasticBeam_m]" )
//{
//uint p_size = moris::par_size();
//
//if( p_size == 1 ) // specify it is a serial test only
//{
//    uint tSizeIndicator = 2;
//    function_TEST_LinElastic_Beam(tSizeIndicator);
//}
//}
////--------------------------------------------------------------------------------------------------------------------//
//TEST_CASE( "LinElasticBeam_l", "[moris],[mdl],[LinElasticBeam_l]" )
//{
//uint p_size = moris::par_size();
//
//if( p_size == 1 ) // specify it is a serial test only
//{
//    uint tSizeIndicator = 3;
//    function_TEST_LinElastic_Beam(tSizeIndicator);
//}
//}


//--------------------------------------------------------------------------------------------------------------------//
}/* namespace mdl */
}/* namespace moris */



















