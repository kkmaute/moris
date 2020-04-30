/*
 * UT_MDL_Transient.cpp
 *
 *  Created on: Apr 30, 2020
 *      Author: noel
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

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

#include "cl_FEM_IWG_Factory.hpp"             //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"             //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"           //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src
#include "cl_HMR_Database.hpp"

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

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_GEN_Geometry_Field_HMR.hpp"

#include "fn_norm.hpp"


namespace moris
{

// define free function for properties
void tPropConstFunc_MDLTransient
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

bool tSolverOutputCriteria_MDLTransient( moris::tsa::Time_Solver * )
{
    return true;
}

TEST_CASE("MDL Transient","[MDL_Transient]")
{
    if(par_size() == 1)
    {
        uint tLagrangeMeshIndex = 0;

        // empty container for B-Spline meshes
        moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

        // create settings object
        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {1}, {1} } );
        tParameters.set_domain_dimensions( 1, 1 );
        tParameters.set_domain_offset( 0.0, 0.0 );
        tParameters.set_side_sets({ {1}, {2}, {3}, {4} });

        tParameters.set_bspline_truncation( true );
        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns( { {0} });
        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_staircase_buffer( 1 );
        tParameters.set_initial_refinement( 0 );
        tParameters.set_number_aura( true );

        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        // create the HMR object by passing the settings to the constructor
        moris::hmr::HMR tHMR( tParameters );

        tHMR.perform_initial_refinement( 0 );

        tHMR.finalize();

        // construct a mesh manager for the fem
        moris::hmr::Interpolation_Mesh_HMR * tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        moris::hmr::Integration_Mesh_HMR *   tIGMesh = tHMR.create_integration_mesh(1, 0, *tIPMesh );

       // place the pair in mesh manager
       mtk::Mesh_Manager tMeshManager;
       tMeshManager.register_mesh_pair( tIPMesh, tIGMesh );

       //------------------------------------------------------------------------------
       // create the properties
       std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
       tPropConductivity->set_parameters( { {{ 1.0 }} } );
       tPropConductivity->set_val_function( tPropConstFunc_MDLTransient );

       std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
       tPropDensity->set_parameters( { {{ 1.0 }} } );
       tPropDensity->set_val_function( tPropConstFunc_MDLTransient );

       std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
       tPropHeatCapacity->set_parameters( { {{ 1.0 }} } );
       tPropHeatCapacity->set_val_function( tPropConstFunc_MDLTransient );

       std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
       tPropNeumann->set_parameters( { {{ 20.0 }} } );
       tPropNeumann->set_val_function( tPropConstFunc_MDLTransient );

       std::shared_ptr< fem::Property > tPropInitCondition = std::make_shared< fem::Property >();
       tPropInitCondition->set_parameters( { {{ 0.0 }} } );
       tPropInitCondition->set_val_function( tPropConstFunc_MDLTransient );

       std::shared_ptr< fem::Property > tPropWeightCurrent = std::make_shared< fem::Property >();
       tPropWeightCurrent->set_parameters( { {{ 1.0 }} } );
       tPropWeightCurrent->set_val_function( tPropConstFunc_MDLTransient );

       std::shared_ptr< fem::Property > tPropWeightPrevious = std::make_shared< fem::Property >();
       tPropWeightPrevious->set_parameters( { {{ 1.0 }} } );
       tPropWeightPrevious->set_val_function( tPropConstFunc_MDLTransient );

       // define constitutive models
       fem::CM_Factory tCMFactory;

       std::shared_ptr< fem::Constitutive_Model > tCMDiffusion
       = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
       tCMDiffusion->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
       tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
       tCMDiffusion->set_space_dim( 2 );

       // define the IWGs
       fem::IWG_Factory tIWGFactory;

       std::shared_ptr< fem::IWG > tIWGDiffusionBulk
       = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
       tIWGDiffusionBulk->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
       tIWGDiffusionBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
       tIWGDiffusionBulk->set_constitutive_model( tCMDiffusion, "Diffusion", mtk::Master_Slave::MASTER );

       std::shared_ptr< fem::IWG > tIWGNeumann
       = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
       tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
       tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
       tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

       std::shared_ptr< fem::IWG > tIWGTimeContinuity
       = tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
       tIWGTimeContinuity->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
       tIWGTimeContinuity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
       tIWGTimeContinuity->set_property( tPropWeightCurrent, "WeightCurrent", mtk::Master_Slave::MASTER );
       tIWGTimeContinuity->set_property( tPropWeightPrevious, "WeightPrevious", mtk::Master_Slave::MASTER );
       tIWGTimeContinuity->set_property( tPropInitCondition, "InitialCondition", mtk::Master_Slave::MASTER );

       // define the IQIs
       fem::IQI_Factory tIQIFactory;

       std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
       tIQITEMP->set_output_type( vis::Output_Type::TEMP );
       tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP} }, mtk::Master_Slave::MASTER );
       tIQITEMP->set_output_type_index( 0 );

       // define set info
       moris::Cell< fem::Set_User_Info > tSetInfo( 3 );

       tSetInfo( 0 ).set_mesh_index( 0 );
       tSetInfo( 0 ).set_IWGs( { tIWGDiffusionBulk } );
       tSetInfo( 0 ).set_IQIs( { tIQITEMP } );

       tSetInfo( 1 ).set_mesh_index( 2 );
       tSetInfo( 1 ).set_IWGs( { tIWGNeumann } );

       tSetInfo( 2 ).set_mesh_index( 0 );
       tSetInfo( 2 ).set_time_continuity( true );
       tSetInfo( 2 ).set_IWGs( { tIWGTimeContinuity } );

       // create model
       mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo );

       // --------------------------------------------------------------------------------------
       // define outputs
       vis::Output_Manager tOutputData;
       tOutputData.set_outputs( 0,
                                vis::VIS_Mesh_Type::STANDARD,
                                "./",
                                "UT_MDL_Transient.exo",
                                { "HMR_dummy" },
                                { "Temperature" },
                                { vis::Field_Type::NODAL },
                                { vis::Output_Type::TEMP } );
       tModel->set_output_manager( &tOutputData );

       // --------------------------------------------------------------------------------------
       // define linear solver and algorithm
       dla::Solver_Factory  tSolFactory;
       std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
       = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

       dla::Linear_Solver tLinSolver;
       tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

       // --------------------------------------------------------------------------------------
       // define nonlinear solver and algorithm
       NLA::Nonlinear_Solver_Factory tNonlinFactory;
       std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
       = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
       tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

       NLA::Nonlinear_Solver tNonlinearSolver;
       tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

       // --------------------------------------------------------------------------------------
       // define time solver and algorithm
       tsa::Time_Solver_Factory tTimeSolverFactory;
       std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
       = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

       tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
       tTimeSolverAlgorithm->set_param("TSA_Num_Time_Steps")   = 10;
       tTimeSolverAlgorithm->set_param("TSA_Time_Frame")       = 1.0;

       tsa::Time_Solver tTimeSolver;
       tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

       sol::SOL_Warehouse tSolverWarehouse;
       tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

       tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
       tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

       tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
       tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

       tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLTransient );

       // --------------------------------------------------------------------------------------
       // solve and check
       tTimeSolver.solve();

       // clean up
       delete tIPMesh;
       delete tIGMesh;
    }

}/* END_TEST_CASE */
}/* END_MORIS_NAMESPACE */




