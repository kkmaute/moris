/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_FEM_DQ_Dp.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "moris_typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
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

#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
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
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"

#include "../projects/FEM/MSI/test/MSI_Test_Proxy/cl_MSI_Design_Variable_Interface_Proxy.hpp"

#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

namespace moris
{

void tPropValConstFunc_MDLFEM_DQ_DP
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

moris::real tPlane_MDLFEM_DQ_DP( const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tOffset = 2.6;
    return    aPoint(0) - 0.0 * aPoint(1) - tOffset;
//    return    aPoint(0) - 0.317 * aPoint(1) - tOffset;
}

TEST_CASE("MDL FEM Elastic DQ/Dp","[MDL_FEM_DQ_DP]")
{
//    if(par_size() == 1)
//    {
////	std::cout<<"I am proc: "<<par_rank()<<std::endl;
//
//        uint tLagrangeMeshIndex = 0;
//
//        // empty container for B-Spline meshes
//        Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;
//
//        // create settings object
//        moris::hmr::Parameters tParameters;
//
//        tParameters.set_number_of_elements_per_dimension( { {4}, {2}, {2} } );
//        tParameters.set_domain_dimensions( 10, 5, 5 );
//        tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
//        tParameters.set_side_sets({ {1}, {2}, {3}, {4}, {5}, {6} });
//
//        tParameters.set_bspline_truncation( true );
//
//        tParameters.set_lagrange_orders  ( { {1} });
//        tParameters.set_lagrange_patterns( { {0} });
//
//        tParameters.set_bspline_orders   ( { {1} } );
//        tParameters.set_bspline_patterns ( { {0} } );
//
//        tParameters.set_output_meshes( { {0} } );
////        tParameters.set_lagrange_input_mesh( { { 0 } } );
//
//        tParameters.set_staircase_buffer( 1 );
//
//        tParameters.set_initial_refinement( 0 );
//
//        tParameters.set_number_aura( true );
//
//        Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//        tLagrangeToBSplineMesh( 0 ) = { {0} };
//
//        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//        // create the HMR object by passing the settings to the constructor
//        moris::hmr::HMR tHMR( tParameters );
//
//        tHMR.perform_initial_refinement( 0 );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );   // HMR Lagrange mesh
//        //==============================
//        std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex);
//
//        tField->evaluate_scalar_function( tPlane_MDLFEM_DQ_DP );
//
//        Vector< std::shared_ptr< moris::hmr::Field > > tFields( 1, tField );
//
//        for( uint k=0; k<0; ++k )
//        {
//            tHMR.flag_surface_elements_on_working_pattern( tField );
////            tHMR.user_defined_flagging( user_defined_refinement_MDLFEMBench, tFields, tParam, 0 );
//            tHMR.perform_refinement_based_on_working_pattern( 0, true );
//            tField->evaluate_scalar_function( tPlane_MDLFEM_DQ_DP );
//        }
//        tHMR.finalize();
//
////==============================
//        tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );
//
//        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);
//
//        Vector<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//        size_t tModelDimension = 3;
//        moris::ge::GEN_Phase_Table  tPhaseTable( tGeometryVector.size());
//        moris::ge::Geometry_Engine  tGeometryEngine( tGeometryVector,tPhaseTable,tModelDimension );
//
//        xtk::Model tXTKModel( tModelDimension,tInterpMesh.get(),&tGeometryEngine );
//        tXTKModel.mVerbose = false;
//
//        //Specify decomposition Method and Cut Mesh ---------------------------------------
//        Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//        tXTKModel.decompose(tDecompositionMethods);
//
//       tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
//
//       xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//       xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//       //==============================
//
//       // place the pair in mesh manager
//       mtk::Mesh_Manager tMeshManager;
//       tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
//
//       //------------------------------------------------------------------------------
//       // create the properties
//       std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
//       tPropEMod1->set_parameters( { {{ 1.0 }} } );
//       tPropEMod1->set_val_function( tPropValConstFunc_MDLFEM_DQ_DP );
//
//       std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
//       tPropEMod2->set_parameters( { {{ 1.0 }} } );
//       tPropEMod2->set_val_function( tPropValConstFunc_MDLFEM_DQ_DP );
//
//       std::shared_ptr< fem::Property > tPropNu1 = std::make_shared< fem::Property >();
//       tPropNu1->set_parameters( { {{ 0.3 }} } );
//       tPropNu1->set_val_function( tPropValConstFunc_MDLFEM_DQ_DP );
//
//       std::shared_ptr< fem::Property > tPropNu2 = std::make_shared< fem::Property >();
//       tPropNu2->set_parameters( { {{ 0.3 }} } );
//       tPropNu2->set_val_function( tPropValConstFunc_MDLFEM_DQ_DP );
//
//       std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
//       tPropDirichlet->set_parameters( { {{ 0.0 }, { 0.0 }, { 0.0 }} } );
//       tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEM_DQ_DP );
//
//       std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
//       tPropNeumann->set_parameters( {{{ 1.0 }, { 0.0 }, { 0.0 }}} );
//       tPropNeumann->set_val_function( tPropValConstFunc_MDLFEM_DQ_DP );
//
//       // working do types
//       Vector< moris::MSI::Dof_Type > tResDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ };
//
//       // define constitutive models
//       fem::CM_Factory tCMFactory;
//
//       std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
//       tCMStrucLinIso1->set_dof_type_list( { tResDofTypes } );
//       tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
//       tCMStrucLinIso1->set_property( tPropNu1, "PoissonRatio" );
//       tCMStrucLinIso1->set_model_type( fem::Model_Type::FULL );
//       tCMStrucLinIso1->set_space_dim( 3 );
//
//       std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
//       tCMStrucLinIso2->set_dof_type_list( { tResDofTypes } );
//       tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
//       tCMStrucLinIso2->set_property( tPropNu2, "PoissonRatio" );
//       tCMStrucLinIso2->set_model_type( fem::Model_Type::FULL );
//       tCMStrucLinIso2->set_space_dim( 3 );
//
//       // define stabilization parameters
//       fem::SP_Factory tSPFactory;
//       std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
//       tSPDirichletNitsche->set_parameters( { {{ 100.0 }} } );
//       tSPDirichletNitsche->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );
//
//       std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
//       tSPNitscheInterface->set_parameters( { {{ 100.0 }} } );
//       tSPNitscheInterface->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );
//       tSPNitscheInterface->set_property( tPropEMod2, "Material", mtk::Leader_Follower::FOLLOWER );
//
//       // define the IWGs
//       fem::IWG_Factory tIWGFactory;
//
//       std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
//       tIWGBulk1->set_residual_dof_type( tResDofTypes );
//       tIWGBulk1->set_dof_type_list( { tResDofTypes } );
//       tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
//
//       std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
//       tIWGBulk2->set_residual_dof_type( tResDofTypes );
//       tIWGBulk2->set_dof_type_list( { tResDofTypes } );
//       tIWGBulk2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );
//
//       std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
//       tIWGDirichlet->set_residual_dof_type( tResDofTypes );
//       tIWGDirichlet->set_dof_type_list( { tResDofTypes } );
//       tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
//       tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
//       tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
//
//       std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
//       tIWGNeumann->set_residual_dof_type( tResDofTypes );
//       tIWGNeumann->set_dof_type_list( { tResDofTypes } );
//       tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );
//
//       std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
//       tIWGInterface->set_residual_dof_type( tResDofTypes );
//       tIWGInterface->set_dof_type_list( { tResDofTypes } );
//       tIWGInterface->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::FOLLOWER );
//       tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
//       tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
//       tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
//
//       // create the IQIs
//       fem::IQI_Factory tIQIFactory;
//
//       std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
//       tIQIUX->set_output_type( vis::Output_Type::UX );
//       tIQIUX->set_dof_type_list( { { tResDofTypes } }, mtk::Leader_Follower::LEADER );
//       tIQIUX->set_output_type_index( 0 );
//
//       std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
//       tIQIUY->set_output_type( vis::Output_Type::UY );
//       tIQIUY->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
//       tIQIUY->set_output_type_index( 1 );
//
//       std::shared_ptr< fem::IQI > tIQIUZ = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
//       tIQIUZ->set_output_type( vis::Output_Type::UZ );
//       tIQIUZ->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
//       tIQIUZ->set_output_type_index( 2 );
//
//       // define set info
//       fem::Set_User_Info tSetBulk1;
//       tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
//       tSetBulk1.set_IWGs( { tIWGBulk1 } );
//       tSetBulk1.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );
//
//       fem::Set_User_Info tSetBulk2;
//       tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
//       tSetBulk2.set_IWGs( { tIWGBulk1 } );
//       tSetBulk2.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );
//
//       fem::Set_User_Info tSetBulk3;
//       tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
//       tSetBulk3.set_IWGs( { tIWGBulk2 } );
//       tSetBulk3.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );
//
//       fem::Set_User_Info tSetBulk4;
//       tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
//       tSetBulk4.set_IWGs( { tIWGBulk2 } );
//       tSetBulk4.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );
//
//       fem::Set_User_Info tSetDirichlet;
//       tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p0" );
//       tSetDirichlet.set_IWGs( { tIWGDirichlet } );
//
//       fem::Set_User_Info tSetNeumann;
//       tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
//       tSetNeumann.set_IWGs( { tIWGNeumann } );
//
//       fem::Set_User_Info tSetInterface;
//       tSetInterface.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 ) );
//       tSetInterface.set_IWGs( { tIWGInterface } );
//
//       // create a cell of set info
//       Vector< fem::Set_User_Info > tSetInfo( 7 );
//       tSetInfo( 0 ) = tSetBulk1;
//       tSetInfo( 1 ) = tSetBulk2;
//       tSetInfo( 2 ) = tSetBulk3;
//       tSetInfo( 3 ) = tSetBulk4;
//       tSetInfo( 4 ) = tSetDirichlet;
//       tSetInfo( 5 ) = tSetNeumann;
//       tSetInfo( 6 ) = tSetInterface;
//
//       MSI::Design_Variable_Interface * tDesignVariableInterface = new MSI::Design_Variable_Interface_Proxy();
//
//       // create model
//       mdl::Model * tModel = new mdl::Model( &tMeshManager,
//                                              0,
//                                              tSetInfo,
//                                              tDesignVariableInterface );
//
//       // --------------------------------------------------------------------------------------
//       // define outputs
//       vis::Output_Manager tOutputData;
//       tOutputData.set_outputs( 0,
//                                vis::VIS_Mesh_Type::STANDARD,
//                                "UT_MDL_FEM_Benchmark_Elast_Interface_Output.exo",
//                                { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
//                                { "Displacement UX", "Displacement UY", "Displacement UZ" },
//                                { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
//                                { vis::Output_Type::UX, vis::Output_Type::UY, vis::Output_Type::UZ } );
//
//       tModel->set_output_manager( &tOutputData );
//
//       // --------------------------------------------------------------------------------------
//       // define linear solver and algorithm
//       dla::Solver_Factory  tSolFactory;
//       std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
//       = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
//
////       tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_all;
////       tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
////       tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;
//
//       dla::Linear_Solver tLinSolver;
//
//       tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//       // --------------------------------------------------------------------------------------
//       // define nonlinear solver and algorithm
//       NLA::Nonlinear_Solver_Factory tNonlinFactory;
//       std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
//       = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//       tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//       NLA::Nonlinear_Solver tNonlinearSolver;
//       tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//       // --------------------------------------------------------------------------------------
//       // define time solver and algorithm
//       tsa::Time_Solver_Factory tTimeSolverFactory;
//       std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
//       = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//       tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//       tsa::Time_Solver tTimeSolver;
//
//       tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//       sol::SOL_Warehouse tSolverWarehouse;
//
//       tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//       tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//       tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//       tNonlinearSolver.set_dof_type_list( { tResDofTypes } );
//       tTimeSolver.set_dof_type_list( { tResDofTypes } );
//
////       tTimeSolver.set_output( 0, nullptr );
//
//       // --------------------------------------------------------------------------------------
//
//       // solve and check
//       tTimeSolver.solve();
//
//       MSI::MSI_Solver_Interface * tSolverInterface = tModel->get_solver_interface();
//
//       tSolverInterface->set_is_forward( false );
//
//       tTimeSolver.solve();
//
//       Matrix<DDRMat> tFullSol;
//       tTimeSolver.get_full_solution(tFullSol);
//    }

}/* END_TEST_CASE */

}/* END_MORIS_NAMESPACE */

