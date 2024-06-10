/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_VIS_Output_Manager.cpp
 *
 */

#include "catch.hpp"
#include "moris_typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"

#define protected public
#define private public
#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"
#undef protected
#undef private

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_GEN_Line.hpp"

#include "cl_VIS_Factory.hpp"

#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Reader_Exodus.hpp"

#include "cl_FEM_NodeProxy.hpp"          //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"       //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"    //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"      //FEM/INT/src
#include "cl_FEM_Set.hpp"                //FEM/INT/src

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "fn_PRM_VIS_Parameters.hpp"

inline moris::real
PlaneVisTest( const moris::Matrix< moris::DDRMat >& aPoint )
{
    moris::real mXC = 0.11;
    moris::real mYC = 0.11;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return ( mNx * ( aPoint( 0 ) - mXC ) + mNy * ( aPoint( 1 ) - mYC ) );
}

inline void
tConstValFunction_VISOutputManager( moris::Matrix< moris::DDRMat >& aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >&              aParameters,
        moris::fem::Field_Interpolator_Manager*                     aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

inline bool
tSolverOutputCriteria( moris::tsa::Time_Solver* )
{
    return true;
}

namespace moris
{
    namespace vis
    {
        TEST_CASE( " Output Data", "[VIS],[Output_Data]" )
        {
            if ( par_size() == 1 )
            {
                std::string tFieldName = "Geometry";

                moris::uint tLagrangeMeshIndex = 0;

                moris::hmr::Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 } } );
                tParameters.set_domain_dimensions( { { 2 }, { 1 } } );
                tParameters.set_domain_offset( { { -1.0 }, { -0.0 } } );
                tParameters.set_bspline_truncation( true );

                tParameters.set_output_meshes( { { { 0 } } } );

                tParameters.set_lagrange_orders( { { 1 } } );
                tParameters.set_lagrange_patterns( { { 0 } } );

                tParameters.set_bspline_orders( { { 1 } } );
                tParameters.set_bspline_patterns( { { 0 } } );

                tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );

                tParameters.set_refinement_buffer( 1 );
                tParameters.set_staircase_buffer( 1 );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                hmr::HMR tHMR( tParameters );

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                // create field
                std::shared_ptr< moris::hmr::Field > tPlaneField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tPlaneField->evaluate_scalar_function( PlaneVisTest );

                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tPlaneField );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );
                    tPlaneField->evaluate_scalar_function( PlaneVisTest );
                }

                tHMR.finalize();

                //                tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_enr_ip2.e" );

                hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

                auto tPlane = std::make_shared< moris::gen::Line >( 0.11, 0.11, 1.0, 0.0 );
                Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

                size_t                                tModelDimension = 2;
                moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mGeometries = tGeometryVector;
                moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

                xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
                tXTKModel.mVerbose = false;

                // Specify decomposition Method and Cut Mesh ---------------------------------------
                Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
                tXTKModel.decompose( tDecompositionMethods );

                tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

                xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
                xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

                // place the pair in mesh manager
                std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
                tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

                //------------------------------------------------------------------------------
                // create the properties
                // create the properties
                std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
                tPropLeaderEMod->set_parameters( { { { 1.0 } } } );
                tPropLeaderEMod->set_val_function( tConstValFunction_VISOutputManager );

                std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
                tPropLeaderNu->set_parameters( { { { 0.3 } } } );
                tPropLeaderNu->set_val_function( tConstValFunction_VISOutputManager );

                std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
                tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 } } } );
                tPropDirichlet->set_val_function( tConstValFunction_VISOutputManager );

                std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
                tPropTraction->set_parameters( { { { 1.0 }, { 0.0 } } } );
                tPropTraction->set_val_function( tConstValFunction_VISOutputManager );

                // define constitutive models
                fem::CM_Factory tCMFactory;

                std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso =
                        tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                tCMLeaderElastLinIso->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } } );
                tCMLeaderElastLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
                tCMLeaderElastLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
                tCMLeaderElastLinIso->set_space_dim( 2 );
                tCMLeaderElastLinIso->set_local_properties();

                std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso_bis =
                        tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                tCMLeaderElastLinIso_bis->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } } );
                tCMLeaderElastLinIso_bis->set_property( tPropLeaderEMod, "YoungsModulus" );
                tCMLeaderElastLinIso_bis->set_property( tPropLeaderNu, "PoissonRatio" );
                tCMLeaderElastLinIso_bis->set_space_dim( 2 );
                tCMLeaderElastLinIso_bis->set_local_properties();

                // define stabilization parameters
                fem::SP_Factory tSPFactory;

                std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                        tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
                tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
                tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                        tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
                tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
                tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );
                tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::FOLLOWER );

                // define the IWGs
                fem::IQI_Factory tIQIFactory;

                std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::STRAIN_ENERGY );
                tIQI->set_constitutive_model( tCMLeaderElastLinIso, "Elast", mtk::Leader_Follower::LEADER );
                tIQI->set_name( "IQI" );

                // define the IWGs
                fem::IWG_Factory tIWGFactory;

                std::shared_ptr< fem::IWG > tIWGBulk =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
                tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGBulk->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGBulk->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGDirichlet =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
                tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
                tIWGDirichlet->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
                tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGNeumann =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
                tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGInterface =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
                tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } }, mtk::Leader_Follower::FOLLOWER );
                tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
                tIWGInterface->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
                tIWGInterface->set_constitutive_model( tCMLeaderElastLinIso_bis, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

                // define set info
                fem::Set_User_Info tSetBulk1;
                tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
                tSetBulk1.set_IWGs( { tIWGBulk } );
                tSetBulk1.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetBulk2;
                tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
                tSetBulk2.set_IWGs( { tIWGBulk } );
                tSetBulk2.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetBulk3;
                tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
                tSetBulk3.set_IWGs( { tIWGBulk } );
                tSetBulk3.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetBulk4;
                tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
                tSetBulk4.set_IWGs( { tIWGBulk } );
                tSetBulk4.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetDirichlet;
                tSetDirichlet.set_mesh_set_name( "SideSet_2_n_p1" );
                tSetDirichlet.set_IWGs( { tIWGDirichlet } );

                fem::Set_User_Info tSetNeumann;
                tSetNeumann.set_mesh_set_name( "SideSet_4_n_p0" );
                tSetNeumann.set_IWGs( { tIWGNeumann } );

                // create a list of active block-sets
                std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );

                fem::Set_User_Info tSetInterface;
                tSetInterface.set_mesh_set_name( tDblInterfaceSideSetName );
                tSetInterface.set_IWGs( { tIWGInterface } );

                // create a cell of set info
                Vector< fem::Set_User_Info > tSetInfo( 7 );
                tSetInfo( 0 ) = tSetBulk1;
                tSetInfo( 1 ) = tSetBulk2;
                tSetInfo( 2 ) = tSetBulk3;
                tSetInfo( 3 ) = tSetBulk4;
                tSetInfo( 4 ) = tSetDirichlet;
                tSetInfo( 5 ) = tSetNeumann;
                tSetInfo( 6 ) = tSetInterface;

                // create model
                mdl::Model* tModel = new mdl::Model( tMeshManager,
                        0,
                        tSetInfo );

                // --------------------------------------------------------------------------------------
                // Define outputs

                Output_Manager tOutputData;

                tOutputData.set_outputs( 0,
                        // VIS_Mesh_Type::STANDARD,
                        VIS_Mesh_Type::STANDARD_WITH_OVERLAP,
                        "./",
                        "Output_Vis_Mesh_overlapping.exo",
                        "./",
                        "temp.exo",
                        { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                        { "strain energy elemental", "strain energy global", "strain energy nodal IP" },
                        { Field_Type::ELEMENTAL_AVG, Field_Type::GLOBAL, Field_Type::NODAL },
                        { "IQI", "IQI", "IQI" } );

                tModel->set_output_manager( &tOutputData );

                // --------------------------------------------------------------------------------------
                // Define Solver
                Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
                tDofTypesU( 0 ) = MSI::Dof_Type::UX;
                tDofTypesU( 1 ) = MSI::Dof_Type::UY;

                dla::Solver_Factory                             tSolFactory;
                Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
                tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
                tLinearSolverParameterList.set( "AZ_output", AZ_none );
                tLinearSolverParameterList.set( "AZ_max_iter", 10000 );
                tLinearSolverParameterList.set( "AZ_solver", AZ_gmres );
                tLinearSolverParameterList.set( "AZ_subdomain_solve", AZ_ilu );
                tLinearSolverParameterList.set( "AZ_graph_fill", 10 );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

                dla::Linear_Solver tLinSolver;
                tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 2: create nonlinear solver and algorithm
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                NLA::Nonlinear_Solver_Factory               tNonlinFactory;
                Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
                tNonlinearSolverParameterList.set( "NLA_max_iter", 3 );
                std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

                tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

                NLA::Nonlinear_Solver tNonlinearSolverMain;
                tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

                tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

                // Create solver database
                sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

                tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 3: create time Solver and algorithm
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                tsa::Time_Solver_Factory                      tTimeSolverFactory;
                std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

                tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

                tsa::Time_Solver tTimeSolver;
                tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
                tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

                tTimeSolver.set_dof_type_list( tDofTypesU );

                tTimeSolver.set_output( 0, tSolverOutputCriteria );

                //------------------------------------------------------------------------------
                tTimeSolver.solve();

                delete tInterpMesh;
            }
        }

        TEST_CASE( " Output Data input", "[VIS],[Output_Data_input]" )
        {
            if ( par_size() == 1 )
            {
                std::string tFieldName = "Geometry";

                moris::uint tLagrangeMeshIndex = 0;

                moris::hmr::Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 } } );
                tParameters.set_domain_dimensions( { { 2 }, { 1 } } );
                tParameters.set_domain_offset( { { -1.0 }, { -0.0 } } );
                tParameters.set_bspline_truncation( true );

                tParameters.set_output_meshes( { { { 0 } } } );

                tParameters.set_lagrange_orders( { { 1 } } );
                tParameters.set_lagrange_patterns( { { 0 } } );

                tParameters.set_bspline_orders( { { 1 } } );
                tParameters.set_bspline_patterns( { { 0 } } );

                tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );

                tParameters.set_refinement_buffer( 1 );
                tParameters.set_staircase_buffer( 1 );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                hmr::HMR tHMR( tParameters );

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                // create field
                std::shared_ptr< moris::hmr::Field > tPlaneField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

                tPlaneField->evaluate_scalar_function( PlaneVisTest );

                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tPlaneField );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );
                    tPlaneField->evaluate_scalar_function( PlaneVisTest );
                }

                tHMR.finalize();

                //                tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_enr_ip2.e" );

                hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

                auto tPlane = std::make_shared< moris::gen::Line >( 0.11, 0.11, 1.0, 0.0 );
                Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

                size_t                                tModelDimension = 2;
                moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mGeometries = tGeometryVector;
                moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

                xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
                tXTKModel.mVerbose = false;

                // Specify decomposition Method and Cut Mesh ---------------------------------------
                Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
                tXTKModel.decompose( tDecompositionMethods );

                tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

                xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
                xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

                // place the pair in mesh manager
                std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
                tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

                //------------------------------------------------------------------------------
                // create the properties
                // create the properties
                std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
                tPropLeaderEMod->set_parameters( { { { 1.0 } } } );
                tPropLeaderEMod->set_val_function( tConstValFunction_VISOutputManager );

                std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
                tPropLeaderNu->set_parameters( { { { 0.3 } } } );
                tPropLeaderNu->set_val_function( tConstValFunction_VISOutputManager );

                std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
                tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 } } } );
                tPropDirichlet->set_val_function( tConstValFunction_VISOutputManager );

                std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
                tPropTraction->set_parameters( { { { 1.0 }, { 0.0 } } } );
                tPropTraction->set_val_function( tConstValFunction_VISOutputManager );

                // define constitutive models
                fem::CM_Factory tCMFactory;

                std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso =
                        tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                tCMLeaderElastLinIso->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } } );
                tCMLeaderElastLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
                tCMLeaderElastLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
                tCMLeaderElastLinIso->set_space_dim( 2 );
                tCMLeaderElastLinIso->set_local_properties();
                //            std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso_bis = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                //            tCMLeaderElastLinIso_bis->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
                //            tCMLeaderElastLinIso_bis->set_property( tPropLeaderEMod, "YoungsModulus" );
                //            tCMLeaderElastLinIso_bis->set_property( tPropLeaderNu, "PoissonRatio" );
                //            tCMLeaderElastLinIso_bis->set_space_dim( 2 );
                std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso_bis =
                        tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                tCMLeaderElastLinIso_bis->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tCMLeaderElastLinIso_bis->set_property( tPropLeaderEMod, "YoungsModulus" );
                tCMLeaderElastLinIso_bis->set_property( tPropLeaderNu, "PoissonRatio" );
                tCMLeaderElastLinIso_bis->set_space_dim( 2 );
                tCMLeaderElastLinIso_bis->set_local_properties();

                // define stabilization parameters
                fem::SP_Factory tSPFactory;

                std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                        tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
                tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
                tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                        tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
                tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
                tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );
                tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::FOLLOWER );

                // define the IWGs
                fem::IQI_Factory tIQIFactory;

                std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::STRAIN_ENERGY );
                tIQI->set_constitutive_model( tCMLeaderElastLinIso, "Elast", mtk::Leader_Follower::LEADER );
                tIQI->set_name( "IQI" );

                // define the IWGs
                fem::IWG_Factory tIWGFactory;

                std::shared_ptr< fem::IWG > tIWGBulk =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
                tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGBulk->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGBulk->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGDirichlet =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
                tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
                tIWGDirichlet->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
                tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGNeumann =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
                tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGInterface =
                        tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
                tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } }, mtk::Leader_Follower::FOLLOWER );
                tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
                tIWGInterface->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
                tIWGInterface->set_constitutive_model( tCMLeaderElastLinIso_bis, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

                // define set info
                fem::Set_User_Info tSetBulk1;
                tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
                tSetBulk1.set_IWGs( { tIWGBulk } );
                tSetBulk1.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetBulk2;
                tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
                tSetBulk2.set_IWGs( { tIWGBulk } );
                tSetBulk2.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetBulk3;
                tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
                tSetBulk3.set_IWGs( { tIWGBulk } );
                tSetBulk3.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetBulk4;
                tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
                tSetBulk4.set_IWGs( { tIWGBulk } );
                tSetBulk4.set_IQIs( { tIQI } );

                fem::Set_User_Info tSetDirichlet;
                tSetDirichlet.set_mesh_set_name( "SideSet_2_n_p1" );
                tSetDirichlet.set_IWGs( { tIWGDirichlet } );

                fem::Set_User_Info tSetNeumann;
                tSetNeumann.set_mesh_set_name( "SideSet_4_n_p0" );
                tSetNeumann.set_IWGs( { tIWGNeumann } );

                // create a list of active block-sets
                std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );

                fem::Set_User_Info tSetInterface;
                tSetInterface.set_mesh_set_name( tDblInterfaceSideSetName );
                tSetInterface.set_IWGs( { tIWGInterface } );

                // create a cell of set info
                Vector< fem::Set_User_Info > tSetInfo( 7 );
                tSetInfo( 0 ) = tSetBulk1;
                tSetInfo( 1 ) = tSetBulk2;
                tSetInfo( 2 ) = tSetBulk3;
                tSetInfo( 3 ) = tSetBulk4;
                tSetInfo( 4 ) = tSetDirichlet;
                tSetInfo( 5 ) = tSetNeumann;
                tSetInfo( 6 ) = tSetInterface;

                // create model
                mdl::Model* tModel = new mdl::Model( tMeshManager, 0, tSetInfo );

                // --------------------------------------------------------------------------------------
                // Define outputs
                moris::Parameter_List tParameterList = moris::prm::create_vis_parameter_list();

                std::string tMorisOutput = std::getenv( "MORISOUTPUT" );

                MORIS_ERROR( tMorisOutput.size() > 0,
                        "Environment variable MORISOUTPUT not set." );

                tParameterList.set( "File_Name", std::pair< std::string, std::string >( tMorisOutput, "Vis_Test.exo" ) );
                tParameterList.set( "Set_Names", std::string( "HMR_dummy_c_p0,HMR_dummy_c_p1,HMR_dummy_n_p0,HMR_dummy_n_p1" ) );
                tParameterList.set( "Field_Names", std::string( "strain_energy_elemental,strain_energy_global,strain_energy_nodal_IP" ) );
                tParameterList.set( "Field_Type", std::string( "ELEMENTAL_INT,GLOBAL,NODAL" ) );
                tParameterList.set( "IQI_Names", std::string( "IQI,IQI,IQI" ) );

                Output_Manager tOutputData( tParameterList );

                tModel->set_output_manager( &tOutputData );

                // --------------------------------------------------------------------------------------
                // Define Solver
                Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
                tDofTypesU( 0 ) = MSI::Dof_Type::UX;
                tDofTypesU( 1 ) = MSI::Dof_Type::UY;

                dla::Solver_Factory                             tSolFactory;
                Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
                tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
                tLinearSolverParameterList.set( "AZ_output", AZ_none );
                tLinearSolverParameterList.set( "AZ_max_iter", 10000 );
                tLinearSolverParameterList.set( "AZ_solver", AZ_gmres );
                tLinearSolverParameterList.set( "AZ_subdomain_solve", AZ_ilu );
                tLinearSolverParameterList.set( "AZ_graph_fill", 10 );
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );
                //        tLinearSolverParameterList.set( "ml_prec_type", "SA" );

                dla::Linear_Solver tLinSolver;
                tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 2: create nonlinear solver and algorithm
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                NLA::Nonlinear_Solver_Factory               tNonlinFactory;
                Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
                tNonlinearSolverParameterList.set( "NLA_max_iter", 3 );
                std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

                tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

                NLA::Nonlinear_Solver tNonlinearSolverMain;
                tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

                tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

                // Create solver database
                sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

                tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 3: create time Solver and algorithm
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                tsa::Time_Solver_Factory                      tTimeSolverFactory;
                std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

                tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

                tsa::Time_Solver tTimeSolver;
                tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
                tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

                tTimeSolver.set_dof_type_list( tDofTypesU );

                tTimeSolver.set_output( 0, tSolverOutputCriteria );

                //------------------------------------------------------------------------------
                tTimeSolver.solve();

                delete tInterpMesh;
            }
        }
    }    // namespace vis
}    // namespace moris
