/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_FEM_Benchmark.cpp
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

#include "cl_Matrix.hpp"    //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"    // ALG/src

#include "cl_FEM_IWG_Factory.hpp"                   //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"                   //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"                 //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src
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
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "cl_GEN_Line.hpp"

#include "fn_norm.hpp"

namespace moris
{

    // define free function for properties
    inline void
    tPropValConstFunc_MDLFEMBench( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&         aParameters,
            moris::fem::Field_Interpolator_Manager*                aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    inline void
    tPropValFuncL2_MDLFEMBench( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&      aParameters,
            moris::fem::Field_Interpolator_Manager*             aFIManager )
    {
        aPropMatrix = { { 20 * aFIManager->get_IP_geometry_interpolator()->valx()( 0 ) } };
    }

    // define function for cutting plane
    inline moris::real
    tPlane_MDLFEMBench( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real tOffset = 2.6;
        return aPoint( 0 ) - tOffset;
    }

    inline bool
    tSolverOutputCriteria_MDLFEMBench( moris::tsa::Time_Solver* )
    {
        return true;
    }

    TEST_CASE( "MDL FEM Benchmark Diff Block", "[MDL_FEM_Benchmark_Diff_Block]" )
    {
        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 }, { 2 } } );
            tParameters.set_domain_dimensions( 10, 5, 5 );
            tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

            tParameters.set_bspline_truncation( true );
            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );
            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );
            tParameters.set_initial_refinement( { { 0 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );
            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            tHMR.finalize();

            // construct a mesh manager for the fem
            moris::hmr::Interpolation_Mesh_HMR* tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            moris::hmr::Integration_Mesh_HMR*   tIGMesh = tHMR.create_integration_mesh( 1, 0, tIPMesh );

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( tIPMesh, tIGMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
            tPropConductivity1->set_parameters( { { { 1.0 } } } );
            tPropConductivity1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 } } } );
            tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropFlux = std::make_shared< fem::Property >();
            tPropFlux->set_parameters( { { { 20.0 } } } );
            tPropFlux->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
            tPropTempLoad1->set_parameters( { { { 0.0 } } } );
            tPropTempLoad1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropL2Analytic = std::make_shared< fem::Property >();
            tPropL2Analytic->set_val_function( tPropValFuncL2_MDLFEMBench );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 3 );
            tCMDiffLinIso1->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropFlux, "Neumann", mtk::Leader_Follower::LEADER );

            // define the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            std::shared_ptr< fem::IQI > tIQIL2TEMP = tIQIFactory.create_IQI( fem::IQI_Type::L2_ERROR_ANALYTIC );
            tIQIL2TEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQIL2TEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQIL2TEMP->set_property( tPropL2Analytic, "L2Check", mtk::Leader_Follower::LEADER );
            tIQIL2TEMP->set_name( "IQI_L2" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_index( 0 );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );
            tSetBulk1.set_IQIs( { tIQITEMP, tIQIL2TEMP } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_index( 4 );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_index( 2 );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 3 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetDirichlet;
            tSetInfo( 2 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_FEM_Benchmark_Output_Block.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy" },
                    { "Temperature", "L2 error" },
                    { vis::Field_Type::NODAL, vis::Field_Type::GLOBAL },
                    { "IQI_TEMP", "IQI_L2" } );
            tModel->set_output_manager( &tOutputData );

            // --------------------------------------------------------------------------------------
            // define linear solver and algorithm
            dla::Solver_Factory                             tSolFactory;
            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_amesos();
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // --------------------------------------------------------------------------------------
            // define nonlinear solver and algorithm
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // --------------------------------------------------------------------------------------
            // define time solver and algorithm
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::TEMP } );
            tTimeSolver.set_dof_type_list( { MSI::Dof_Type::TEMP } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench );

            // --------------------------------------------------------------------------------------
            // solve and check
            tTimeSolver.solve();
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );
            delete tIPMesh;
            delete tIGMesh;
        }

    } /* END_TEST_CASE */

    TEST_CASE( "MDL FEM Benchmark Diff Interface", "[MDL_FEM_Benchmark_Diff_Interface]" )
    {

        if ( par_size() == 1 )
        {
            //	std::cout<<"I am proc: "<<par_rank()<<std::endl;

            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 }, { 2 } } );
            tParameters.set_domain_dimensions( 10, 5, 5 );
            tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );

            tParameters.set_initial_refinement( { { 0 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );

            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );    // HMR Lagrange mesh
            //==============================
            std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex );

            tField->evaluate_scalar_function( tPlane_MDLFEMBench );

            Vector< std::shared_ptr< moris::hmr::Field > > tFields( 1, tField );

            // FIXME what is the following test about
            if ( false )
            {
                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tField );
                    //            tHMR.user_defined_flagging( user_defined_refinement_MDLFEMBench, tFields, tParam, 0 );
                    tHMR.perform_refinement_based_on_working_pattern( 0, true );
                    tField->evaluate_scalar_function( tPlane_MDLFEMBench );
                }
            }

            tHMR.finalize();

            //==============================
            tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tPlane = std::make_shared< moris::gen::Line >( 2.6, 0.0, 1.0, 0.0 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

            size_t                                tModelDimension = 3;
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            //==============================

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
            tPropConductivity1->set_parameters( { { { 1.0 } } } );
            tPropConductivity1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
            tPropConductivity2->set_parameters( { { { 1.0 } } } );
            tPropConductivity2->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 } } } );
            tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropFlux = std::make_shared< fem::Property >();
            tPropFlux->set_parameters( { { { 20.0 } } } );
            tPropFlux->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
            tPropTempLoad1->set_parameters( { { { 0.0 } } } );
            tPropTempLoad1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
            tPropTempLoad2->set_parameters( { { { 0.0 } } } );
            tPropTempLoad2->set_val_function( tPropValConstFunc_MDLFEMBench );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 3 );
            tCMDiffLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
            tCMDiffLinIso2->set_space_dim( 3 );
            tCMDiffLinIso2->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
            tSPNitscheInterface->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );
            tSPNitscheInterface->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::FOLLOWER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGBulk2->set_property( tPropTempLoad2, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropFlux, "Neumann", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
            tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
            tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
            tIWGInterface->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGInterface->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::FOLLOWER );

            // create the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk2 } );
            tSetBulk1.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk2 } );
            tSetBulk2.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk3;
            tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk3.set_IWGs( { tIWGBulk1 } );
            tSetBulk3.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk4;
            tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk4.set_IWGs( { tIWGBulk1 } );
            tSetBulk4.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p0" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );

            fem::Set_User_Info tSetInterface1;
            tSetInterface1.set_mesh_set_name( tDblInterfaceSideSetName );
            tSetInterface1.set_IWGs( { tIWGInterface } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 7 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetBulk3;
            tSetInfo( 3 ) = tSetBulk4;
            tSetInfo( 4 ) = tSetDirichlet;
            tSetInfo( 5 ) = tSetNeumann;
            tSetInfo( 6 ) = tSetInterface1;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // Define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_FEM_Benchmark_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                    { "Temperature" },
                    { vis::Field_Type::NODAL },
                    { "IQI_TEMP" } );

            tModel->set_output_manager( &tOutputData );

            dla::Solver_Factory                             tSolFactory;
            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
            tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
            tLinearSolverParameterList.set( "AZ_output", AZ_all );
            tLinearSolverParameterList.set( "AZ_solver", AZ_gmres_condnum );
            tLinearSolverParameterList.set( "AZ_precond", AZ_none );
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );
            // tLinearSolverParameterList.set( "AZ_kspace", 500 );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create time Solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::TEMP } );
            tTimeSolver.set_dof_type_list( { MSI::Dof_Type::TEMP } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tTimeSolver.solve();
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );

            delete tInterpMesh;
        }

    } /* END_TEST_CASE */

    TEST_CASE( "MDL FEM Benchmark Diff Ghost", "[MDL_FEM_Benchmark_Diff_Ghost]" )
    {
        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 }, { 2 } } );
            tParameters.set_domain_dimensions( 10, 5, 5 );
            tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );

            tParameters.set_initial_refinement( { { 0 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );

            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );    // HMR Lagrange mesh
            //==============================
            std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex );

            tField->evaluate_scalar_function( tPlane_MDLFEMBench );

            Vector< std::shared_ptr< moris::hmr::Field > > tFields( 1, tField );

            // FIXME what is the following test about
            if ( false )
            {
                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tField );
                    //            tHMR.user_defined_flagging( user_defined_refinement_MDLFEMBench, tFields, tParam, 0 );
                    tHMR.perform_refinement_based_on_working_pattern( 0, true );
                    tField->evaluate_scalar_function( tPlane_MDLFEMBench );
                }
            }

            tHMR.finalize();

            //==============================
            tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tPlane = std::make_shared< moris::gen::Line >( 2.6, 0.0, 1.0, 0.0 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

            size_t                                tModelDimension = 3;
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
            tXTKModel.construct_face_oriented_ghost_penalization_cells();

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            moris_index tSSIndex = tEnrIntegMesh.create_side_set_from_dbl_side_set( 1, "ghost_ss_p0" );
            tEnrIntegMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_p0", mtk::CellTopology::QUAD4 );

            //==============================

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
            tPropConductivity1->set_parameters( { { { 1.0 } } } );
            tPropConductivity1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
            tPropConductivity2->set_parameters( { { { 1.0 } } } );
            tPropConductivity2->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 } } } );
            tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { 20.0 } } } );
            tPropNeumann->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
            tPropTempLoad1->set_parameters( { { { 0.0 } } } );
            tPropTempLoad1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
            tPropTempLoad2->set_parameters( { { { 0.0 } } } );
            tPropTempLoad2->set_val_function( tPropValConstFunc_MDLFEMBench );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 3 );
            tCMDiffLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 =
                    tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
            tCMDiffLinIso2->set_space_dim( 3 );
            tCMDiffLinIso2->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
            tSPNitscheInterface->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );
            tSPNitscheInterface->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::FOLLOWER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPGhost =
                    tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
            tSPGhost->set_parameters( { { { 0.1 } } } );
            tSPGhost->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGBulk2->set_property( tPropTempLoad2, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
            tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
            tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
            tIWGInterface->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGInterface->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::FOLLOWER );

            // Ghost stabilization
            std::shared_ptr< fem::IWG > tIWGGhost = tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
            tIWGGhost->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGGhost->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGGhost->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
            tIWGGhost->set_stabilization_parameter( tSPGhost, "GhostSP" );

            // create the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk2 } );
            tSetBulk1.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk2 } );
            tSetBulk2.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk3;
            tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk3.set_IWGs( { tIWGBulk1 } );
            tSetBulk3.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk4;
            tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk4.set_IWGs( { tIWGBulk1 } );
            tSetBulk4.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p0" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );

            fem::Set_User_Info tSetInterface1;
            tSetInterface1.set_mesh_set_name( tDblInterfaceSideSetName );
            tSetInterface1.set_IWGs( { tIWGInterface } );

            fem::Set_User_Info tSetGhost;
            tSetGhost.set_mesh_set_name( "ghost_p0" );
            tSetGhost.set_IWGs( { tIWGGhost } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 8 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetBulk3;
            tSetInfo( 3 ) = tSetBulk4;
            tSetInfo( 4 ) = tSetDirichlet;
            tSetInfo( 5 ) = tSetNeumann;
            tSetInfo( 6 ) = tSetInterface1;
            tSetInfo( 7 ) = tSetGhost;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // Define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_FEM_Benchmark_Ghost_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                    { "Temperature" },
                    { vis::Field_Type::NODAL },
                    { "IQI_TEMP" } );

            tModel->set_output_manager( &tOutputData );

            dla::Solver_Factory                             tSolFactory;
            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
            tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
            tLinearSolverParameterList.set( "AZ_output", AZ_all );
            tLinearSolverParameterList.set( "AZ_solver", AZ_gmres_condnum );
            tLinearSolverParameterList.set( "AZ_precond", AZ_none );
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 3: create time Solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::TEMP } );
            tTimeSolver.set_dof_type_list( { MSI::Dof_Type::TEMP } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            tTimeSolver.solve();
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );

            delete tInterpMesh;
        }

    } /* END_TEST_CASE */

    TEST_CASE( "MDL FEM Benchmark Elast Block", "[MDL_FEM_Benchmark_Elast_Block]" )
    {
        if ( par_size() == 1 )
        {
            //	std::cout<<"I am proc: "<<par_rank()<<std::endl;

            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 }, { 2 } } );
            tParameters.set_domain_dimensions( 10, 5, 5 );
            tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

            tParameters.set_bspline_truncation( true );
            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );
            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );
            tParameters.set_initial_refinement( { { 0 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );
            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            tHMR.finalize();

            // construct a mesh manager for the fem
            moris::hmr::Interpolation_Mesh_HMR* tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            moris::hmr::Integration_Mesh_HMR*   tIGMesh = tHMR.create_integration_mesh( 1, 0, tIPMesh );

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( tIPMesh, tIGMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
            tPropEMod1->set_parameters( { { { 1.0 } } } );
            tPropEMod1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropNu1 = std::make_shared< fem::Property >();
            tPropNu1->set_parameters( { { { 0.3 } } } );
            tPropNu1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 }, { 0.0 } } } );
            tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
            tPropTraction->set_parameters( { { { 1.0 }, { 0.0 }, { 0.0 } } } );
            tPropTraction->set_val_function( tPropValConstFunc_MDLFEMBench );

            // working do types
            Vector< moris::MSI::Dof_Type > tResDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ };

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { tResDofTypes } );
            tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropNu1, "PoissonRatio" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::FULL );
            tCMStrucLinIso1->set_space_dim( 3 );
            tCMStrucLinIso1->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk1->set_residual_dof_type( { tResDofTypes } );
            tIWGBulk1->set_dof_type_list( { tResDofTypes } );
            tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { tResDofTypes } );
            tIWGDirichlet->set_dof_type_list( { tResDofTypes } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { tResDofTypes } );
            tIWGNeumann->set_dof_type_list( { tResDofTypes } );
            tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

            // create the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUX->set_quantity_dof_type( tResDofTypes );
            tIQIUX->set_dof_type_list( { { tResDofTypes } }, mtk::Leader_Follower::LEADER );
            tIQIUX->set_output_type_index( 0 );
            tIQIUX->set_name( "IQI_UX" );

            std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUY->set_quantity_dof_type( tResDofTypes );
            tIQIUY->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
            tIQIUY->set_output_type_index( 1 );
            tIQIUY->set_name( "IQI_UY" );

            std::shared_ptr< fem::IQI > tIQIUZ = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUZ->set_quantity_dof_type( tResDofTypes );
            tIQIUZ->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
            tIQIUZ->set_output_type_index( 2 );
            tIQIUZ->set_name( "IQI_UZ" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_index( 0 );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );
            tSetBulk1.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_index( 4 );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_index( 2 );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 3 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetDirichlet;
            tSetInfo( 2 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_FEM_Benchmark_Elast_Block_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy" },
                    { "Displacement UX", "Displacement UY", "Displacement UZ" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_UX", "IQI_UY", "IQI_UZ" } );

            tModel->set_output_manager( &tOutputData );

            // --------------------------------------------------------------------------------------
            // define linear solver and algorithm
            dla::Solver_Factory                             tSolFactory;
            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_amesos();
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // --------------------------------------------------------------------------------------
            // define nonlinear solver and algorithm
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // --------------------------------------------------------------------------------------
            // define time solver and algorithm
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { tResDofTypes } );
            tTimeSolver.set_dof_type_list( { tResDofTypes } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench );

            // --------------------------------------------------------------------------------------
            // solve and check
            tTimeSolver.solve();
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );

            delete tIPMesh;
            delete tIGMesh;
        }

    } /* END_TEST_CASE */

    TEST_CASE( "MDL FEM Benchmark Elast Interface", "[MDL_FEM_Benchmark_Elast_Interface]" )
    {
        if ( par_size() == 1 )
        {
            //	std::cout<<"I am proc: "<<par_rank()<<std::endl;

            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 }, { 2 } } );
            tParameters.set_domain_dimensions( 10, 5, 5 );
            tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );

            tParameters.set_initial_refinement( { { 0 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );

            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );    // HMR Lagrange mesh
            //==============================
            std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex );

            tField->evaluate_scalar_function( tPlane_MDLFEMBench );

            Vector< std::shared_ptr< moris::hmr::Field > > tFields( 1, tField );

            // FIXME what is the following test about
            if ( false )
            {
                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tField );
                    //            tHMR.user_defined_flagging( user_defined_refinement_MDLFEMBench, tFields, tParam, 0 );
                    tHMR.perform_refinement_based_on_working_pattern( 0, true );
                    tField->evaluate_scalar_function( tPlane_MDLFEMBench );
                }
            }

            tHMR.finalize();

            //==============================
            tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tPlane = std::make_shared< moris::gen::Line >( 2.6, 0.0, 1.0, 0.0 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

            size_t                                tModelDimension = 3;
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            //==============================

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
            tPropEMod1->set_parameters( { { { 1.0 } } } );
            tPropEMod1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
            tPropEMod2->set_parameters( { { { 1.0 } } } );
            tPropEMod2->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropNu1 = std::make_shared< fem::Property >();
            tPropNu1->set_parameters( { { { 0.3 } } } );
            tPropNu1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropNu2 = std::make_shared< fem::Property >();
            tPropNu2->set_parameters( { { { 0.3 } } } );
            tPropNu2->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 }, { 0.0 } } } );
            tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
            tPropTraction->set_parameters( { { { 1.0 }, { 0.0 }, { 0.0 } } } );
            tPropTraction->set_val_function( tPropValConstFunc_MDLFEMBench );

            // working do types
            Vector< moris::MSI::Dof_Type > tResDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ };

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { tResDofTypes } );
            tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropNu1, "PoissonRatio" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::FULL );
            tCMStrucLinIso1->set_space_dim( 3 );
            tCMStrucLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso2->set_dof_type_list( { tResDofTypes } );
            tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
            tCMStrucLinIso2->set_property( tPropNu2, "PoissonRatio" );
            tCMStrucLinIso2->set_model_type( fem::Model_Type::FULL );
            tCMStrucLinIso2->set_space_dim( 3 );
            tCMStrucLinIso2->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
            tSPNitscheInterface->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );
            tSPNitscheInterface->set_property( tPropEMod2, "Material", mtk::Leader_Follower::FOLLOWER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk1->set_residual_dof_type( { tResDofTypes } );
            tIWGBulk1->set_dof_type_list( { tResDofTypes } );
            tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk2->set_residual_dof_type( { tResDofTypes } );
            tIWGBulk2->set_dof_type_list( { tResDofTypes } );
            tIWGBulk2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { tResDofTypes } );
            tIWGDirichlet->set_dof_type_list( { tResDofTypes } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { tResDofTypes } );
            tIWGNeumann->set_dof_type_list( { tResDofTypes } );
            tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGInterface =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
            tIWGInterface->set_residual_dof_type( { tResDofTypes } );
            tIWGInterface->set_dof_type_list( { tResDofTypes } );
            tIWGInterface->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::FOLLOWER );
            tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
            tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

            // create the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUX->set_quantity_dof_type( tResDofTypes );
            tIQIUX->set_dof_type_list( { { tResDofTypes } }, mtk::Leader_Follower::LEADER );
            tIQIUX->set_output_type_index( 0 );
            tIQIUX->set_name( "IQI_UX" );

            std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUY->set_quantity_dof_type( tResDofTypes );
            tIQIUY->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
            tIQIUY->set_output_type_index( 1 );
            tIQIUY->set_name( "IQI_UY" );

            std::shared_ptr< fem::IQI > tIQIUZ = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUZ->set_quantity_dof_type( tResDofTypes );
            tIQIUZ->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
            tIQIUZ->set_output_type_index( 2 );
            tIQIUZ->set_name( "IQI_UZ" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );
            tSetBulk1.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk1 } );
            tSetBulk2.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetBulk3;
            tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk3.set_IWGs( { tIWGBulk2 } );
            tSetBulk3.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetBulk4;
            tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk4.set_IWGs( { tIWGBulk2 } );
            tSetBulk4.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p0" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            fem::Set_User_Info tSetInterface;
            tSetInterface.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 ) );
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
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_FEM_Benchmark_Elast_Interface_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                    { "Displacement UX", "Displacement UY", "Displacement UZ" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_UX", "IQI_UY", "IQI_UZ" } );

            tModel->set_output_manager( &tOutputData );

            // --------------------------------------------------------------------------------------
            // define linear solver and algorithm
            dla::Solver_Factory                             tSolFactory;
            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_amesos();
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // --------------------------------------------------------------------------------------
            // define nonlinear solver and algorithm
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // --------------------------------------------------------------------------------------
            // define time solver and algorithm
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { tResDofTypes } );
            tTimeSolver.set_dof_type_list( { tResDofTypes } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench );

            // --------------------------------------------------------------------------------------
            // solve and check
            tTimeSolver.solve();
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );

            delete tInterpMesh;
        }

    } /* END_TEST_CASE */

    TEST_CASE( "MDL FEM Benchmark Elast Ghost", "[MDL_FEM_Benchmark_Elast_Ghost]" )
    {
        if ( par_size() == 1 )
        {
            //	std::cout<<"I am proc: "<<par_rank()<<std::endl;

            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 2 }, { 2 } } );
            tParameters.set_domain_dimensions( 10, 5, 5 );
            tParameters.set_domain_offset( 0.0, 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );
            //        tParameters.set_lagrange_input_mesh( { { 0 } } );

            tParameters.set_staircase_buffer( 1 );

            tParameters.set_initial_refinement( { { 0 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );

            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );    // HMR Lagrange mesh
            //==============================
            std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex );

            tField->evaluate_scalar_function( tPlane_MDLFEMBench );

            Vector< std::shared_ptr< moris::hmr::Field > > tFields( 1, tField );

            // FIXME what is the following test about
            if ( false )
            {
                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tField );
                    //            tHMR.user_defined_flagging( user_defined_refinement_MDLFEMBench, tFields, tParam, 0 );
                    tHMR.perform_refinement_based_on_working_pattern( 0, true );
                    tField->evaluate_scalar_function( tPlane_MDLFEMBench );
                }
            }

            tHMR.finalize();

            //==============================
            tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tPlane = std::make_shared< moris::gen::Line >( 2.6, 0.0, 1.0, 0.0 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

            size_t tModelDimension = 3;

            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
            tXTKModel.construct_face_oriented_ghost_penalization_cells();

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            moris_index tSSIndex = tEnrIntegMesh.create_side_set_from_dbl_side_set( 1, "ghost_ss_p0" );
            tEnrIntegMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_p0", mtk::CellTopology::QUAD4 );

            //==============================

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
            tPropEMod1->set_parameters( { { { 1.0 } } } );
            tPropEMod1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
            tPropEMod2->set_parameters( { { { 1.0 } } } );
            tPropEMod2->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropNu1 = std::make_shared< fem::Property >();
            tPropNu1->set_parameters( { { { 0.3 } } } );
            tPropNu1->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropNu2 = std::make_shared< fem::Property >();
            tPropNu2->set_parameters( { { { 0.3 } } } );
            tPropNu2->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 }, { 0.0 } } } );
            tPropDirichlet->set_val_function( tPropValConstFunc_MDLFEMBench );

            std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
            tPropTraction->set_parameters( { { { 1.0 }, { 0.0 }, { 0.0 } } } );
            tPropTraction->set_val_function( tPropValConstFunc_MDLFEMBench );

            // working do types
            Vector< moris::MSI::Dof_Type > tResDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ };

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { tResDofTypes } );
            tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropNu1, "PoissonRatio" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::FULL );
            tCMStrucLinIso1->set_space_dim( 3 );
            tCMStrucLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso2->set_dof_type_list( { tResDofTypes } );
            tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
            tCMStrucLinIso2->set_property( tPropNu2, "PoissonRatio" );
            tCMStrucLinIso2->set_model_type( fem::Model_Type::FULL );
            tCMStrucLinIso2->set_space_dim( 3 );
            tCMStrucLinIso2->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
            tSPNitscheInterface->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );
            tSPNitscheInterface->set_property( tPropEMod2, "Material", mtk::Leader_Follower::FOLLOWER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPGhost =
                    tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
            tSPGhost->set_parameters( { { { 0.01 } } } );
            tSPGhost->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk1->set_residual_dof_type( { tResDofTypes } );
            tIWGBulk1->set_dof_type_list( { tResDofTypes } );
            tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulk2 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk2->set_residual_dof_type( { tResDofTypes } );
            tIWGBulk2->set_dof_type_list( { tResDofTypes } );
            tIWGBulk2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { tResDofTypes } );
            tIWGDirichlet->set_dof_type_list( { tResDofTypes } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { tResDofTypes } );
            tIWGNeumann->set_dof_type_list( { tResDofTypes } );
            tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGInterface =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
            tIWGInterface->set_residual_dof_type( { tResDofTypes } );
            tIWGInterface->set_dof_type_list( { tResDofTypes } );
            tIWGInterface->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::FOLLOWER );
            tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
            tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

            // Ghost stabilization
            std::shared_ptr< fem::IWG > tIWGGhost =
                    tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
            tIWGGhost->set_residual_dof_type( { tResDofTypes } );
            tIWGGhost->set_dof_type_list( { tResDofTypes } );
            tIWGGhost->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::FOLLOWER );
            tIWGGhost->set_stabilization_parameter( tSPGhost, "GhostSP" );

            // create the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUX->set_quantity_dof_type( tResDofTypes );
            tIQIUX->set_dof_type_list( { { tResDofTypes } }, mtk::Leader_Follower::LEADER );
            tIQIUX->set_output_type_index( 0 );
            tIQIUX->set_name( "IQI_UX" );

            std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUY->set_quantity_dof_type( tResDofTypes );
            tIQIUY->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
            tIQIUY->set_output_type_index( 1 );
            tIQIUY->set_name( "IQI_UY" );

            std::shared_ptr< fem::IQI > tIQIUZ = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUZ->set_quantity_dof_type( tResDofTypes );
            tIQIUZ->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
            tIQIUZ->set_output_type_index( 2 );
            tIQIUZ->set_name( "IQI_UZ" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );
            tSetBulk1.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk1 } );
            tSetBulk2.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetBulk3;
            tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk3.set_IWGs( { tIWGBulk2 } );
            tSetBulk3.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetBulk4;
            tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk4.set_IWGs( { tIWGBulk2 } );
            tSetBulk4.set_IQIs( { tIQIUX, tIQIUY, tIQIUZ } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p0" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            fem::Set_User_Info tSetInterface;
            tSetInterface.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 ) );
            tSetInterface.set_IWGs( { tIWGInterface } );

            fem::Set_User_Info tSetGhost;
            tSetGhost.set_mesh_set_name( "ghost_p0" );
            tSetGhost.set_IWGs( { tIWGGhost } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 8 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetBulk3;
            tSetInfo( 3 ) = tSetBulk4;
            tSetInfo( 4 ) = tSetDirichlet;
            tSetInfo( 5 ) = tSetNeumann;
            tSetInfo( 6 ) = tSetInterface;
            tSetInfo( 7 ) = tSetGhost;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_FEM_Benchmark_Elast_Ghost_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                    { "Displacement UX", "Displacement UY", "Displacement UZ" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_UX", "IQI_UY", "IQI_UZ" } );

            tModel->set_output_manager( &tOutputData );

            // --------------------------------------------------------------------------------------
            // define linear solver and algorithm
            dla::Solver_Factory                             tSolFactory;
            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_amesos();
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // --------------------------------------------------------------------------------------
            // define nonlinear solver and algorithm
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // --------------------------------------------------------------------------------------
            // define time solver and algorithm
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { tResDofTypes } );
            tTimeSolver.set_dof_type_list( { tResDofTypes } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench );

            // --------------------------------------------------------------------------------------
            // solve and check
            tTimeSolver.solve();
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );

            delete tInterpMesh;
        }

    } /* END_TEST_CASE */

}    // namespace moris
