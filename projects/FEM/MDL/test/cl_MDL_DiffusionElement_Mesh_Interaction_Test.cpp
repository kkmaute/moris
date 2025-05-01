/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MDL_DiffusionElement_Mesh_Interaction_Test.cpp
 *
 */

#include "catch.hpp"

#include "paths.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"

#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

//#include "cl_VIS_Factory.hpp"
//#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

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

#include "fn_norm.hpp"

moris::real LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 0.5;
}

namespace moris::mdl
{

void tConstValFunction_MDLDIFF
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

TEST_CASE( "Diffusion_2x2x2", "[moris],[mdl],[Diffusion_2x2x2]" )
{

    if(par_size() == 1 )
    {
        // Create a 3D mesh of HEX8 using MTK ------------------------------------------
        std::string tPrefix = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "projects/FEM/INT/test/data/Cube_with_side_sets.g";

        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // Declare some supplementary fields
        mtk::MtkMeshData tMeshData;
        tMeshData.FieldsInfo = &tFieldsInfo;

        // construct the mesh data
        mtk::Interpolation_Mesh* tInterpMesh1 = mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tMeshData );
        mtk::Integration_Mesh*   tIntegMesh1  = mtk::create_integration_mesh_from_interpolation_mesh( mtk::MeshType::STK, tInterpMesh1 );

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(tInterpMesh1,tIntegMesh1);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
        tPropTempLoad->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad->set_val_function( tConstValFunction_MDLDIFF );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 3 );
        tCMDiffLinIso->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        // define the IQIs
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIQITEMP->set_output_type_index( 0 );

        // define set info
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_index( 0 ); // index 0
        tSetBulk.set_IWGs( { tIWGBulk } );
        tSetBulk.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( 4 ); // index 4
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( 6 ); //index 6
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 3 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetDirichlet;
        tSetInfo( 2 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               1,
                                               tSetInfo );

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 10 );
        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        moris::Matrix< DDRMat > tSolution11;
        tTimeSolver.get_full_solution( tSolution11 );

        // Expected solution
        Matrix< DDRMat > tExpectedSolution = {{ 25.0, 25.0, 25.0,
                                                25.0,  5.0, 25.0,
                                                45.0, 25.0,  5.0,
                                                25.0, 45.0, 25.0,
                                                5.0, 25.0, 45.0,
                                                5.0, 45.0,  5.0,
                                                45.0,  5.0, 45.0,
                                                5.0, 45.0,  5.0,
                                                45.0,  5.0, 45.0 }};

        // define an epsilon environment
        real tEpsilon = 1E-6;

        // define a bool for solution check
        bool tCheckNodalSolution = true;

        // number of mesh nodes
        uint tNumOfNodes = tInterpMesh1->get_num_nodes();

        // loop over the node and check solution
        for ( uint i = 0; i < tNumOfNodes; i++ )
        {
            // check solution
            tCheckNodalSolution = tCheckNodalSolution
                    && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
        }
        // check bool is true
        REQUIRE( tCheckNodalSolution );

        // clean up
        delete tModel;
        delete tIntegMesh1;
        delete tInterpMesh1;

    }/* if( par_size() */
}

TEST_CASE( "Element_Diffusion_3", "[moris],[mdl],[Diffusion_block_7x8x9]" )
{

    if(par_size() == 1 )
    {
        std::string tPrefix = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "projects/FEM/MDL/test/data/Block_7x8x9.g";

        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // Declare some supplementary fields
        mtk::MtkMeshData tMeshData;
        tMeshData.FieldsInfo = &tFieldsInfo;

        // construct the mesh data
        mtk::Interpolation_Mesh* tInterpMesh1 = mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tMeshData );
        mtk::Integration_Mesh*   tIntegMesh1  = mtk::create_integration_mesh_from_interpolation_mesh(mtk::MeshType::STK,tInterpMesh1);

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(tInterpMesh1,tIntegMesh1);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
        tPropTempLoad->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad->set_val_function( tConstValFunction_MDLDIFF );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 3 );
        tCMDiffLinIso->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        // define set info
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_index( 0 );
        tSetBulk.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( 4 );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( 6 );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 3 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetDirichlet;
        tSetInfo( 2 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               1,
                                               tSetInfo );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 10 );
        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        moris::Matrix< DDRMat > tSolution11;
        tTimeSolver.get_full_solution( tSolution11 );

        // Expected solution
        Matrix< DDRMat > tExpectedSolution( 1, 25, 2.5e+01 );

        // define an epsilon environment
        real tEpsilon = 1E-6;

        // define a bool for solution check
        bool tCheckNodalSolution = true;

        // loop over the node and check solution
        for ( uint i = 0; i < 25; i++ )
        {
            // check solution
            tCheckNodalSolution = tCheckNodalSolution
                    && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
        }
        // check bool is true
        REQUIRE( tCheckNodalSolution );

        // clean up
        delete tIntegMesh1;
        delete tInterpMesh1;
        delete tModel;

    }/* if( par_size() */
}

//-------------------------------------------------------------------------------------------------------
TEST_CASE( "Diffusion_hmr_10x4x4", "[moris],[mdl],[Diffusion_hmr_10x4x4]" )
{

    if( par_size() == 2 )
    {
        //------------------------------------------------------------------------------

        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( 10, 4, 4 );
        tParameters.set_domain_dimensions( 10, 4, 4 );
        tParameters.set_domain_offset( -10, -2, -2 );
        tParameters.set_bspline_truncation( true );
        tParameters.set_create_side_sets( true );

        tParameters.set_output_meshes( {{ 0 }} );

        tParameters.set_lagrange_orders  ( { 1 });
        tParameters.set_lagrange_patterns({ 0 });

        tParameters.set_bspline_orders   ( { 1 } );
        tParameters.set_bspline_patterns ( { 0 } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

//        tParameters.set_number_aura( true );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeMeshIndex );

        for( uint k=0; k<3; ++k )
        {
            tField->evaluate_scalar_function( LevelSetFunction );
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
        }

        tHMR.finalize();

        // construct a mesh manager for the fem
        moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);
        moris::hmr::Integration_Mesh_HMR *   tIntegrationMesh   = tHMR.create_integration_mesh(1, 0, tInterpolationMesh);

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( tInterpolationMesh, tIntegrationMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property > ();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property > ();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property > ();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property > ();
        tPropTempLoad->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad->set_val_function( tConstValFunction_MDLDIFF );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 3 );
        tCMDiffLinIso->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        // define set info
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_index( 0 );
        tSetBulk.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( 4 );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( 6 );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 3 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetDirichlet;
        tSetInfo( 2 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               tBSplineMeshIndex,
                                               tSetInfo );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 10 );
        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        moris::Matrix< DDRMat > tSolution11;
        tTimeSolver.get_full_solution( tSolution11 );

        // Expected solution
        Matrix< DDRMat > tExpectedSolution;

        if ( par_rank() == 0 )
        {
            // Expected solution for first processor
            tExpectedSolution = {{ 
                    +5.190827153138305e+00, +4.682782578636699e+01, +8.632909250554140e+01, +1.224604658451722e+02, +1.543627066812462e+02, 
                    +1.816024902272104e+02, +5.190827153402125e+00, +4.682782578458072e+01, +8.632909249584647e+01, +1.224604657663189e+02, 
                    +1.543627060330361e+02, +1.816024855001556e+02, +5.190827153683495e+00, +4.682782578206277e+01, +8.632909248106306e+01, 
                    +1.224604656778186e+02, +1.543627053930827e+02, +1.816024806081994e+02, +5.190827153559872e+00, +4.682782578185444e+01, 
                    +8.632909248634425e+01, +1.224604657408669e+02, +1.543627059794248e+02, +1.816024853324490e+02, +5.190827153434689e+00 }};
        }
        else if ( par_rank() == 1 )
        {
            // Expected solution for second processor
            tExpectedSolution = {{
                    +1.816024902272104e+02, +2.039890594364211e+02, +2.214404864961653e+02, +2.339232490921326e+02, +2.386405535437454e+02,
                    +2.386406610739002e+02, +2.392004587257753e+02, +2.392002833591539e+02, +2.417340664637686e+02, +2.436172967529151e+02,
                    +2.417352168660505e+02, +2.436198791258595e+02, +2.423633165461652e+02, +2.442398541418869e+02, +2.423638901847680e+02,
                    +2.442436149253762e+02, +2.442441375582408e+02, +2.442473271171616e+02, +2.448681914413924e+02, +2.448732020403027e+02,
                    +1.816024855001556e+02, +2.039890202246969e+02, +2.214401928865578e+02, +2.339230300671687e+02, +2.386409485619482e+02 }};
            print(tSolution11,"Processor_TWO");
        }

        // define an epsilon environment
        double tEpsilon = 1E-6;

        // define a bool for solution check
        bool tCheckNodalSolution = true;

        // loop over the node and chyeck solution
        for ( uint i = 0; i < 25; i++ )
        {
            // check solution
            tCheckNodalSolution = tCheckNodalSolution
                    && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
        }

        // check bool is true
        REQUIRE( tCheckNodalSolution );

        // clean up
        delete tModel;
        delete tInterpolationMesh;
        delete tIntegrationMesh;

    }/* if( par_size() */
}

TEST_CASE( "Diffusion_hmr3_10x4x4", "[moris],[mdl],[Diffusion_hmr3_10x4x4]" )
{

    if( (par_size() == 1) || (par_size() == 2) )
    {
        //------------------------------------------------------------------------------

        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( 4, 2, 2 );
        tParameters.set_domain_dimensions( 4, 2, 2 );
        tParameters.set_domain_offset( -2, 0, 0 );
        tParameters.set_bspline_truncation( true );
        tParameters.set_create_side_sets( true );

        tParameters.set_output_meshes( {{ 0 }} );

        tParameters.set_lagrange_orders  ( { 2 });
        tParameters.set_lagrange_patterns({ 0 });

        tParameters.set_bspline_orders   ( { 2 } );
        tParameters.set_bspline_patterns ( { 0 } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 1 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeMeshIndex );

        for( uint k=0; k<3; ++k )
        {
            tField->evaluate_scalar_function( LevelSetFunction );
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
        }

        tHMR.finalize();

        // construct a mesh manager for the fem
        moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        moris::hmr::Integration_Mesh_HMR *   tIntegrationMesh   = tHMR.create_integration_mesh( 2, 0, tInterpolationMesh );

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( tInterpolationMesh,tIntegrationMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property > ();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property > ();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property > ();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property > ();
        tPropTempLoad->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad->set_val_function( tConstValFunction_MDLDIFF );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 3 );
        tCMDiffLinIso->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        // define set info
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_index( 0 );
        tSetBulk.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( 4 );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( 6 );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 3 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetDirichlet;
        tSetInfo( 2 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               tBSplineMeshIndex,
                                               tSetInfo );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 10 );
        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        moris::Matrix< DDRMat > tSolution11;
        tTimeSolver.get_full_solution( tSolution11 );

        // Expected solution when running in serial
        Matrix< DDRMat > tExpectedSolution = {{ 
                -2.998789793682302e+00, +1.282912197972422e+01, -2.990108456510009e+00, +1.282479340728375e+01, -2.727515379981700e+00,
                +1.276753440802280e+01, -2.726961961458098e+00, +1.276722880234914e+01, +2.796549355435167e+01, +3.146278456764225e+01,
                +3.837636418457888e+01, +3.146178009510643e+01, +3.837415552376613e+01, +3.146216179259339e+01, +3.837359215288387e+01,
                +3.146246688695902e+01, +3.837433884926831e+01, +4.483623755822933e+01, +4.633269836608036e+01, +4.931554740561488e+01,
                +4.633246928240963e+01, +4.931568705745656e+01, +4.633241987224589e+01, +4.931564286552350e+01, +4.633285989583719e+01 }};

        // expected solutions when running in parallel
        if (par_size() == 2)
        {
            if ( par_rank() == 0 )
            {
                // Expected solution for first processor
                tExpectedSolution = {{ 
                        -2.998789935777626e+00, +1.282912205236343e+01, -2.990108399856919e+00, +1.282479338003131e+01, -2.727515342468292e+00, 
                        +1.276753439020575e+01, -2.726962005061556e+00, +1.276722882490287e+01, +2.796549348107517e+01, +3.146278462203448e+01, 
                        +3.837636420761202e+01, +3.146178012302775e+01, +3.837415558256843e+01, +3.146216182742653e+01, +3.837359220480018e+01, 
                        +3.146246692320263e+01, +3.837433889477731e+01, +4.483623765560167e+01, +4.633269842004368e+01, +4.931554750475398e+01, 
                        +4.633246935266567e+01, +4.931568712228110e+01, +4.633241994480259e+01, +4.931564293169782e+01, +4.633285995928019e+01 }};
            }
            else if ( par_rank() == 1 )
            {
                // Expected solution for second processor
                tExpectedSolution = {{
                        +5.623957216613431e+01, +6.111112047899337e+01, +5.624098542310917e+01, +6.110966639970183e+01, +5.624058907310716e+01,
                        +6.110952500726074e+01, +5.624071044186167e+01, +6.110967527872084e+01, +6.539470546310521e+01, +6.631870453975621e+01,
                        +6.816339220548488e+01, +6.631861896236025e+01, +6.816340876099639e+01, +6.631861378394017e+01, +6.816339510715262e+01,
                        +6.631868535188821e+01, +6.816333175900422e+01, +6.908659478022295e+01, +6.539600520925121e+01, +6.908611715265627e+01,
                        +6.539590999060719e+01, +6.908607267216748e+01, +6.539591583595094e+01, +6.908603017351150e+01, +5.624040154355884e+01 }};
            }
            else {} //do nothing
        }
        else {} // end expected solutions for parallel

        // define an epsilon environment
        real tEpsilon = 1E-6;

        // define a bool for solution check
        bool tCheckNodalSolution = true;

        // loop over the node and check solution
        for ( uint i = 0; i < 25; i++ )
        {
            // check solution
            tCheckNodalSolution = tCheckNodalSolution
                    && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
        }
        // check bool is true
        REQUIRE( tCheckNodalSolution );

        delete tIntegrationMesh;
        delete tInterpolationMesh;

    }/* if( par_size() */
}

//-------------------------------------------------------------------------------------------------------

TEST_CASE( "Diffusion_hmr_cubic_10x4x4", "[moris],[mdl],[Diffusion_hmr_cubic_10x4x4]" )
{

    if( par_size() == 1 )
    {
        //------------------------------------------------------------------------------

        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( 4, 2, 2 );
        tParameters.set_domain_dimensions( 4, 2, 2 );
        tParameters.set_domain_offset( -2, 0, 0 );
        tParameters.set_bspline_truncation( true );
        tParameters.set_create_side_sets( true );

        tParameters.set_output_meshes( {{ 0 }} );

        tParameters.set_lagrange_orders  ( { 3 });
        tParameters.set_lagrange_patterns({ 0 });

        tParameters.set_bspline_orders   ( { 3 } );
        tParameters.set_bspline_patterns ( { 0 } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeMeshIndex );

        for( uint k=0; k<2; ++k )
        {
            tField->evaluate_scalar_function( LevelSetFunction );
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
        }

        tHMR.finalize();

        // construct a mesh manager for the fem
        moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        moris::hmr::Integration_Mesh_HMR *   tIntegrationMesh   = tHMR.create_integration_mesh( 3, 0, tInterpolationMesh );

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(tInterpolationMesh,tIntegrationMesh);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property > ();
        tPropConductivity->set_parameters( { {{ 1.0 }} } );
        tPropConductivity->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property > ();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property > ();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction_MDLDIFF );

        std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property > ();
        tPropTempLoad->set_parameters( { {{ 0.0 }} } );
        tPropTempLoad->set_val_function( tConstValFunction_MDLDIFF );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_space_dim( 3 );
        tCMDiffLinIso->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Leader_Follower::LEADER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche , "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        // define set info
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_index( 0 );
        tSetBulk.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( 4 );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( 6 );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 3 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetDirichlet;
        tSetInfo( 2 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               tBSplineMeshIndex,
                                               tSetInfo );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algortihm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 10 );
        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        moris::Matrix< DDRMat > tSolution11;
        tTimeSolver.get_full_solution( tSolution11 );

        // Expected solution when running in serial
        Matrix< DDRMat > tExpectedSolution = {{ 
                -6.098822088327460e+00, +4.585718808771261e+00, -6.644452157658328e+00, +4.711372418219552e+00, -1.137258126847988e+01,
                +5.032341577176382e+00, -1.136540426268445e+01, +5.030584956857357e+00, +2.103739166755025e+01, +3.559502876141929e+01,
                +4.807817199377635e+01, +4.808686482256833e+01, +4.783808188916422e+01, +4.783238249532101e+01, +2.093931702174379e+01,
                +3.574069062951020e+01, +2.043249172806436e+01, +3.493808144516003e+01, +2.043375158956196e+01, +3.493632365610171e+01,
                +4.926338936207750e+01, +5.383366943645356e+01, +5.895920244464836e+01, +5.380440335245414e+01, +5.896536165925041e+01 }};

        // define an epsilon environment
        real tEpsilon = 1E-6;

        // define a bool for solution check
        bool tCheckNodalSolution = true;

        // loop over the node and check solution of the first 25 values
        for ( uint i = 0; i < 25; i++ )
        {
            // check solution
            tCheckNodalSolution = tCheckNodalSolution
                    && ( std::abs( tSolution11( i ) - tExpectedSolution( i ) ) < tEpsilon );
        }
        // check bool is true
        REQUIRE( tCheckNodalSolution );

        // clean up
        delete tModel;
        delete tIntegrationMesh;
        delete tInterpolationMesh;

    }/* if( par_size() */
}

}    // namespace moris::mdl
