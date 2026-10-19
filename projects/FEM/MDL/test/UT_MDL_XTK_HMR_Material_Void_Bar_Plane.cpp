/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_Material_Void_Bar_Plane.cpp
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
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_Matrix.hpp"    //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"    // ALG/src

#include "cl_FEM_NodeProxy.hpp"          //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"       //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"    //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"      //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src

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
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

#include "cl_GEN_Line.hpp"

moris::real
PlaneMatVoidMDL( const moris::Matrix< moris::DDRMat >& aPoint )
{
    moris::real mXC = 0.9545459;
    moris::real mYC = 0.11;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return ( mNx * ( aPoint( 0 ) - mXC ) + mNy * ( aPoint( 1 ) - mYC ) );
}

void tConstValFunctionMatVoidMDL( moris::Matrix< moris::DDRMat >& aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >&                 aParameters,
        moris::fem::Field_Interpolator_Manager*                   aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

TEST_CASE( "XTK HMR Material Void Bar Intersected By Plane", "[XTK_HMR_PLANE_BAR_MAT_VOID_2D]" )
{
    if ( par_size() == 1 )
    {
        std::string tFieldName = "Geometry";

        moris::uint tLagrangeMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( 11, 4 );
        tParameters.set_domain_dimensions( 6, 2 );
        tParameters.set_domain_offset( -3, -1 );
        tParameters.set_bspline_truncation( true );

        tParameters.set_output_meshes( { { 0 } } );

        tParameters.set_lagrange_orders( { 1 } );
        tParameters.set_lagrange_patterns( { 0 } );

        tParameters.set_bspline_orders( { 1 } );
        tParameters.set_bspline_patterns( { 0 } );

        tParameters.set_create_side_sets( true );
        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tPlaneField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tPlaneField->evaluate_scalar_function( PlaneMatVoidMDL );

        for ( uint k = 0; k < 2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tPlaneField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
            tPlaneField->evaluate_scalar_function( PlaneMatVoidMDL );
        }

        tPlaneField->evaluate_scalar_function( PlaneMatVoidMDL );
        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_mat_void_ip.e" );

        hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        auto tPlane = std::make_shared< moris::gen::Line >( 0.9545459, 0.11, 1.0, 0.0 );

        Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

        size_t                                 tModelDimension = 2;
        moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometryVector;
        moris::gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

        moris::xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
        tXTKModel.mVerbose = false;

        // Specify decomposition Method and Cut Mesh ---------------------------------------
        Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
        tXTKModel.decompose( tDecompositionMethods );

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        moris::xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        moris::xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // tEnrIntegMesh.print();

        moris_index tSSIndex = tEnrIntegMesh.create_side_set_from_dbl_side_set( 1, "ghost_ss_p0" );
        tEnrIntegMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_p0", moris::mtk::CellTopology::QUAD4 );

        // Write mesh
        moris::mtk::Writer_Exodus writer( &tEnrIntegMesh );

        writer.write_mesh( "./mdl_exo", "xtk_hmr_bar_plane_mat_void_integ_2d_ghost.e", "./", "temp.exo" );

        // Write the fields
        writer.set_time( 0.0 );
        writer.close_file();

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
        tPropConductivity1->set_parameters( { { { 1.0 } } } );
        tPropConductivity1->set_val_function( tConstValFunctionMatVoidMDL );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { { { 5.0 } } } );
        tPropDirichlet->set_val_function( tConstValFunctionMatVoidMDL );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { { { 20.0 } } } );
        tPropNeumann->set_val_function( tConstValFunctionMatVoidMDL );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { { { 100.0 } } } );
        tPropTempLoad1->set_val_function( tConstValFunctionMatVoidMDL );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
        tCMDiffLinIso1->set_space_dim( 2 );
        tCMDiffLinIso1->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory                                 tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
        tSPDirichletNitsche->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPGhost = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
        tSPGhost->set_parameters( { { { 0.1 } } } );
        tSPGhost->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );

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
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGGhost = tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
        tIWGGhost->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGGhost->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGGhost->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGGhost->set_stabilization_parameter( tSPGhost, "GhostSP" );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
        tSetBulk1.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
        tSetBulk2.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p0" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        std::string        tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
        tSetNeumann.set_mesh_set_name( tInterfaceSideSetName );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetGhost;
        tSetGhost.set_mesh_set_name( "ghost_p0" );
        tSetGhost.set_IWGs( { tIWGGhost } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 5 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetDirichlet;
        tSetInfo( 3 ) = tSetNeumann;
        tSetInfo( 4 ) = tSetGhost;

        // create model
        mdl::Model* tModel = new mdl::Model( tMeshManager,
                0,
                tSetInfo,
                0,
                false );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory tSolFactory;
        Parameter_List      tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_all );
        tLinearSolverParameterList.set( "AZ_solver", AZ_gmres_condnum );
        tLinearSolverParameterList.set( "AZ_precond", AZ_none );
        tLinearSolverParameterList.set( "AZ_kspace", 500 );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;

        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();

        tNonlinearSolverParameterList.set( "NLA_max_iter", 3 );

        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm =
                tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

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

        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
        tTimeSolver.set_dof_type_list( tDofTypes1 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        //
        //        // verify solution
        //        moris::Matrix<DDRMat> tGoldSolution =
        //        {{ 5.000000e+00},
        //         { 2.500000e+01},
        //         { 4.500000e+01},
        //         { 6.500000e+01},
        //         { 5.000000e+00},
        //         { 2.500000e+01},
        //         { 4.500000e+01},
        //         { 6.500000e+01}};
        //
        //        Matrix<DDRMat> tFullSol;
        //        tTimeSolver.get_full_solution(tFullSol);
        //
        //        CHECK(norm(tFullSol - tGoldSolution)<1e-08);

        //    delete tInterpMesh1;
        delete tModel;
        delete tInterpMesh;
    }
}
