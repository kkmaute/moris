/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_Linear_Struc_2D.cpp
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

#include "cl_FEM_NodeProxy.hpp"          //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"       //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"    //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"         //FEM/INT/src
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

#include "cl_GEN_Circle.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris
{
    inline moris::real
    LvlSetLin( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real tOffset = 200;
        return aPoint( 0 ) - 0.317 * aPoint( 1 ) - tOffset;
    }

    inline moris::real
    LvlSetCircle_2D( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        return std::sqrt( aPoint( 0 ) * aPoint( 0 ) + aPoint( 1 ) * aPoint( 1 ) ) - 0.2505;
    }

    inline moris::real
    LvlSetCircle_2D_outsideDomain( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        Matrix< DDRMat > tCenter{ { -10, -10 } };
        return norm( aPoint - tCenter ) - 0.001;
    }

    inline void
    tConstValFunction( moris::Matrix< moris::DDRMat >&     aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    inline void
    tMValFunction( moris::Matrix< moris::DDRMat >&         aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { aParameters( 0 )( 0 ), 0.0 },
            { 0.0, aParameters( 0 )( 1 ) } };
    }

    inline void
    tMValFunction_3D( moris::Matrix< moris::DDRMat >&      aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { aParameters( 0 )( 0 ), 0.0, 0.0 },
            { 0.0, aParameters( 0 )( 1 ), 0.0 },
            { 0.0, 0.0, aParameters( 0 )( 2 ) } };
    }

    inline moris::real
    LevelSetFunction_star1( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real tPhi          = std::atan2( aPoint( 0 ), aPoint( 1 ) );
        moris::real tLevelSetVaue = 0.501 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aPoint( 0 ), 2 ) + std::pow( aPoint( 1 ), 2 ) );
        return -tLevelSetVaue;
    }

    inline moris::real
    Plane4MatMDL1( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real mXC = 0.1;
        moris::real mYC = 0.1;
        moris::real mNx = 1.0;
        moris::real mNy = 0.0;
        return ( mNx * ( aPoint( 0 ) - mXC ) + mNy * ( aPoint( 1 ) - mYC ) );
    }

    inline moris::real
    Circle4MatMDL( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real mXCenter = 0.01;
        moris::real mYCenter = 0.01;
        moris::real mRadius  = 0.47334;

        return ( aPoint( 0 ) - mXCenter ) * ( aPoint( 0 ) - mXCenter )
             + ( aPoint( 1 ) - mYCenter ) * ( aPoint( 1 ) - mYCenter )
             - ( mRadius * mRadius );
    }

    TEST_CASE( "2D XTK WITH HMR Struc Interface 2D", "[XTK_HMR_Struc_Interface_2D]" )
    {
        if ( par_size() <= 1 )
        {
            uint        tLagrangeMeshIndex = 0;
            std::string tFieldName         = "Cylinder";

            moris::ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "20, 20" );
            tParameters.set( "domain_dimensions", "2, 2" );
            tParameters.set( "domain_offset", "-1.0, -1.0" );
            tParameters.set( "domain_sidesets", "1,2,3,4" );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "1" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            //  create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LvlSetCircle_2D );

            for ( uint k = 0; k < 2; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LvlSetCircle_2D );
            }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tCircle = std::make_shared< moris::gen::Circle >( 0.0, 0.0, 0.2501 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };

            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( 2, tInterpolationMesh, &tGeometryEngine );

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

            // get meshes
            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
            tPropEMod1->set_parameters( { { { 1.0 } } } );
            tPropEMod1->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
            tPropEMod2->set_parameters( { { { 1.0 } } } );
            tPropEMod2->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
            tPropNu->set_parameters( { { { 0.0 } } } );
            tPropNu->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet2 = std::make_shared< fem::Property >();
            tPropDirichlet2->set_parameters( { { { 1.0, 1.0 } } } );
            tPropDirichlet2->set_val_function( tMValFunction );

            std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
            tPropTraction->set_parameters( { { { 1.0 }, { 0.0 } } } );
            tPropTraction->set_val_function( tConstValFunction );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropNu, "PoissonRatio" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso1->set_space_dim( 2 );
            tCMStrucLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
            tCMStrucLinIso2->set_property( tPropNu, "PoissonRatio" );
            tCMStrucLinIso2->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso2->set_space_dim( 2 );
            tCMStrucLinIso2->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitsche->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
            tSPNitscheInterface->set_property( tPropEMod2, "Material", mtk::Leader_Follower::LEADER );
            tSPNitscheInterface->set_property( tPropEMod1, "Material", mtk::Leader_Follower::FOLLOWER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulk2 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk2->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet2, "Select", mtk::Leader_Follower::LEADER );

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
            tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

            // create a list of active block-sets
            std::string tInterfaceSideSetName    = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
            std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk2 } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk2 } );

            fem::Set_User_Info tSetBulk3;
            tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk3.set_IWGs( { tIWGBulk1 } );

            fem::Set_User_Info tSetBulk4;
            tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk4.set_IWGs( { tIWGBulk1 } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

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
                    tSetInfo,
                    0,
                    false );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
            tDofTypesU( 0 ) = MSI::Dof_Type::UX;
            tDofTypesU( 1 ) = MSI::Dof_Type::UY;

            dla::Solver_Factory                             tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param( "AZ_diagnostics" )     = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_output" )          = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_max_iter" )        = 10000;
            tLinearSolverAlgorithm->set_param( "AZ_solver" )          = AZ_gmres;
            tLinearSolverAlgorithm->set_param( "AZ_subdomain_solve" ) = AZ_ilu;
            tLinearSolverAlgorithm->set_param( "AZ_graph_fill" )      = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param( "NLA_max_iter" ) = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

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

            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            // output solution and meshes
            // FIXME if needed

            // clean up
            delete tModel;
            delete tInterpolationMesh;
        }
    }

    TEST_CASE( "2D XTK WITH HMR Struc Interface 3D", "[XTK_HMR_Struc_Interface_3D]" )    // FIXME
    {
        //    if(par_size()<=1)
        //    {
        //        uint tLagrangeMeshIndex = 0;
        //        std::string tFieldName = "Cylinder";
        //
        //        moris::ParameterList tParameters = prm::create_hmr_parameter_list();
        //
        //        tParameters.set( "number_of_elements_per_dimension", "2, 2, 2");
        //        tParameters.set( "domain_dimensions", "2, 2, 2" );
        //        tParameters.set( "domain_offset", "-1.0, -1.0, -1.0" );
        //        tParameters.set( "domain_sidesets", "1,2,3,4,5, 6" );
        //        tParameters.set( "lagrange_output_meshes", "0" );
        //
        //        tParameters.set( "lagrange_orders", "1" );
        //        tParameters.set( "lagrange_pattern",std::string( "0") );
        //        tParameters.set( "bspline_orders", "1" );
        //        tParameters.set( "bspline_pattern", "0" );
        //
        //        tParameters.set( "lagrange_to_bspline", "0" );
        //
        //        tParameters.set( "truncate_bsplines", 1 );
        //        tParameters.set( "refinement_buffer", 3 );
        //        tParameters.set( "staircase_buffer", 3 );
        //        tParameters.set( "initial_refinement", 0 );
        //
        //        tParameters.set( "use_multigrid", 0 );
        //        tParameters.set( "severity_level", 2 );
        //
        //        hmr::HMR tHMR( tParameters );
        //
        //        //initial refinement
        //        tHMR.perform_initial_refinement( 0 );
        //
        //        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
        //
        //        //  create field
        //        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
        //
        //        tField->evaluate_scalar_function( LvlSetCircle_2D );
        //
        //        tHMR.finalize();
        //
        //        tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_3d.e" );
        //
        //        moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);
        //
        //        Vector<std::shared_ptr<moris::gen::Level_Set_Geometry>> tGeometry(1);
        //        tGeometry(0) = std::make_shared<moris::gen::Circle>(100.0, 0.0, 0.2501);
        //
        //        moris::gen::Phase_Table tPhaseTable (1);
        //        moris::gen::Geometry_Engine tGeometryEngine(tGeometry, tPhaseTable, 3);
        //
        //         xtk::Model tXTKModel(3, tInterpolationMesh, &tGeometryEngine);
        //
        //        //Specify decomposition Method and Cut Mesh ---------------------------------------
        //        Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        //        tXTKModel.decompose(tDecompositionMethods);
        //
        //        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);
        //
        //        // get meshes
        //        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        //        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
        //
        //        // place the pair in mesh manager
        //        mtk::Mesh_Manager tMeshManager;
        //        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
        //
        //        //------------------------------------------------------------------------------
        //        // create the properties
        //        std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
        //        tPropEMod1->set_parameters( { {{ 1.0 }} } );
        //        tPropEMod1->set_val_function( tConstValFunction );
        //
        //        std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
        //        tPropEMod2->set_parameters( { {{ 1.0 }} } );
        //        tPropEMod2->set_val_function( tConstValFunction );
        //
        //        std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
        //        tPropNu->set_parameters( { {{ 0.3 }} } );
        //        tPropNu->set_val_function( tConstValFunction );
        //
        //        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        //        tPropDirichlet->set_parameters( { {{ 0.0 }, { 0.0 }, { 0.0 }} } );
        //        tPropDirichlet->set_val_function( tConstValFunction );
        //
        //        std::shared_ptr< fem::Property > tPropDirichlet2 = std::make_shared< fem::Property >();
        //        tPropDirichlet2->set_parameters( { {{ 1.0, 1.0, 1.0 }} } );
        //        tPropDirichlet2->set_val_function( tMValFunction_3D );
        //
        //        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        //        tPropNeumann->set_parameters( {{{ 1.0 } , { 0.0 }, { 0.0 }}} );
        //        tPropNeumann->set_val_function( tConstValFunction );
        //
        //        // define constitutive models
        //        fem::CM_Factory tCMFactory;
        //
        //        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        //        tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        //        tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
        //        tCMStrucLinIso1->set_property( tPropNu, "PoissonRatio" );
        //        tCMStrucLinIso1->set_space_dim( 3 );
        //
        //        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        //        tCMStrucLinIso2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        //        tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
        //        tCMStrucLinIso2->set_property( tPropNu, "PoissonRatio" );
        //        tCMStrucLinIso2->set_space_dim( 3 );
        //
        //        // define stabilization parameters
        //        fem::SP_Factory tSPFactory;
        //        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        //        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        //        tSPDirichletNitsche->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );
        //
        //        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        //        tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
        //        tSPNitscheInterface->set_property( tPropEMod2, "Material", mtk::Leader_Follower::LEADER );
        //        tSPNitscheInterface->set_property( tPropEMod1, "Material", mtk::Leader_Follower::FOLLOWER );
        //
        //        // define the IWGs
        //        fem::IWG_Factory tIWGFactory;
        //
        //        std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        //        tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  } );
        //        tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  }} );
        //        tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
        //
        //        std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        //        tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  } );
        //        tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  }} );
        //        tIWGBulk2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );
        //
        //        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        //        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  } );
        //        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  }} );
        //        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        //        tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
        //        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
        //        tIWGDirichlet->set_property( tPropDirichlet2, "Select", mtk::Leader_Follower::LEADER );
        //
        //        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        //        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  } );
        //        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  }} );
        //        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );
        //
        //        std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
        //        tIWGInterface->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  } );
        //        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  }} );
        //        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ  }},mtk::Leader_Follower::FOLLOWER );
        //        tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
        //        tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );
        //        tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
        //
        //        // create a list of active block-sets
        //        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
        //        std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        //
        //        // define set info
        //        fem::Set_User_Info tSetBulk1;
        //        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
        //        tSetBulk1.set_IWGs( { tIWGBulk2 } );
        //
        //        fem::Set_User_Info tSetBulk2;
        //        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
        //        tSetBulk2.set_IWGs( { tIWGBulk2 } );
        //
        //        fem::Set_User_Info tSetBulk3;
        //        tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
        //        tSetBulk3.set_IWGs( { tIWGBulk1 } );
        //
        //        fem::Set_User_Info tSetBulk4;
        //        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        //        tSetBulk4.set_IWGs( { tIWGBulk1 } );
        //
        //        fem::Set_User_Info tSetDirichlet;
        //        tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
        //        tSetDirichlet.set_IWGs( { tIWGDirichlet } );
        //
        //        fem::Set_User_Info tSetNeumann;
        //        tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
        //        tSetNeumann.set_IWGs( { tIWGNeumann } );
        //
        //        fem::Set_User_Info tSetInterface;
        //        tSetInterface.set_mesh_set_name( tDblInterfaceSideSetName );
        //        tSetInterface.set_IWGs( { tIWGInterface } );
        //
        //        // create a cell of set info
        //        Vector< fem::Set_User_Info > tSetInfo( 7 );
        //        tSetInfo( 0 ) = tSetBulk1;
        //        tSetInfo( 1 ) = tSetBulk2;
        //        tSetInfo( 2 ) = tSetBulk3;
        //        tSetInfo( 3 ) = tSetBulk4;
        //        tSetInfo( 4 ) = tSetDirichlet;
        //        tSetInfo( 5 ) = tSetNeumann;
        //        tSetInfo( 6 ) = tSetInterface;
        //
        //        // create model
        //        mdl::Model * tModel = new mdl::Model( &tMeshManager,
        //                                              0,
        //                                              tSetInfo,
        //                                              0, false );
        //
        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //        // STEP 1: create linear solver and algorithm
        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //
        //        Vector< enum MSI::Dof_Type > tDofTypesU( 3 );
        //        tDofTypesU( 0 ) = MSI::Dof_Type::UX;
        //        tDofTypesU( 1 ) = MSI::Dof_Type::UY;
        //        tDofTypesU( 2 ) = MSI::Dof_Type::UZ;
        //
        //        dla::Solver_Factory  tSolFactory;
        //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
        //
        //        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        //        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
        //        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10000;
        //        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
        //        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
        //        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 10;
        //        //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";
        //
        //        dla::Linear_Solver tLinSolver;
        //        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
        //
        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //        // STEP 2: create nonlinear solver and algorithm
        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
        //        //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
        //
        //        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
        //        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
        //        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
        //        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;
        //
        //        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
        //        //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );
        //
        //        NLA::Nonlinear_Solver tNonlinearSolverMain;
        //        tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
        //
        //
        //        tNonlinearSolverMain       .set_dof_type_list( tDofTypesU );
        //
        //        // Create solver database
        //        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
        //
        //        tNonlinearSolverMain       .set_solver_warehouse( &tSolverWarehouse );
        //
        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //        // STEP 3: create time Solver and algorithm
        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //        tsa::Time_Solver_Factory tTimeSolverFactory;
        //        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
        //
        //        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );
        //
        //        tsa::Time_Solver tTimeSolver;
        //        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        //        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
        //
        //        tTimeSolver.set_dof_type_list( tDofTypesU );
        //
        //        //------------------------------------------------------------------------------
        //        tTimeSolver.solve();
        //
        //        // output solution and meshes
        //        // FIXME if needed
        //
        //        delete tModel;
        //        delete tInterpolationMesh;
        //    }
    }

    TEST_CASE( "2D XTK WITH HMR Struc 2D first", "[XTK_HMR_Struc_2D_01]" )
    {
        if ( par_size() <= 1 )
        {
            uint        tLagrangeMeshIndex = 0;
            std::string tFieldName         = "Cylinder";

            moris::ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "20, 20" );
            tParameters.set( "domain_dimensions", "2, 2" );
            tParameters.set( "domain_offset", "-1.0, -1.0" );
            tParameters.set( "domain_sidesets", "1,2,3,4" );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", std::string( "1" ) );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            //// create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LvlSetCircle_2D );
            //
            for ( uint k = 0; k < 2; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LvlSetCircle_2D );
            }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tCircle = std::make_shared< moris::gen::Circle >( 0.0, 0.0, 0.2505 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };

            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( 2, tInterpolationMesh, &tGeometryEngine );

            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            bool                            tSuccess              = tXTKModel.decompose( tDecompositionMethods );
            CHECK( tSuccess );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::NODE, 0 );

            // get meshes
            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
            tPropEMod->set_parameters( { { { 1000000.0 } } } );
            tPropEMod->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
            tPropNu->set_parameters( { { { 0.3 } } } );
            tPropNu->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet_ss4 = std::make_shared< fem::Property >();
            tPropDirichlet_ss4->set_parameters( { { { 0.0 }, { 0.0 } } } );
            tPropDirichlet_ss4->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet2_ss4 = std::make_shared< fem::Property >();
            tPropDirichlet2_ss4->set_parameters( { { { 1.0, 1.0 } } } );
            tPropDirichlet2_ss4->set_val_function( tMValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 }, { 1.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet2 = std::make_shared< fem::Property >();
            tPropDirichlet2->set_parameters( { { { 1.0, 1.0 } } } );
            tPropDirichlet2->set_val_function( tMValFunction );

            std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
            tPropTraction->set_parameters( { { { 1000.0 }, { 100.0 } } } );
            tPropTraction->set_val_function( tConstValFunction );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tCMStrucLinIso1->set_property( tPropEMod, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropNu, "PoissonRatio" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso1->set_space_dim( 2 );
            tCMStrucLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tCMStrucLinIso2->set_property( tPropEMod, "YoungsModulus" );
            tCMStrucLinIso2->set_property( tPropNu, "PoissonRatio" );
            tCMStrucLinIso2->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso2->set_space_dim( 2 );
            tCMStrucLinIso2->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitsche->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
            tSPNitscheInterface->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );
            tSPNitscheInterface->set_property( tPropEMod, "Material", mtk::Leader_Follower::FOLLOWER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulk2 =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk2->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );

            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet2, "Select", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichletFixed =
                    tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichletFixed->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichletFixed->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichletFixed->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichletFixed->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichletFixed->set_property( tPropDirichlet_ss4, "Dirichlet", mtk::Leader_Follower::LEADER );
            tIWGDirichletFixed->set_property( tPropDirichlet2_ss4, "Select", mtk::Leader_Follower::LEADER );

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
            tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

            // create a list of active block-sets
            std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk2.set_IWGs( { tIWGBulk2 } );

            fem::Set_User_Info tSetBulk3;
            tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk3.set_IWGs( { tIWGBulk1 } );

            fem::Set_User_Info tSetBulk4;
            tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk4.set_IWGs( { tIWGBulk1 } );

            fem::Set_User_Info tSetDirichletFixed;
            tSetDirichletFixed.set_mesh_set_name( "SideSet_4_n_p1" );
            tSetDirichletFixed.set_IWGs( { tIWGDirichletFixed } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            //        fem::Set_User_Info tSetInterface;
            //        tSetInterface.set_mesh_set_name( tDblInterfaceSideSetName );
            //        tSetInterface.set_IWGs( { tIWGInterface } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 4 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet;
            tSetInfo( 3 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0,
                    false );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
            tDofTypesU( 0 ) = MSI::Dof_Type::UX;
            tDofTypesU( 1 ) = MSI::Dof_Type::UY;

            dla::Solver_Factory                             tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param( "AZ_diagnostics" )     = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_output" )          = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_max_iter" )        = 10000;
            tLinearSolverAlgorithm->set_param( "AZ_solver" )          = AZ_gmres;
            tLinearSolverAlgorithm->set_param( "AZ_subdomain_solve" ) = AZ_ilu;
            tLinearSolverAlgorithm->set_param( "AZ_graph_fill" )      = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param( "NLA_max_iter" ) = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

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

            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            // output solution and meshes
            // FIXME if needed

            // clean up
            delete tModel;
            delete tInterpolationMesh;
        }
    }

    TEST_CASE( "2D XTK WITH HMR Struc 2D second", "[XTK_HMR_Struc_2D_02]" )
    {
        if ( par_size() <= 1 )
        {
            uint        tLagrangeMeshIndex = 0;
            std::string tFieldName         = "Cylinder";

            moris::ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "20, 20" );
            tParameters.set( "domain_dimensions", "2, 2" );
            tParameters.set( "domain_offset", "-1.0, -1.0" );
            tParameters.set( "domain_sidesets", "1,2,3,4" );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "1" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            //// create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LvlSetCircle_2D );
            //
            for ( uint k = 0; k < 2; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LvlSetCircle_2D );
            }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto tCircle = std::make_shared< moris::gen::Circle >( 0.0, 0.0, 0.2505 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };

            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( 2, tInterpolationMesh, &tGeometryEngine );

            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::NODE, 0 );

            // get meshes
            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
            tPropEMod->set_parameters( { { { 1000000.0 } } } );
            tPropEMod->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
            tPropNu->set_parameters( { { { 0.5 } } } );
            tPropNu->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 0.0 }, { 0.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichlet2 = std::make_shared< fem::Property >();
            tPropDirichlet2->set_parameters( { { { 1.0 }, { 1.0 } } } );
            tPropDirichlet2->set_val_function( tMValFunction );

            std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
            tPropTraction->set_parameters( { { { 1000.0 }, { 100.0 } } } );
            tPropTraction->set_val_function( tConstValFunction );

            //        std::shared_ptr< fem::Property > tPropDirichletUX = std::make_shared< fem::Property >();
            //        tPropDirichletUX->set_parameters( { {{ 0.0 }} } );
            //        tPropDirichletUX->set_val_function( tConstValFunction );
            //        tPropDirichletUX->set_dof_type( MSI::Dof_Type::UX );
            //        std::shared_ptr< fem::Property > tPropDirichletUY = std::make_shared< fem::Property >();
            //        tPropDirichletUY->set_parameters( { {{ 0.0 }} } );  // specify UY displacement
            //        tPropDirichletUY->set_val_function( tConstValFunction );
            //        tPropDirichletUY->set_dof_type( MSI::Dof_Type::UY );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tCMStrucLinIso->set_property( tPropEMod, "YoungsModulus" );
            tCMStrucLinIso->set_property( tPropNu, "PoissonRatio" );
            tCMStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso->set_space_dim( 2 );
            tCMStrucLinIso->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitsche->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulk->set_constitutive_model( tCMStrucLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMStrucLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
            tIWGDirichlet->set_property( tPropDirichlet2, "Select", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

            // create a list of active block-sets
            std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk1.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk2.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 4 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet;
            tSetInfo( 3 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0,
                    false );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
            tDofTypesU( 0 ) = MSI::Dof_Type::UX;
            tDofTypesU( 1 ) = MSI::Dof_Type::UY;

            dla::Solver_Factory                             tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param( "AZ_diagnostics" )     = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_output" )          = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_max_iter" )        = 10000;
            tLinearSolverAlgorithm->set_param( "AZ_solver" )          = AZ_gmres;
            tLinearSolverAlgorithm->set_param( "AZ_subdomain_solve" ) = AZ_ilu;
            tLinearSolverAlgorithm->set_param( "AZ_graph_fill" )      = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param( "NLA_max_iter" ) = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

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

            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            // output solution and meshes
            // FIXME if needed

            // clean up
            delete tModel;
            delete tInterpolationMesh;
        }
    }

}    // namespace moris
