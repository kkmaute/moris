/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_Symm_BCs.cpp
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
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
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

#include "cl_GEN_Circle.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris
{

void tConstValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tMValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                   0.0 },
            { 0.0,                   aParameters( 0 )( 1 ) }};
}

moris::real LvlSetCircle_2D_outsideDomain(const moris::Matrix< moris::DDRMat > & aPoint )
{
    Matrix< DDRMat > tCenter{ { -10, -10 } };
    return    norm( aPoint - tCenter ) - 0.001;
}

TEST_CASE("2D XTK WITH HMR SYMM BCs","[XTK_HMR_2D_Symm_BCs]")
{
    if(par_size()<=1)
    {
        uint tLagrangeMeshIndex = 0;
        std::string tFieldName = "Cylinder";

        moris::Parameter_List tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", 600, 6 );
        tParameters.set( "domain_dimensions", 100.0, 2.0 );
//        tParameters.set( "domain_offset", "-50.0, -1.0" );
        tParameters.set( "lagrange_output_meshes", "0" );

        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_to_bspline", "0" );

        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", "0" );

        tParameters.set( "severity_level", 2 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement();

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //  create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( LvlSetCircle_2D_outsideDomain );

        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( LvlSetCircle_2D_outsideDomain );
        }

        tHMR.finalize();

        moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        //-----------------------------------------------------------------------------------------------

        auto tCircle = std::make_shared< moris::gen::Circle >( -100.0, -100.0, 0.001 );
        Vector<std::shared_ptr<moris::gen::Geometry>> tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };

        moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::gen::Geometry_Engine tGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);

         xtk::Model tXTKModel(2, tInterpolationMesh, &tGeometryEngine);

        tXTKModel.mVerbose = false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(mtk::EntityRank::BSPLINE,0);

        // get meshes
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

//        tEnrInterpMesh.print_enriched_cells();
//        tEnrIntegMesh.print_double_side_sets(2);

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
        tPropEMod->set_parameters( { {{ 100.0 }} } );
        tPropEMod->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
        tPropNu->set_parameters( { {{ 0.0 }} } );
        tPropNu->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropDirichletUXUY_ss4 = std::make_shared< fem::Property >();    // fix displacement at side-set 4 to be zero in x and y
        tPropDirichletUXUY_ss4->set_parameters( { {{ 0.0 }, {0.0}} } );
        tPropDirichletUXUY_ss4->set_val_function( tConstValFunction );
//        tPropDirichletUX_ss4->set_dof_type( MSI::Dof_Type::UX );

        std::shared_ptr< fem::Property > tPropDirichletUXUY_ss4_select = std::make_shared< fem::Property >();    // fix displacement at side-set 4 to be zero in x and y
        tPropDirichletUXUY_ss4_select->set_parameters( { {{ 1.0 }, {1.0}} } );
        tPropDirichletUXUY_ss4_select->set_val_function( tMValFunction );
//        tPropDirichletUX_ss4->set_dof_type( MSI::Dof_Type::UX );

        std::shared_ptr< fem::Property > tPropDirichletUX_ss2 = std::make_shared< fem::Property >();        // allow only for y-displacement at other end (side-set 2)
        tPropDirichletUX_ss2->set_parameters( { {{ 0.0 }, {1.0}} } );                                       // fix x-displ. to zero
        tPropDirichletUX_ss2->set_val_function( tConstValFunction );
//        tPropDirichletUX->set_dof_type( MSI::Dof_Type::UX );

        std::shared_ptr< fem::Property > tPropDirichletUX_ss2_select = std::make_shared< fem::Property >();        // allow only for y-displacement at other end (side-set 2)
        tPropDirichletUX_ss2_select->set_parameters( { {{ 1.0 }, {0.0}} } );                                       // fix x-displ. to zero
        tPropDirichletUX_ss2_select->set_val_function( tMValFunction );
//        tPropDirichletUX->set_dof_type( MSI::Dof_Type::UX );

        std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
        tPropTraction->set_parameters( {{{ 0.0 } , { 10.0 }}} );
        tPropTraction->set_val_function( tConstValFunction );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tCMStrucLinIso1->set_property( tPropEMod, "YoungsModulus" );
        tCMStrucLinIso1->set_property( tPropNu, "PoissonRatio" );
        tCMStrucLinIso1->set_space_dim( 2 );
        tCMStrucLinIso1->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 100.0 }} } );
        tSPDirichletNitsche->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGBulk1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet_ss2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet_ss2->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGDirichlet_ss2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGDirichlet_ss2->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet_ss2->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIWGDirichlet_ss2->set_property( tPropDirichletUX_ss2, "Dirichlet", mtk::Leader_Follower::LEADER );
        tIWGDirichlet_ss2->set_property( tPropDirichletUX_ss2_select, "Select", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichletFixed_ss4 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichletFixed_ss4->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGDirichletFixed_ss4->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGDirichletFixed_ss4->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichletFixed_ss4->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIWGDirichletFixed_ss4->set_property( tPropDirichletUXUY_ss4, "Dirichlet", mtk::Leader_Follower::LEADER );
        tIWGDirichletFixed_ss4->set_property( tPropDirichletUXUY_ss4_select, "Select", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

        // create a list of active block-sets
        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
        std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);

        // define set info
        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk4.set_IWGs( { tIWGBulk1 } );

        fem::Set_User_Info tSetDirichletFixed;
        tSetDirichletFixed.set_mesh_set_name( "SideSet_4_n_p1" );
        tSetDirichletFixed.set_IWGs( { tIWGDirichletFixed_ss4 } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_2_n_p1" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet_ss2 } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk4;
        tSetInfo( 1 ) = tSetDirichletFixed;
        tSetInfo( 2 ) = tSetDirichlet;
        tSetInfo( 3 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                              0,
                                              tSetInfo,
                                              0,
                                              false );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
        tDofTypesU( 0 ) = MSI::Dof_Type::UX;    tDofTypesU( 1 ) = MSI::Dof_Type::UY;

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );
        tLinearSolverParameterList.set( "AZ_max_iter", 100 );
        tLinearSolverParameterList.set( "AZ_solver", AZ_gmres );
        tLinearSolverParameterList.set( "AZ_subdomain_solve", AZ_ilu );
        tLinearSolverParameterList.set( "AZ_graph_fill", 10 );
        //        tLinearSolverParameterList.set( "ml_prec_type", "SA" );

        dla::Linear_Solver tLinSolver;
        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 10 );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
        //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolverMain;
        tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        tNonlinearSolverMain       .set_dof_type_list( tDofTypesU );

        // Create solver database
        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

        tNonlinearSolverMain       .set_solver_warehouse( &tSolverWarehouse );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

        tsa::Time_Solver tTimeSolver;
        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tTimeSolver.set_dof_type_list( tDofTypesU );

        //------------------------------------------------------------------------------
        tTimeSolver.solve();
        //-----  print solution vector -----
//        Matrix< DDRMat > tSolVec;
//        tNonlinearSolverMain.get_full_solution(tSolVec);
//        print( tSolVec, " Solution Vector ");

        delete tModel;
        delete tInterpolationMesh;
    }
}
}   // end moris namespace

