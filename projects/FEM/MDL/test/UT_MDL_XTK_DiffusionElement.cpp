/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_DiffusionElement.cpp
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

#include "cl_FEM_NodeProxy.hpp"                  //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"               //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                  //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"            //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"                //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"                //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"                //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp" //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

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

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Plane.hpp"

namespace moris
{

    void tConstValFunction_MDL_XTK
    ( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    TEST_CASE("XTK Cut Diffusion Model","[XTK_DIFF]")
    {

        if(par_size() == 1)
        {
            moris::Matrix<moris::DDRMat> tCenters = {{ 1.0,1.0,3.1 }};
            moris::Matrix<moris::DDRMat> tNormals = {{ 0.0,0.0,1.0 }};
            Cell< std::shared_ptr< moris::ge::Geometry > > tGeometry( 1 );
            auto tField = std::make_shared<moris::ge::Plane>(tCenters(0), tCenters(1), tCenters(2), tNormals(0), tNormals(1), tNormals(2));
            tGeometry( 0 ) = std::make_shared< ge::Level_Set_Geometry >( tField );

            // Initialize field information container
            mtk::MtkFieldsInfo tFieldsInfo;

            // Declare some supplementary fields
            mtk::MtkMeshData tMeshData;
            tMeshData.FieldsInfo = &tFieldsInfo;

            // Create Mesh --------------------------------------------------------------------
            std::string tMeshFileName = "generated:1x1x4|sideset:z";
            moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tMeshData );

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGeometryEngine(tInterpMesh1, tGeometryEngineParameters);

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            xtk::Model tXTKModel(tModelDimension,tInterpMesh1,&tGeometryEngine);

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::NODE );

            // get meshes
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

            std::string tDirchletSideName      = "surface_1_n_p0";
            std::string tNeumannSideName       = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
            std::string tBulkBlockNamesChild   = "block_1_c_p0";
            std::string tBulkBlockNamesNoChild = "block_1_n_p0";

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
            tPropConductivity->set_parameters( { {{ 1.0 }} } );
            tPropConductivity->set_val_function( tConstValFunction_MDL_XTK );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { {{ 5.0 }} } );
            tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { {{ 20.0 }} } );
            tPropNeumann->set_val_function( tConstValFunction_MDL_XTK );

            std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
            tPropTempLoad->set_parameters( { {{ 0.0 }} } );
            tPropTempLoad->set_val_function( tConstValFunction_MDL_XTK );

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
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( tBulkBlockNamesNoChild );
            tSetBulk1.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( tBulkBlockNamesChild );
            tSetBulk2.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( tDirchletSideName );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( tNeumannSideName );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet;
            tSetInfo( 3 ) = tSetNeumann;

            // create model
            mdl::Model * tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD, //STANDARD_WITH_OVERLAP
                    "./",
                    "UT_Output_xtk_mdl_enr_integ.exo",
                    "./",
                     "temp.exo",
                    { "block_1_c_p0", "block_1_n_p0" },
                    { "TEMP" },
                    { vis::Field_Type::NODAL },
                    { "IQI_temp" } );
            tModel->set_output_manager( &tOutputData );

            moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
            tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
            tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
            tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

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

            Matrix<DDRMat> tGoldSolution = {
            {2.499999999999997e+01},
            {2.499999999999997e+01},
            {2.499999999999996e+01},
            {2.499999999999997e+01},
            {4.499999999999995e+01},
            {4.499999999999996e+01},
            {4.499999999999994e+01},
            {4.499999999999995e+01},
            {6.499999999999994e+01},
            {6.499999999999993e+01},
            {6.499999999999993e+01},
            {6.499999999999996e+01},
            {8.500000000000014e+01},
            {8.500000000000011e+01},
            {8.500000000000013e+01},
            {8.499999999999947e+01},
            {5.000000000000000e+00},
            {5.000000000000005e+00},
            {5.000000000000003e+00},
            {5.000000000000001e+00}};

            // verify solution
            // define a bool for solution check
            bool tCheckNodalSolution = true;

            // number of mesh nodes
            uint tNumOfNodes = tInterpMesh1->get_num_nodes();

            // loop over the node and check solution
            for ( uint i = 0; i < tNumOfNodes; i++ )
            {
                // check solution
                tCheckNodalSolution = tCheckNodalSolution
                        && ( std::abs( tSolution11( i ) - tGoldSolution( i ) ) < 1e-08 );
            }
            // check bool is true
            REQUIRE( tCheckNodalSolution );

            // clean up
            delete tInterpMesh1;
            delete tModel;
        }
    }
}

