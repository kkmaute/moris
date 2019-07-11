/*
 * UT_MDL_XTK_DiffusionElement.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_Plane.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Multi_Cylinder.hpp"
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

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

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

#include "fn_norm.hpp"


namespace moris
{
TEST_CASE("XTK Cut Diffusion Model","[XTK_DIFF]")
                {
    if(par_size() == 1)
    {
        //        moris::real tXCenter = 1.0;
        //        moris::real tYCenter = 1.0;
        //        moris::real tZCenter = 3.51;
        //        moris::real tXNorm  = 0.0;
        //        moris::real tYNorm  = 0.0;
        //        moris::real tZNorm  = 1.0;
        //        xtk::Plane tPlane(tXCenter, tYCenter, tZCenter, tXNorm, tYNorm, tZNorm);


        Cell<Cell<moris::real>> aCenter = {{5.0, 5.0, 10.0}};
        Cell<Cell<moris::real>> aAxis   = {{0.0,1.0,0.0}};
        Cell<moris::real> aRad = {2.21};
        Cell<moris::real> aLength = {100};

        xtk::Multi_Cylinder tPlane (aCenter,aRad,aLength,aAxis);

        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        xtk::Geometry_Engine tGeometryEngine(tPlane,tPhaseTable);

        // declare solution field on mesh
        // Declare scalar node field
        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
        std::string tFieldName1 = "temp";
        tNodeField1.set_field_name(tFieldName1);
        tNodeField1.set_field_entity_rank(EntityRank::NODE);

        // Initialize field information container
        mtk::MtkFieldsInfo tFieldsInfo;

        // Place the node field into the field info container
        mtk::add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

        // Declare some supplementary fields
        mtk::MtkMeshData tMeshData;
        tMeshData.FieldsInfo = &tFieldsInfo;

        // Create Mesh --------------------------------------------------------------------
        std::string tMeshFileName = "generated:10x10x20|sideset:zZ";
        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tMeshData );

        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 3;
        xtk::Model tXTKModel(tModelDimension,tInterpMesh1,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        // output to exodus file ----------------------------------------------------------
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets       = true;
        tOutputOptions.mAddSideSets       = true;
        tOutputOptions.mAddClusters       = true;
        tOutputOptions.mAddParallelFields = true;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tXTKIntegFile = "./mdl_exo/xtk_mdl_bar_hole_no_sol.e";
        tIntegMesh1->create_output_mesh(tXTKIntegFile);


        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh1, tIntegMesh1);

        // create a list of IWG type
        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

        // create a list of active block-sets
        moris::Cell< moris_index >  tBlocksetList = { 4, 5 };

        // create a list of active side-sets
        moris::Cell< moris_index >  tSidesetList = { 1, 3 };

        // create a list of BC type for the side-sets
        moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
                fem::BC_Type::NEUMANN };

        // create a list of active double side-sets
        moris::Cell< moris_index >  tDoubleSidesetList = {  };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, 1, tIWGTypeList,
                tBlocksetList, tSidesetList,
                tSidesetBCTypeList,
                tDoubleSidesetList );

        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );

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

        NLA::SOL_Warehouse tSolverWarehouse;

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


        Matrix<DDRMat> tGoldSolution =  {{+2.50e+01},
                {+2.50e+01},
                {+2.50e+01},
                {+2.50e+01},
                {+4.50e+01},
                {+4.50e+01},
                {+4.50e+01},
                {+4.50e+01},
                {+6.50e+01},
                {+6.50e+01},
                {+6.50e+01},
                {+6.50e+01},
                {+8.50e+01},
                {+8.50e+01},
                {+8.50e+01},
                {+8.50e+01},
                {+5.00e+00},
                {+5.00e+00},
                {+5.00e+00},
                {+5.00e+00}};

        moris::print(tSolution11,"tSolution11");

        std::cout<<"Min = "<<tSolution11.min()<<std::endl;
        // verify solution
        //        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);

        // output solution and meshes
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);


        std::string tMeshOutputFile = "./mdl_exo/xtk_bar_mesh.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tInterpMesh1;
        delete tIntegMesh1;
    }
                }

//TEST_CASE("STK Interp XTK Cut Diffusion Model Petsc","[XTK_DIFF_PETSC]")
//{
//    if(par_size() == 1)
//    {
//        Cell<Cell<moris::real>> aCenter = {{5.0, 5.0, 10.0}};
//        Cell<Cell<moris::real>> aAxis   = {{0.0,1.0,0.0}};
//        Cell<moris::real> aRad = {2.21};
//        Cell<moris::real> aLength = {100};
//
//        xtk::Multi_Cylinder tPlane (aCenter,aRad,aLength,aAxis);
//
//        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//        xtk::Geometry_Engine tGeometryEngine(tPlane,tPhaseTable);
//
//        // declare solution field on mesh
//        // Declare scalar node field
//        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
//        std::string tFieldName1 = "temp";
//        tNodeField1.set_field_name(tFieldName1);
//        tNodeField1.set_field_entity_rank(EntityRank::NODE);
//
//        // Initialize field information container
//        mtk::MtkFieldsInfo tFieldsInfo;
//
//        // Place the node field into the field info container
//        mtk::add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
//
//        // Declare some supplementary fields
//        mtk::MtkMeshData tMeshData;
//        tMeshData.FieldsInfo = &tFieldsInfo;
//
//        // Create Mesh --------------------------------------------------------------------
//        std::string tMeshFileName = "generated:10x10x20|sideset:zZ";
//        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tMeshData );
//
//        // Setup XTK Model ----------------------------------------------------------------
//        size_t tModelDimension = 3;
//        xtk::Model tXTKModel(tModelDimension,tInterpMesh1,tGeometryEngine);
//        tXTKModel.mVerbose = true;
//
//        //Specify decomposition Method and Cut Mesh ---------------------------------------
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//        tXTKModel.decompose(tDecompositionMethods);
//
//        // output to exodus file ----------------------------------------------------------
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets       = true;
//        tOutputOptions.mAddSideSets       = true;
//        tOutputOptions.mAddClusters       = true;
//        tOutputOptions.mAddParallelFields = true;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//        std::string tXTKIntegFile = "./mdl_exo/xtk_mdl_bar_hole_no_sol.e";
//        tIntegMesh1->create_output_mesh(tXTKIntegFile);
//
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair(tInterpMesh1, tIntegMesh1);
//
//        // create a list of IWG type
//        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
//        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
//        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
//        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );
//
//        // create a list of active block-sets
//        moris::Cell< moris_index >  tBlocksetList = { 4, 5 };
//
//        // create a list of active side-sets
//        moris::Cell< moris_index >  tSidesetList = { 1, 3 };
//
//        // create a list of BC type for the side-sets
//        moris::Cell< fem::BC_Type > tSidesetBCTypeList = { fem::BC_Type::DIRICHLET,
//                fem::BC_Type::NEUMANN };
//
//        // create a list of active double side-sets
//        moris::Cell< moris_index >  tDoubleSidesetList = {  };
//
//        // create model
//        mdl::Model * tModel = new mdl::Model( &tMeshManager, 1, tIWGTypeList,
//                tBlocksetList, tSidesetList,
//                tSidesetBCTypeList,
//                tDoubleSidesetList );
//
//        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 1: create linear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        dla::Solver_Factory  tSolFactory;
//
//        // create linear solver
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );
//
//        //        tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPFGMRES );
//        //        //tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCMG );
//        //        tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCILU );
//        //        tLinearSolverAlgorithm->set_param("ILUFill")  = 3;
//
//        dla::Linear_Solver tLinSolver;
//
//        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 2: create nonlinear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );
//
//        // create factory for nonlinear solver
//        NLA::Nonlinear_Solver_Factory tNonlinFactory;
//
//        // create nonlinear solver
//        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        NLA::Nonlinear_Solver  tNonlinearSolver;
//
//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//        // set manager and settings
//        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//        tNonlinearSolver.solve( tNonlinearProblem );
//
//        //        // temporary array for solver
//        //        Matrix< DDRMat > tSolution;
//        //        tNonlinearSolver.get_full_solution( tSolution );
//
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 3: create time Solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        //        tsa::Time_Solver_Factory tTimeSolverFactory;
//        //        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//        //
//        //        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//        //
//        //        tsa::Time_Solver tTimeSolver;
//        //
//        //        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//        //
//        //        NLA::SOL_Warehouse tSolverWarehouse;
//        //
//        //        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//        //
//        //        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//        //        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//        //
//        //        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
//        //        tTimeSolver.set_dof_type_list( tDofTypes1 );
//        //
//        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        //        // STEP 4: Solve and check
//        //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        //
//        //        tTimeSolver.solve();
//
//
//        // TODO: add gold solution data for this problem
//
//        // Write to Integration mesh for visualization
//        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//
//        // add solution field to integration mesh
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);
//
//
//        // verify solution
//        //    CHECK(norm(tSolution11 - tGoldSolution)<1e-08);
//
//        // output solution and meshes
//        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ.e";
//        tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//        //    delete tInterpMesh1;
//        delete tModel;
//        delete tIntegMesh1;
//    }
//}


}
