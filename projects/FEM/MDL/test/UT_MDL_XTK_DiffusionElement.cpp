/*
 * UT_MDL_XTK_DiffusionElement.cpp
 *
 *  Created on: Jun 18, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
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
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src

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

Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                    moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator )
{
    return aCoeff( 0 );
}

TEST_CASE("XTK Cut Diffusion Model","[XTK_DIFF]")
                {
    if(par_size() == 1)
    {
        moris::real tXCenter = 1.0;
        moris::real tYCenter = 1.0;
        moris::real tZCenter = 3.1;
        moris::real tXNorm  = 0.0;
        moris::real tYNorm  = 0.0;
        moris::real tZNorm  = 1.0;
        xtk::Plane tPlane(tXCenter, tYCenter, tZCenter, tXNorm, tYNorm, tZNorm);
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
        std::string tMeshFileName = "generated:1x1x4|sideset:z";
        moris::mtk::Interpolation_Mesh* tInterpMesh1 = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tMeshData );

        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 3;
        xtk::Model tXTKModel(tModelDimension,tInterpMesh1,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::NODE);

        // get meshes
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        std::string tDirchletSideName      = "surface_1_n_p0";
        std::string tNeumannSideName       = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
        std::string tBulkBlockNamesChild   = "block_1_c_p0";
        std::string tBulkBlockNamesNoChild = "block_1_n_p0";

        std::cout<<"Dirichlet = "<<tDirchletSideName<<std::endl;
        std::cout<<"Neumann = "<<tNeumannSideName<<std::endl;


        // output to exodus file ----------------------------------------------------------
        //        xtk::Output_Options tOutputOptions;
        //        tOutputOptions.mAddNodeSets       = true;
        //        tOutputOptions.mAddSideSets       = true;
        //        tOutputOptions.mAddClusters       = true;
        //        tOutputOptions.mAddParallelFields = true;
        //
        //        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);


        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create a list of IWG type
        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGTypeList( 3 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

        // list of property type
        Cell< fem::Property_Type > tPropertyTypeList = {{ fem::Property_Type::CONDUCTIVITY   },
                                                        { fem::Property_Type::TEMP_DIRICHLET },
                                                        { fem::Property_Type::TEMP_NEUMANN   }};

        // list of property dependencies
        Cell< Cell< Cell< MSI::Dof_Type > > > tPropertyDofList( 3 );

        // list of the property coefficients
        Cell< Cell< Matrix< DDRMat > > > tCoeffList( 3 );
        tCoeffList( 0 ).resize( 1 );
        tCoeffList( 0 )( 0 )= {{ 1.0 }};
        tCoeffList( 1 ).resize( 1 );
        tCoeffList( 1 )( 0 )= {{ 5.0 }};
        tCoeffList( 2 ).resize( 1 );
        tCoeffList( 2 )( 0 )= {{ 20.0 }};

        // cast free function into std::function
        fem::PropertyFunc tValFunction0 = tConstValFunction;

        // create the list with function pointers for the value
        Cell< fem::PropertyFunc > tValFuncList( 3, tValFunction0 );

        // create the list with cell of function pointers for the derivatives
        Cell< Cell< fem::PropertyFunc > > tDerFuncList( 3 );

        // collect properties info
        fem::Property_User_Defined_Info tPropertyUserDefinedInfo( tPropertyTypeList,
                                                                  tPropertyDofList,
                                                                  tCoeffList,
                                                                  tValFuncList,
                                                                  tDerFuncList );

        // create a list of active block-sets
        moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index(tBulkBlockNamesNoChild),
                                                 tEnrIntegMesh.get_block_set_index(tBulkBlockNamesChild),
                                                 tEnrIntegMesh.get_side_set_index(tDirchletSideName),
                                                 tEnrIntegMesh.get_side_set_index(tNeumannSideName) };

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, 1, tIWGTypeList,
                                              tSetList, tSetTypeList,
                                              &tPropertyUserDefinedInfo );

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


        // verify solution
        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);

        // output solution and meshes

        tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName1,EntityRank::NODE,tSolution11);

        std::string tOutputInterp = "./mdl_exo/xtk_mdl_interp.exo";
        tInterpMesh1->create_output_mesh(tOutputInterp);

        //        std::string tMeshOutputFile = "./mdl_exo/xtk_bar_mesh.e";
        //        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tInterpMesh1;
        delete tModel;
        //        delete tIntegMesh1;
    }
                }

TEST_CASE("XTK STK Cut Diffusion Model","[XTK_STK_DIFF]")
{

    if(par_size() == 1)
    {
        moris::real tXCenter = 1.0;
        moris::real tYCenter = 1.0;
        moris::real tZCenter = 3.51;
        moris::real tXNorm  = 0.0;
        moris::real tYNorm  = 0.0;
        moris::real tZNorm  = 1.0;
        xtk::Plane tPlane(tXCenter, tYCenter, tZCenter, tXNorm, tYNorm, tZNorm);
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
        std::string tMeshFileName = "generated:1x1x4|sideset:z";
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

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);


        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(tInterpMesh1, tIntegMesh1);

        // create a list of IWG type
        Cell< Cell< fem::IWG_Type > >tIWGTypeList( 4 );
        tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
        tIWGTypeList( 3 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

        // list of property type
        Cell< fem::Property_Type > tPropertyTypeList = {{ fem::Property_Type::CONDUCTIVITY   },
                                                        { fem::Property_Type::TEMP_DIRICHLET },
                                                        { fem::Property_Type::TEMP_NEUMANN   }};

        // list of property dependencies
        Cell< Cell< Cell< MSI::Dof_Type > > > tPropertyDofList( 3 );

        // list of the property coefficients
        Cell< Cell< Matrix< DDRMat > > > tCoeffList( 3 );
        tCoeffList( 0 ).resize( 1 );
        tCoeffList( 0 )( 0 )= {{ 1.0 }};
        tCoeffList( 1 ).resize( 1 );
        tCoeffList( 1 )( 0 )= {{ 5.0 }};
        tCoeffList( 2 ).resize( 1 );
        tCoeffList( 2 )( 0 )= {{ 20.0 }};

        // cast free function into std::function
        fem::PropertyFunc tValFunction0 = tConstValFunction;

        // create the list with function pointers for the value
        Cell< fem::PropertyFunc > tValFuncList( 3, tValFunction0 );

        // create the list with cell of function pointers for the derivatives
        Cell< Cell< fem::PropertyFunc > > tDerFuncList( 3 );

        // collect properties info
        fem::Property_User_Defined_Info tPropertyUserDefinedInfo( tPropertyTypeList,
                                                                  tPropertyDofList,
                                                                  tCoeffList,
                                                                  tValFuncList,
                                                                  tDerFuncList );


        // create a list of active block-sets
        moris::Cell< moris_index >  tSetList = { 4, 5, 1, 0 };

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, 1, tIWGTypeList,
                                              tSetList, tSetTypeList,
                                              &tPropertyUserDefinedInfo );

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

        //moris::print(tSolution11,"tSolution11");

        // verify solution
        CHECK(norm(tSolution11 - tGoldSolution) < 1e-08);

        // output solution and meshes

        tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName1,EntityRank::NODE,tSolution11);

        std::string tOutputInterp = "./mdl_exo/xtk_mdl_interp.exo";
        tInterpMesh1->create_output_mesh(tOutputInterp);

        std::string tMeshOutputFile = "./mdl_exo/xtk_bar_mesh.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tInterpMesh1;
        delete tIntegMesh1;
        delete tModel;
    }
}


}
