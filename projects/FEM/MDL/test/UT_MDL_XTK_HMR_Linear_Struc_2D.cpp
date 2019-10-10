/*
 * UT_MDL_XTK_HMR_2D.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: schmidt
 */

#include "catch.hpp"
#include "cl_Star.hpp"
#include "cl_Circle.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

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
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src

#include "cl_FEM_Constitutive_User_Defined_Info.hpp"

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

Matrix< DDRMat >
exactTempFunc(moris::Cell< Matrix< DDRMat > >         & aCoeff,
              moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator,
              fem::Geometry_Interpolator             * aGeometryInterpolator )
{
    Matrix< DDRMat > tCoord = aGeometryInterpolator->valx();
    real xcoord = tCoord(0);
    real ycoord = tCoord(1);

    real rad = std::pow (  std::pow( xcoord - 0, 2.0)
                         + std::pow( ycoord - 0, 2.0), 0.5);

    return {{(1.0/3.0)*(1.0/rad-0.501)}};
}

moris::real LvlSetLin(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tOffset = 200;

    return    aPoint(0) - 0.317 * aPoint(1) - tOffset;
}

moris::real LvlSetCircle_1(const moris::Matrix< moris::DDRMat > & aPoint )
{
    return    norm( aPoint ) - 0.501;
}


Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                    moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator,
                                    fem::Geometry_Interpolator             * aGeometryInterpolator )
{
    return aCoeff( 0 );
}

moris::real
LevelSetFunction_star( const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tPhi = std::atan2( aPoint( 0 ), aPoint( 2 ) );

    moris::real tLevelSetVaue = 0.501 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aPoint( 0 ), 2 ) + std::pow( aPoint( 2 ), 2 ) );

    return tLevelSetVaue;
}


TEST_CASE("2D XTK WITH HMR Struc 2D","[XTK_HMR_Struc_2D]")
{
if(par_size()<=1)
{
    uint tLagrangeMeshIndex = 0;
    std::string tFieldName = "Cylinder";

hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

tParameters.set( "number_of_elements_per_dimension", "2, 2");
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
tParameters.set( "initial_refinement", 0 );

tParameters.set( "use_multigrid", 0 );
tParameters.set( "severity_level", 2 );

hmr::HMR tHMR( tParameters );

// initial refinement
tHMR.perform_initial_refinement( 0 );

std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//// create field
std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

tField->evaluate_scalar_function( LvlSetLin );
//
//    for( uint k=0; k<1; ++k )
//    {
//        tHMR.flag_surface_elements_on_working_pattern( tField );
//        tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//        tField->evaluate_scalar_function( LevelSetFunction_star );
//    }

tHMR.finalize();

tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);


//    xtk::Geom_Field tFieldAsGeom(tField);
//
//    moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//    // Tell the geometry engine about the discrete field mesh and how to interpret phases
//    xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//    xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);

    xtk::Circle tCircle( 0.45, 2.0, 0.0 );

    xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    xtk::Geometry_Engine tGeometryEngine(tCircle,tPhaseTable, 2);

     xtk::Model tXTKModel(2, tInterpolationMesh.get(), tGeometryEngine);

     tXTKModel.mVerbose = true;

    //Specify decomposition Method and Cut Mesh ---------------------------------------
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
    tXTKModel.decompose(tDecompositionMethods);

    tXTKModel.perform_basis_enrichment(EntityRank::NODE,0);

    // get meshes
    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);


    //---------------------------------------------------------------------------------------


//        // output solution and meshes
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = false;
//        tOutputOptions.mAddClusters = false;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

//        // Write to Integration mesh for visualization
//        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//
//        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
//        {
//            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//            tSTKIntegSol(i) = tIntegSol(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//        }

        // crate field in integration mesh
//        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("solution",EntityRank::NODE);
//        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);
//
//        // add solution field to integration mesh
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);

//        std::string tMeshOutputFile = "./mdl_exo/stk_xtk_inv_ilu_quad_bspline.e";

//        tIntegMesh1->create_output_mesh(tMeshOutputFile);

    //-------------------------------------------------------------------------------------------

    uint tSpatialDimension = 2;

    // create IWG user defined info
    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 4 );
    tIWGUserDefinedInfo( 0 ).resize( 1 );
    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_BULK,
                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY },
                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }},
                                                                moris::Cell< fem::Property_Type >( 0 ),
                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
    tIWGUserDefinedInfo( 1 ).resize( 1 );
    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_BULK,
                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY },
                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }},
                                                                moris::Cell< fem::Property_Type >( 0 ),
                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );

    tIWGUserDefinedInfo( 2 ).resize( 1 );
    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_DIRICHLET,
                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY },
                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }},
                                                                { fem::Property_Type::STRUC_DIRICHLET },
                                                                { fem::Constitutive_Type::STRUC_LIN_ISO } );
    tIWGUserDefinedInfo( 3 ).resize( 1 );
    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_NEUMANN,
                                                                { MSI::Dof_Type::UX, MSI::Dof_Type::UY },
                                                                {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY },},
                                                                { fem::Property_Type::STRUC_NEUMANN },
                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );

    // create property user defined info
    fem::Property_User_Defined_Info tYoungs_Modulus( fem::Property_Type::YOUNGS_MODULUS,
                                                   Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                   {{{ 1000000.0 }}},
                                                   tConstValFunction,
                                                   Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tPoissons_Ratio( fem::Property_Type::POISSONS_RATIO,
                                                   Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                   {{{ 0.0 }}},
                                                   tConstValFunction,
                                                   Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tStrucDirichlet( fem::Property_Type::STRUC_DIRICHLET,
                                                    Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                    {{{ 0.0, 0.0 }}},
                                                    tConstValFunction,
                                                    Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tStrucNeumann( fem::Property_Type::STRUC_NEUMANN,
                                                  Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                  {{{ 1000.0, 0.0 }}},
                                                  tConstValFunction,
                                                  Cell< fem::PropertyFunc >( 0 ) );

    // create property user defined info
    Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 4 );
    tPropertyUserDefinedInfo( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 0 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tYoungs_Modulus;
    tPropertyUserDefinedInfo( 0 )( 0 )( 1 ) = tPoissons_Ratio;
    tPropertyUserDefinedInfo( 1 ).resize( 1 );
    tPropertyUserDefinedInfo( 1 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tYoungs_Modulus;
    tPropertyUserDefinedInfo( 1 )( 0 )( 1 ) = tPoissons_Ratio;
    tPropertyUserDefinedInfo( 2 ).resize( 1 );
    tPropertyUserDefinedInfo( 2 )( 0 ).resize( 3 );
    tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tYoungs_Modulus;
    tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tPoissons_Ratio;
    tPropertyUserDefinedInfo( 2 )( 0 )( 2 ) = tStrucDirichlet;
    tPropertyUserDefinedInfo( 3 ).resize( 1 );
    tPropertyUserDefinedInfo( 3 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tStrucNeumann;

    fem::Constitutive_User_Defined_Info tStrucLinIso( fem::Constitutive_Type::STRUC_LIN_ISO,
                                                     {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY  }},
                                                     { fem::Property_Type::YOUNGS_MODULUS, fem::Property_Type::POISSONS_RATIO } );

    // create constitutive user defined info
    Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 4 );
    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 0 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 0 )( 0 )( 0 ) = tStrucLinIso;
    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 1 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 1 )( 0 )( 0 ) = tStrucLinIso;
    tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 2 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 2 )( 0 )( 0 ) = tStrucLinIso;
    tConstitutiveUserDefinedInfo( 3 ).resize( 1 );


    // create a list of active block-sets
    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );


    // create a list of active block-sets
     moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1")};


    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                      fem::Element_Type::BULK,
                                                      fem::Element_Type::SIDESET,
                                                      fem::Element_Type::SIDESET};


    uint tBSplineMeshIndex = 0;
    // create model
    mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                          tBSplineMeshIndex,
                                          tSetList, tSetTypeList,
                                          tIWGUserDefinedInfo,
                                          tPropertyUserDefinedInfo,
                                          tConstitutiveUserDefinedInfo );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 1: create linear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );

    tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
    tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
    tLinearSolverAlgorithm->set_param("AZ_max_iter") = 1000;
    tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
    tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
    tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 5;

//    tLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;

    dla::Linear_Solver tLinSolver;

    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 2: create nonlinear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Epetra );
//    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );

    NLA::Nonlinear_Solver_Factory tNonlinFactory;
    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

    NLA::Nonlinear_Solver tNonlinearSolver;

    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

    tNonlinearSolver.solve( tNonlinearProblem );


//        Matrix<DDRMat> tFullSol;
        tNonlinearProblem->get_full_vector()->print();

//        print(tFullSol,"tFullSol");


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 3: create time Solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    tsa::Time_Solver_Factory tTimeSolverFactory;
//    std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//    tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//    tsa::Time_Solver tTimeSolver;
//
//    tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//    NLA::SOL_Warehouse tSolverWarehouse;
//
//    tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//    tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//    tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//    tNonlinearSolver.set_dof_type_list( tDofTypes1 );
//    tTimeSolver.set_dof_type_list( tDofTypes1 );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 4: Solve and check
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//    tTimeSolver.solve();

//    MSI::MSI_Solver_Interface * tSolver_Interface = tModel->get_solver_interface();

//    tSolver_Interface->write_solution_to_hdf5_file( "Exact_Sol_petsc.h5" );
//
//    char SolVector[100];
//    std::strcpy( SolVector, "Exact_Sol_petsc.h5" );
//
//    tSolver_Interface->get_exact_solution_from_hdf5_and_calculate_error( SolVector );
//
////    tSolver_Interface->get_residual_vector_for_output("Res_vec.h5");
//
    // output solution and meshes
    xtk::Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = false;
    tOutputOptions.mAddClusters = false;

    // add solution field to integration mesh
    std::string tIntegSolFieldName = "solution";
    tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

    moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

    // Write to Integration mesh for visualization
    Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );

    Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

    for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
    {
        moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
        tSTKIntegSol(i) = tIntegSol(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
    }

    // crate field in integration mesh
    moris::moris_index tFieldIndexUX = tEnrIntegMesh.create_field("UX",EntityRank::NODE);
    moris::moris_index tFieldIndexUY = tEnrIntegMesh.create_field("UY",EntityRank::NODE);
    tEnrIntegMesh.add_field_data(tFieldIndexUX,EntityRank::NODE,tSTKIntegSol);
    tEnrIntegMesh.add_field_data(tFieldIndexUY,EntityRank::NODE,tSTKIntegSol);

    // add solution field to integration mesh
    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);

//    Matrix<DDRMat> tFullSol;
//    tNonlinearSolver.get_full_solution(tFullSol);
//
//    print(tFullSol,"tFullSol");


    std::string tMeshOutputFile = "./mdl_exo/stk_xtk_inv_ilu_quad_bspline.e";

    tIntegMesh1->create_output_mesh(tMeshOutputFile);

    delete tIntegMesh1;

    delete tModel;

}
}
}



