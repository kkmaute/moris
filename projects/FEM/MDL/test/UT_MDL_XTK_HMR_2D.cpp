/*
 * UT_MDL_XTK_HMR_2D.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: schmidt
 */

#include "catch.hpp"
#include "cl_Star.hpp"

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
    moris::real tOffset = -0.6;

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

TEST_CASE("2D XTK WITH HMR 2D","[XTK_HMR_2D]")
{
//    if(par_size()<=1)
//    {
//        std::string tFieldName = "Cylinder";
//
//         moris::uint tLagrangeMeshIndex = 0;
//         moris::uint tBSplineMeshIndex = 0;
//
//         hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//
//         tParameters.set( "number_of_elements_per_dimension", "5, 2" );
//         tParameters.set( "domain_dimensions", "5, 2" );
//         tParameters.set( "domain_offset", "-2.5, -1.0" );
//         tParameters.set( "domain_sidesets", "2, 4" );
//         tParameters.set( "lagrange_output_meshes", "0" );
//
//         tParameters.set( "lagrange_orders", "1" );
//         tParameters.set( "lagrange_pattern", "0" );
//         tParameters.set( "bspline_orders", "1" );
//         tParameters.set( "bspline_pattern", "0" );
//
//         tParameters.set( "lagrange_to_bspline", "0" );
//
//         tParameters.set( "truncate_bsplines", 1 );
//         tParameters.set( "refinement_buffer", 2 );
//         tParameters.set( "staircase_buffer", 2 );
//         tParameters.set( "initial_refinement", 2 );
//
//         tParameters.set( "use_multigrid", 0 );
//         tParameters.set( "severity_level", 2 );
//
//         hmr::HMR tHMR( tParameters );
//
//         // initial refinement
//         tHMR.perform_initial_refinement( 0 );
//
//         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//         // create field
//         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//         tField->evaluate_scalar_function( LvlSetLin );
//
////         for( uint k=0; k<2; ++k )
////         {
////             tHMR.flag_surface_elements_on_working_pattern( tField );
////             tHMR.perform_refinement_based_on_working_pattern( 0 );
////
////             tField->evaluate_scalar_function( CircleFunc );
////         }
//
//         tHMR.finalize();
//
//         tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );
//
//         std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//         xtk::Geom_Field tFieldAsGeom(tField);
//
//         moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//         size_t tModelDimension = 2;
//         xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//         xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
//         xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
//         tXTKModel.mVerbose = true;
//
//        //Specify decomposition Method and Cut Mesh ---------------------------------------
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//        tXTKModel.decompose(tDecompositionMethods);
//
//        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);
//
//        // get meshes
//        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//        // output to exodus file ----------------------------------------------------------
//        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();
//
//         // Declare the fields related to enrichment strategy in output options
//         Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
//
//        // output solution and meshes
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = false;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh1);
//
//        std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_cut.e";
//        tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//        // end output -------------------------------------------------------
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
//
//        uint tSpatialDimension = 2;
//
//        // create IWG user defined info
//        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 6 );
//        tIWGUserDefinedInfo( 0 ).resize( 1 );
//        tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
//                                                                    {{ MSI::Dof_Type::TEMP }},
//                                                                    { fem::Property_Type::CONDUCTIVITY },
//                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
//        tIWGUserDefinedInfo( 1 ).resize( 1 );
//        tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
//                                                                    {{ MSI::Dof_Type::TEMP }},
//                                                                    { fem::Property_Type::CONDUCTIVITY },
//                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//        tIWGUserDefinedInfo( 2 ).resize( 1 );
//        tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
//                                                                    {{ MSI::Dof_Type::TEMP }},
//                                                                    { fem::Property_Type::CONDUCTIVITY },
//                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
//        tIWGUserDefinedInfo( 3 ).resize( 1 );
//        tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
//                                                                    {{ MSI::Dof_Type::TEMP }},
//                                                                    { fem::Property_Type::CONDUCTIVITY },
//                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//        tIWGUserDefinedInfo( 4 ).resize( 1 );
//        tIWGUserDefinedInfo( 4 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, { MSI::Dof_Type::TEMP },
//                                                                    {{ MSI::Dof_Type::TEMP }},
//                                                                    { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
//                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );
//        tIWGUserDefinedInfo( 5 ).resize( 1 );
//        tIWGUserDefinedInfo( 5 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN, { MSI::Dof_Type::TEMP },
//                                                                    {{ MSI::Dof_Type::TEMP }},
//                                                                    { fem::Property_Type::TEMP_NEUMANN },
//                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );
//
//        // create property user defined info
//        Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 6 );
//        tPropertyUserDefinedInfo( 0 ).resize( 1 );
//        tPropertyUserDefinedInfo( 0 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 1.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//        tPropertyUserDefinedInfo( 1 ).resize( 1 );
//        tPropertyUserDefinedInfo( 1 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 1.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//
//        tPropertyUserDefinedInfo( 2 ).resize( 1 );
//        tPropertyUserDefinedInfo( 2 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 1.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//        tPropertyUserDefinedInfo( 3 ).resize( 1 );
//        tPropertyUserDefinedInfo( 3 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 1.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//
//        tPropertyUserDefinedInfo( 4 ).resize( 2 );
//        tPropertyUserDefinedInfo( 4 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 1.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//        tPropertyUserDefinedInfo( 4 )( 1 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_DIRICHLET,
//                                                                             Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 5.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//        tPropertyUserDefinedInfo( 5 ).resize( 1 );
//        tPropertyUserDefinedInfo( 5 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_NEUMANN,
//                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                              {{{ 20.0 }}},
//                                                                              tConstValFunction,
//                                                                              Cell< fem::PropertyFunc >( 0 ) );
//
//
//        // create constitutive user defined info
//        Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 6 );
//        tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
//        tConstitutiveUserDefinedInfo( 0 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                      {{ MSI::Dof_Type::TEMP }},
//                                                                                      { fem::Property_Type::CONDUCTIVITY },
//                                                                                      tSpatialDimension );
//        tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
//        tConstitutiveUserDefinedInfo( 1 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                      {{ MSI::Dof_Type::TEMP }},
//                                                                                      { fem::Property_Type::CONDUCTIVITY },
//                                                                                      tSpatialDimension );
//        tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
//        tConstitutiveUserDefinedInfo( 2 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                      {{ MSI::Dof_Type::TEMP }},
//                                                                                      { fem::Property_Type::CONDUCTIVITY },
//                                                                                      tSpatialDimension );
//        tConstitutiveUserDefinedInfo( 3 ).resize( 1 );
//        tConstitutiveUserDefinedInfo( 4 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                      {{ MSI::Dof_Type::TEMP }},
//                                                                                      { fem::Property_Type::CONDUCTIVITY },
//                                                                                      tSpatialDimension );
//
//        // create a list of active block-sets
////        moris::Cell< moris_index >  tSetList = { 4, 5, 1, 3 };
//
//        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
//
//        // create a list of active block-sets
//         moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
//                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
//                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
//                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
//                                                  tEnrIntegMesh.get_side_set_index(tInterfaceSideSetName),
//                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p0")};
//
//        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
//                                                          fem::Element_Type::BULK,
//                                                          fem::Element_Type::BULK,
//                                                          fem::Element_Type::BULK,
//                                                          fem::Element_Type::SIDESET,
//                                                          fem::Element_Type::SIDESET };
//
//        // create model
//        mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex,
//                                              tIWGUserDefinedInfo,
//                                              tSetList, tSetTypeList,
//                                              tPropertyUserDefinedInfo,
//                                              tConstitutiveUserDefinedInfo );
//
//        moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 1: create linear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        dla::Solver_Factory  tSolFactory;
//        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//
//        dla::Linear_Solver tLinSolver;
//
//        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 2: create nonlinear solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        NLA::Nonlinear_Solver_Factory tNonlinFactory;
//        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
////        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
////        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
////        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
////        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//        NLA::Nonlinear_Solver tNonlinearSolver;
//
//        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 3: create time Solver and algorithm
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        tsa::Time_Solver_Factory tTimeSolverFactory;
//        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//        tsa::Time_Solver tTimeSolver;
//
//        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//        NLA::SOL_Warehouse tSolverWarehouse;
//
//        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//        tNonlinearSolver.set_dof_type_list( tDofTypes1 );
//        tTimeSolver.set_dof_type_list( tDofTypes1 );
//
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        // STEP 4: Solve and check
//        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//        tTimeSolver.solve();
//
//
//        // TODO: add gold solution data for this problem
//
//        // Write to Integration mesh for visualization
//        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//
////        moris::print_fancy(tIntegSol,"tIntegSol");
//
//        // add solution field to integration mesh
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);
//
//
////        Matrix<DDRMat> tFullSol;
////        tTimeSolver.get_full_solution(tFullSol);
////
////        print_fancy(tFullSol,"Full Solution");
//
//        // verify solution
////        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);
//        tModel->output_solution( "Circle" );
//        tField->put_scalar_values_on_field( tModel->get_mSolHMR() );
//        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_stk_bar_plane_interp_l2_b2.e" );
//
//
//        // output solution and meshes
////        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_stk_bar_hole_integ.e";
////        tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//        //    delete tInterpMesh1;
//        delete tModel;
//
//
//
//
//        delete tIntegMesh1;
//    }
}


TEST_CASE("2D XTK WITH HMR Error","[XTK_HMR_2D_Error]")
{
if(par_size()<=1)
{
//    std::string tFieldName = "Cylinder";
//
//     moris::uint tLagrangeMeshIndex = 0;
//     moris::uint tBSplineMeshIndex = 0;
//
//     hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//
//     tParameters.set( "number_of_elements_per_dimension", "5, 2" );
//     tParameters.set( "domain_dimensions", "5, 2" );
//     tParameters.set( "domain_offset", "-2.5, -1.0" );
//     tParameters.set( "domain_sidesets", "2, 4" );
//     tParameters.set( "lagrange_output_meshes", "0" );
//
//     tParameters.set( "lagrange_orders", "1" );
//     tParameters.set( "lagrange_pattern", "0" );
//     tParameters.set( "bspline_orders", "1" );
//     tParameters.set( "bspline_pattern", "0" );
//
//     tParameters.set( "lagrange_to_bspline", "0" );
//
//     tParameters.set( "truncate_bsplines", 1 );
//     tParameters.set( "refinement_buffer", 2 );
//     tParameters.set( "staircase_buffer", 2 );
//     tParameters.set( "initial_refinement", 2 );
//
//     tParameters.set( "use_multigrid", 0 );
//     tParameters.set( "severity_level", 2 );
//
//     hmr::HMR tHMR( tParameters );
//
//     // initial refinement
//     tHMR.perform_initial_refinement( 0 );
//
//     std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//     // create field
//     std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//     tField->evaluate_scalar_function( LvlSetCircle_1 );
//
////         for( uint k=0; k<2; ++k )
////         {
////             tHMR.flag_surface_elements_on_working_pattern( tField );
////             tHMR.perform_refinement_based_on_working_pattern( 0 );
////
////             tField->evaluate_scalar_function( CircleFunc );
////         }
//
//     tHMR.finalize();
//
//     tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );
//
//     std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//     xtk::Geom_Field tFieldAsGeom(tField);
//
//     moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//     size_t tModelDimension = 2;
//     xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//     xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
//     xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
//     tXTKModel.mVerbose = true;
//
//    //Specify decomposition Method and Cut Mesh ---------------------------------------
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//    tXTKModel.decompose(tDecompositionMethods);
//
//    tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);
//
//    // get meshes
//    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//    // output to exodus file ----------------------------------------------------------
////    xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();
////
////     // Declare the fields related to enrichment strategy in output options
////     Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
////
////    // output solution and meshes
////    xtk::Output_Options tOutputOptions;
////    tOutputOptions.mAddNodeSets = false;
////    tOutputOptions.mAddSideSets = true;
////    tOutputOptions.mAddClusters = false;
////
////    // add solution field to integration mesh
////    std::string tIntegSolFieldName = "solution";
////    tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
////    tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;
////
////    moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
////
////    tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh1);
////
////    std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_cut.e";
////    tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//    // end output -------------------------------------------------------
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
//
//    uint tSpatialDimension = 2;
//
//    // create IWG user defined info
//    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 4 );
//    tIWGUserDefinedInfo( 0 ).resize( 1 );
//    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::CONDUCTIVITY },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 1 ).resize( 1 );
//    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::CONDUCTIVITY },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//    tIWGUserDefinedInfo( 2 ).resize( 1 );
//    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
//                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );
//    tIWGUserDefinedInfo( 3 ).resize( 1 );
//    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_NEUMANN },
//                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );
//
//    // create property user defined info
//    Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 4 );
//    tPropertyUserDefinedInfo( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 0 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 1.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//    tPropertyUserDefinedInfo( 1 ).resize( 1 );
//    tPropertyUserDefinedInfo( 1 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 1.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//
//    tPropertyUserDefinedInfo( 2 ).resize( 2 );
//    tPropertyUserDefinedInfo( 2 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 1.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//    tPropertyUserDefinedInfo( 2 )( 1 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_DIRICHLET,
//                                                                         Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 5.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//    tPropertyUserDefinedInfo( 3 ).resize( 1 );
//    tPropertyUserDefinedInfo( 3 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_NEUMANN,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 20.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//
//
//    // create constitutive user defined info
//    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 4 );
//    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 0 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                  {{ MSI::Dof_Type::TEMP }},
//                                                                                  { fem::Property_Type::CONDUCTIVITY },
//                                                                                  tSpatialDimension );
//    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 1 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                  {{ MSI::Dof_Type::TEMP }},
//                                                                                  { fem::Property_Type::CONDUCTIVITY },
//                                                                                  tSpatialDimension );
//
//
//    // create a list of active block-sets
//    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
//
//    // create a list of active block-sets
//     moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
//                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_1_n_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1")};
//
//
//    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
//                                                      fem::Element_Type::BULK,
//                                                      fem::Element_Type::SIDESET,
//                                                      fem::Element_Type::SIDESET };
//
//    // create model
//    mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex,
//                                          tIWGUserDefinedInfo,
//                                          tSetList, tSetTypeList,
//                                          tPropertyUserDefinedInfo,
//                                          tConstitutiveUserDefinedInfo );
//
//    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 1: create linear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    dla::Solver_Factory  tSolFactory;
//    std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
//
//    tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//    tLinearSolverAlgorithm->set_param("AZ_max_iter") = 2;
//
//    dla::Linear_Solver tLinSolver;
//
//    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 2: create nonlinear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    NLA::Nonlinear_Solver_Factory tNonlinFactory;
//    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 1;
////        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
////        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
////        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//    NLA::Nonlinear_Solver tNonlinearSolver;
//
//    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 3: create time Solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 4: Solve and check
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    tTimeSolver.solve();
//
//    MSI::MSI_Solver_Interface * tSolver_Interface = tModel->get_solver_interface();
//
////    tSolver_Interface->write_solution_to_hdf5_file( "Exact_Sol.h5" );
//
//    char SolVector[100];
//    std::strcpy( SolVector, "Exact_Sol.h5" );
//
//    tSolver_Interface->get_exact_solution_from_hdf5_and_calculate_error( SolVector );
//
//    // output solution and meshes
//    xtk::Output_Options tOutputOptions;
//    tOutputOptions.mAddNodeSets = false;
//    tOutputOptions.mAddSideSets = true;
//    tOutputOptions.mAddClusters = false;
//
//    // add solution field to integration mesh
//    std::string tIntegSolFieldName = "solution";
//    tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//    moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//    // Write to Integration mesh for visualization
//    Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//    Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//
//    for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
//    {
//        moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//        tSTKIntegSol(i) = tIntegSol(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//    }
//
//    // crate field in integration mesh
//    moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("solution",EntityRank::NODE);
//    tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);
//
//    // add solution field to integration mesh
//    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);
//
////    Matrix<DDRMat> tFullSol;
////    tTimeSolver.get_full_solution(tFullSol);
//
//    std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d.e";
//    tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//    delete tIntegMesh1;
//
//    delete tModel;
//
}
}


TEST_CASE("2D XTK WITH HMR Error 2","[XTK_HMR_2D_ERROR_2]")
{
//if(par_size()<=1)
//{
//    std::string tFieldName = "Cylinder";
//
//     moris::uint tLagrangeMeshIndex = 0;
//     moris::uint tBSplineMeshIndex = 0;
//
//     hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
//
//     tParameters.set( "number_of_elements_per_dimension", "5, 2" );
//     tParameters.set( "domain_dimensions", "5, 2" );
//     tParameters.set( "domain_offset", "-2.5, -1.0" );
//     tParameters.set( "domain_sidesets", "2, 4" );
//     tParameters.set( "lagrange_output_meshes", "0" );
//
//     tParameters.set( "lagrange_orders", "1" );
//     tParameters.set( "lagrange_pattern", "0" );
//     tParameters.set( "bspline_orders", "1" );
//     tParameters.set( "bspline_pattern", "0" );
//
//     tParameters.set( "lagrange_to_bspline", "0" );
//
//     tParameters.set( "truncate_bsplines", 1 );
//     tParameters.set( "refinement_buffer", 2 );
//     tParameters.set( "staircase_buffer", 2 );
//     tParameters.set( "initial_refinement", 2 );
//
//     tParameters.set( "use_multigrid", 0 );
//     tParameters.set( "severity_level", 2 );
//
//     hmr::HMR tHMR( tParameters );
//
//     // initial refinement
//     tHMR.perform_initial_refinement( 0 );
//
//     std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//     // create field
//     std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//     tField->evaluate_scalar_function( LvlSetCircle_1 );
//
////         for( uint k=0; k<4; ++k )
////         {
////             tHMR.flag_surface_elements_on_working_pattern( tField );
////             tHMR.perform_refinement_based_on_working_pattern( 0 );
////
////             tField->evaluate_scalar_function( LvlSetCircle_1 );
////         }
//
//     tHMR.finalize();
//
//     tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );
//
//     std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//     xtk::Geom_Field tFieldAsGeom(tField);
//
//     moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//     size_t tModelDimension = 2;
//     xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//     xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
//     xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
//     tXTKModel.mVerbose = true;
//
//    //Specify decomposition Method and Cut Mesh ---------------------------------------
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//    tXTKModel.decompose(tDecompositionMethods);
//
//    tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);
//
//    // get meshes
//    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
//
//    uint tSpatialDimension = 2;
//
//    // create IWG user defined info
//    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 4 );
//    tIWGUserDefinedInfo( 0 ).resize( 1 );
//    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::CONDUCTIVITY },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//    tIWGUserDefinedInfo( 1 ).resize( 1 );
//    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::CONDUCTIVITY },
//                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
//
//    tIWGUserDefinedInfo( 2 ).resize( 1 );
//    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
//                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );
//    tIWGUserDefinedInfo( 3 ).resize( 1 );
//    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN, tSpatialDimension, { MSI::Dof_Type::TEMP },
//                                                                {{ MSI::Dof_Type::TEMP }},
//                                                                { fem::Property_Type::TEMP_NEUMANN },
//                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );
//
//    // create property user defined info
//    Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 4 );
//    tPropertyUserDefinedInfo( 0 ).resize( 1 );
//    tPropertyUserDefinedInfo( 0 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 1.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//    tPropertyUserDefinedInfo( 1 ).resize( 1 );
//    tPropertyUserDefinedInfo( 1 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 1.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//
//    tPropertyUserDefinedInfo( 2 ).resize( 2 );
//    tPropertyUserDefinedInfo( 2 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 1.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//    tPropertyUserDefinedInfo( 2 )( 1 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_DIRICHLET,
//                                                                         Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 5.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//    tPropertyUserDefinedInfo( 3 ).resize( 1 );
//    tPropertyUserDefinedInfo( 3 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_NEUMANN,
//                                                                          Cell< Cell< MSI::Dof_Type > >( 0 ),
//                                                                          {{{ 20.0 }}},
//                                                                          tConstValFunction,
//                                                                          Cell< fem::PropertyFunc >( 0 ) );
//
//
//    // create constitutive user defined info
//    Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 4 );
//    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 0 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                  {{ MSI::Dof_Type::TEMP }},
//                                                                                  { fem::Property_Type::CONDUCTIVITY },
//                                                                                  tSpatialDimension );
//    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
//    tConstitutiveUserDefinedInfo( 1 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
//                                                                                  {{ MSI::Dof_Type::TEMP }},
//                                                                                  { fem::Property_Type::CONDUCTIVITY },
//                                                                                  tSpatialDimension );
//
//
//    // create a list of active block-sets
//    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
//
//    // create a list of active block-sets
//     moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
//                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_1_n_p1"),
//                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1")};
//
//
//    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
//                                                      fem::Element_Type::BULK,
//                                                      fem::Element_Type::SIDESET,
//                                                      fem::Element_Type::SIDESET };
//
//    // create model
//    mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex,
//                                          tIWGUserDefinedInfo,
//                                          tSetList, tSetTypeList,
//                                          tPropertyUserDefinedInfo,
//                                          tConstitutiveUserDefinedInfo );
//
//    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 1: create linear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    dla::Solver_Factory  tSolFactory;
//            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );
//
////            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPPREONLY );
//            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPRICHARDSON );
////            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPGMRES );
////            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCMG );
////            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCSOR );
//            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCJACOBI );
////            tLinearSolverAlgorithm->set_param("ILUFill")  = 3;
//            tLinearSolverAlgorithm->set_param("KSPMaxits")  = 1000;
//
//    dla::Linear_Solver tLinSolver;
//
//    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 2: create nonlinear solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );
//
//    NLA::Nonlinear_Solver_Factory tNonlinFactory;
//    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 1;
////        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
////        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
////        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
//
//    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//    NLA::Nonlinear_Solver tNonlinearSolver;
//
//    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//    tNonlinearSolver.solve( tNonlinearProblem );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 3: create time Solver and algorithm
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
////    tsa::Time_Solver_Factory tTimeSolverFactory;
////    std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
////
////    tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
////
////    tsa::Time_Solver tTimeSolver;
////
////    tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
////
////    NLA::SOL_Warehouse tSolverWarehouse;
////
////    tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
////
////    tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
////    tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
////
////    tNonlinearSolver.set_dof_type_list( tDofTypes1 );
////    tTimeSolver.set_dof_type_list( tDofTypes1 );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // STEP 4: Solve and check
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
////    tTimeSolver.solve();
//
//    MSI::MSI_Solver_Interface * tSolver_Interface = tModel->get_solver_interface();
//
////    tSolver_Interface->write_solution_to_hdf5_file( "Exact_Sol_petsc.h5" );
//
//    char SolVector[100];
//    std::strcpy( SolVector, "Exact_Sol_petsc.h5" );
//
//    tSolver_Interface->get_exact_solution_from_hdf5_and_calculate_error( SolVector );
//
//    // output solution and meshes
//    xtk::Output_Options tOutputOptions;
//    tOutputOptions.mAddNodeSets = false;
//    tOutputOptions.mAddSideSets = true;
//    tOutputOptions.mAddClusters = false;
//
//    // add solution field to integration mesh
//    std::string tIntegSolFieldName = "solution";
//    tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//    moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//    // Write to Integration mesh for visualization
//    Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
//
//    Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//
//    for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
//    {
//        moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//        tSTKIntegSol(i) = tIntegSol(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//    }
//
//    // crate field in integration mesh
//    moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("solution",EntityRank::NODE);
//    tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);
//
//    // add solution field to integration mesh
//    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);
//
////    Matrix<DDRMat> tFullSol;
////    tTimeSolver.get_full_solution(tFullSol);
//
//    std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d_petsc_no_hole_1000.e";
//    tIntegMesh1->create_output_mesh(tMeshOutputFile);
//
//    delete tIntegMesh1;
//
//    delete tModel;
//
//}
}



TEST_CASE("2D XTK WITH HMR Error Analytical","[XTK_HMR_2D_ERROR_Analytical]")
{
if(par_size()<=1)
{
    std::string tFieldName = "Cylinder";

     moris::uint tLagrangeMeshIndex = 0;
     moris::uint tBSplineMeshIndex = 0;

     hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

     tParameters.set( "number_of_elements_per_dimension", "8, 8");
     tParameters.set( "domain_dimensions", "4, 4" );
     tParameters.set( "domain_offset", "-2.0, -2.0" );
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
     tParameters.set( "initial_refinement", 2 );

     tParameters.set( "use_multigrid", 0 );
     tParameters.set( "severity_level", 2 );

     hmr::HMR tHMR( tParameters );

     // initial refinement
     tHMR.perform_initial_refinement( 0 );

     std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

     // create field
     std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

     tField->evaluate_scalar_function( LvlSetCircle_1 );

         for( uint k=0; k<4; ++k )
         {
             tHMR.flag_surface_elements_on_working_pattern( tField );
             tHMR.perform_refinement_based_on_working_pattern( 0 );

             tField->evaluate_scalar_function( LvlSetCircle_1 );
         }

     tHMR.finalize();

     tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

     std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

     xtk::Geom_Field tFieldAsGeom(tField);

     moris::Cell<xtk::Geometry*> tGeometryVector = {&tFieldAsGeom};

     size_t tModelDimension = 2;
     xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
     xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
     xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
     tXTKModel.mVerbose = true;

    //Specify decomposition Method and Cut Mesh ---------------------------------------
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
    tXTKModel.decompose(tDecompositionMethods);

    tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

    // get meshes
    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

    uint tSpatialDimension = 2;

    // create IWG user defined info
    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 7 );
    tIWGUserDefinedInfo( 0 ).resize( 1 );
    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                moris::Cell< fem::Property_Type >( 0 ),
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 1 ).resize( 1 );
    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                moris::Cell< fem::Property_Type >( 0 ),
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );

    tIWGUserDefinedInfo( 2 ).resize( 1 );
    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_DIRICHLET },
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 3 ).resize( 1 );
    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_DIRICHLET },
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 4 ).resize( 1 );
    tIWGUserDefinedInfo( 4 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_DIRICHLET },
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 5 ).resize( 1 );
    tIWGUserDefinedInfo( 5 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_DIRICHLET },
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 6 ).resize( 1 );
    tIWGUserDefinedInfo( 6 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN,  { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_NEUMANN },
                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );

    // create property user defined info
    fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
                                                   Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                   {{{ 1.0 }}},
                                                   tConstValFunction,
                                                   Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
                                                    Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                    {{{ 0 }}},
                                                    exactTempFunc,
                                                    Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tTempNeumann( fem::Property_Type::TEMP_NEUMANN,
                                                  Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                  {{{ 1.0 }}},
                                                  tConstValFunction,
                                                  Cell< fem::PropertyFunc >( 0 ) );

    // create property user defined info
    Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 7 );
    tPropertyUserDefinedInfo( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 0 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 1 ).resize( 1 );
    tPropertyUserDefinedInfo( 1 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 2 ).resize( 1 );
    tPropertyUserDefinedInfo( 2 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tTempDirichlet;
    tPropertyUserDefinedInfo( 3 ).resize( 1 );
    tPropertyUserDefinedInfo( 3 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 3 )( 0 )( 1 ) = tTempDirichlet;
    tPropertyUserDefinedInfo( 4 ).resize( 1 );
    tPropertyUserDefinedInfo( 4 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 4 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 4 )( 0 )( 1 ) = tTempDirichlet;
    tPropertyUserDefinedInfo( 5 ).resize( 1 );
    tPropertyUserDefinedInfo( 5 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 5 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 5 )( 0 )( 1 ) = tTempDirichlet;
    tPropertyUserDefinedInfo( 6 ).resize( 1 );
    tPropertyUserDefinedInfo( 6 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 6 )( 0 )( 0 ) = tTempNeumann;


    fem::Constitutive_User_Defined_Info tDiffLinIso( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                     {{ MSI::Dof_Type::TEMP }},
                                                     { fem::Property_Type::CONDUCTIVITY } );
    // create constitutive user defined info
    Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 7 );
    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 0 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 0 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 1 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 1 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 2 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 2 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 3 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 3 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 3 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 4 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 4 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 4 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 5 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 5 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 5 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 6 ).resize( 1 );




    // create a list of active block-sets
    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 1, 0 );

    // create a list of active block-sets
     moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_1_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_3_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1"),
                                              tEnrIntegMesh.get_side_set_index( tInterfaceSideSetName )
                                              };


    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                      fem::Element_Type::BULK,
                                                      fem::Element_Type::SIDESET,
                                                      fem::Element_Type::SIDESET,
                                                      fem::Element_Type::SIDESET,
                                                      fem::Element_Type::SIDESET,
                                                      fem::Element_Type::SIDESET
     };

    // create model
    mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex,
                                          tSetList,
                                          tSetTypeList,
                                          tIWGUserDefinedInfo,
                                          tPropertyUserDefinedInfo,
                                          tConstitutiveUserDefinedInfo );

    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 1: create linear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );

//            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPPREONLY );
//            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPRICHARDSON );
            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPGMRES );
//            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCMG );
//            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCSOR );
            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCILU );
            tLinearSolverAlgorithm->set_param("ILUFill")  = 3;
            tLinearSolverAlgorithm->set_param("KSPMaxits")  = 1000;

    dla::Linear_Solver tLinSolver;

    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 2: create nonlinear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );

    NLA::Nonlinear_Solver_Factory tNonlinFactory;
    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

    NLA::Nonlinear_Solver tNonlinearSolver;

    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

    tNonlinearSolver.solve( tNonlinearProblem );

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

    MSI::MSI_Solver_Interface * tSolver_Interface = tModel->get_solver_interface();

    tSolver_Interface->write_solution_to_hdf5_file( "Exact_Sol_petsc.h5" );

//    char SolVector[100];
//    std::strcpy( SolVector, "Exact_Sol_petsc.h5" );
//
//    tSolver_Interface->get_exact_solution_from_hdf5_and_calculate_error( SolVector );

    // output solution and meshes
    xtk::Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = true;
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
    moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("solution",EntityRank::NODE);
    tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

    // add solution field to integration mesh
    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);

//    Matrix<DDRMat> tFullSol;
//    tTimeSolver.get_full_solution(tFullSol);

    std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d_pipe.e";
    tIntegMesh1->create_output_mesh(tMeshOutputFile);

    delete tIntegMesh1;

    delete tModel;

}
}


TEST_CASE("2D XTK WITH HMR Residual","[XTK_HMR_2D_Residual]")
{
if(par_size()<=1)
{
    std::string tFieldName = "Cylinder";

//    std::string tPrefix = std::getenv("MORISROOT");
    std::string tMeshFileName = "/home/schmidt/Work/Star_residual/Star_1x1.g";
    std::cout<<"Mesh input name = "<<tMeshFileName<<std::endl;

    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
    std::string tFieldName1 = "Temp_Field";
    tNodeField1.set_field_name( tFieldName1 );
    tNodeField1.set_field_entity_rank( EntityRank::NODE );

    // Initialize field information container
    moris::mtk::MtkFieldsInfo tFieldsInfo;

    // Place the node field into the field info container
    add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

    // Declare some supplementary fields
    mtk::MtkMeshData tMeshData;
    tMeshData.FieldsInfo = &tFieldsInfo;


    // construct the mesh data
    mtk::Interpolation_Mesh* tInterpMesh = mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tMeshData );

    moris::Matrix<moris::DDRMat> tCenters = {{ 1.0,1.0,3.1 }};
    moris::Matrix<moris::DDRMat> tNormals = {{ 0.0,0.0,1.0 }};
    xtk::Star tStar;

    xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    xtk::Geometry_Engine tGeometryEngine(tStar,tPhaseTable);

     xtk::Model tXTKModel(2, tInterpMesh, tGeometryEngine);
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

    uint tSpatialDimension = 2;

    // create IWG user defined info
    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 4 );
    tIWGUserDefinedInfo( 0 ).resize( 1 );
    tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                moris::Cell< fem::Property_Type >( 0 ),
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 1 ).resize( 1 );
    tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                moris::Cell< fem::Property_Type >( 0 ),
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );

    tIWGUserDefinedInfo( 2 ).resize( 1 );
    tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_DIRICHLET },
                                                                { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tIWGUserDefinedInfo( 3 ).resize( 1 );
    tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN,  { MSI::Dof_Type::TEMP },
                                                                {{ MSI::Dof_Type::TEMP }},
                                                                { fem::Property_Type::TEMP_NEUMANN },
                                                                moris::Cell< fem::Constitutive_Type >( 0 ) );

    // create property user defined info
    fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
                                                   Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                   {{{ 1.0 }}},
                                                   tConstValFunction,
                                                   Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
                                                    Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                    {{{ 5.0 }}},
                                                    tConstValFunction,
                                                    Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tTempNeumann( fem::Property_Type::TEMP_NEUMANN,
                                                  Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                  {{{ 20.0 }}},
                                                  tConstValFunction,
                                                  Cell< fem::PropertyFunc >( 0 ) );

    // create property user defined info
    Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 4 );
    tPropertyUserDefinedInfo( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 0 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 1 ).resize( 1 );
    tPropertyUserDefinedInfo( 1 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 2 ).resize( 1 );
    tPropertyUserDefinedInfo( 2 )( 0 ).resize( 2 );
    tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tConductivity;
    tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tTempDirichlet;
    tPropertyUserDefinedInfo( 3 ).resize( 1 );
    tPropertyUserDefinedInfo( 3 )( 0 ).resize( 1 );
    tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tTempNeumann;

    fem::Constitutive_User_Defined_Info tDiffLinIso( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                     {{ MSI::Dof_Type::TEMP }},
                                                     { fem::Property_Type::CONDUCTIVITY } );

    // create constitutive user defined info
    Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 4 );
    tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 0 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 0 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 1 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 1 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 2 )( 0 ).resize( 1 );
    tConstitutiveUserDefinedInfo( 2 )( 0 )( 0 ) = tDiffLinIso;
    tConstitutiveUserDefinedInfo( 3 ).resize( 1 );


    // create a list of active block-sets
    std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );

    // create a list of active block-sets
     moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                              tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_1_n_p1"),
                                              tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1")};


    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                      fem::Element_Type::BULK,
                                                      fem::Element_Type::SIDESET,
                                                      fem::Element_Type::SIDESET };

    uint tBSplineMeshIndex = 0;
    // create model
    mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                          tBSplineMeshIndex,
                                          tSetList, tSetTypeList,
                                          tIWGUserDefinedInfo,
                                          tPropertyUserDefinedInfo,
                                          tConstitutiveUserDefinedInfo );

    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 1: create linear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );

//            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPPREONLY );
            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPRICHARDSON );
//            tLinearSolverAlgorithm->set_param("KSPType") = std::string( KSPGMRES );
//            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCMG );
//            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCSOR );
            tLinearSolverAlgorithm->set_param("PCType")  = std::string( PCJACOBI );
//            tLinearSolverAlgorithm->set_param("ILUFill")  = 3;
            tLinearSolverAlgorithm->set_param("KSPMaxits")  = 1000;

    dla::Linear_Solver tLinSolver;

    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 2: create nonlinear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );

    NLA::Nonlinear_Solver_Factory tNonlinFactory;
    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 1;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

    NLA::Nonlinear_Solver tNonlinearSolver;

    tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

    tNonlinearSolver.solve( tNonlinearProblem );

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

    MSI::MSI_Solver_Interface * tSolver_Interface = tModel->get_solver_interface();

//    tSolver_Interface->write_solution_to_hdf5_file( "Exact_Sol_petsc.h5" );

//    char SolVector[100];
//    std::strcpy( SolVector, "Exact_Sol_petsc.h5" );
//
//    tSolver_Interface->get_exact_solution_from_hdf5_and_calculate_error( SolVector );

    tSolver_Interface->get_residual_vecor_for_output(SolVector);

    // output solution and meshes
    xtk::Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = true;
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
    moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("solution",EntityRank::NODE);
    tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

    // add solution field to integration mesh
    tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);

//    Matrix<DDRMat> tFullSol;
//    tTimeSolver.get_full_solution(tFullSol);

    std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_hole_integ_2d_petsc_no_hole_1000.e";
    tIntegMesh1->create_output_mesh(tMeshOutputFile);

    delete tIntegMesh1;

    delete tModel;

}
}


}

