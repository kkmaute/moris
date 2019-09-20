/*
 * UT_MDL_XTK_HMR_2D.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: schmidt
 */

#include "catch.hpp"

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

moris::real LvlSetLin(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tOffset = -0.6;

    return    aPoint(0) - 0.317 * aPoint(1) - tOffset;
}

Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                    moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator,
                                    fem::Geometry_Interpolator             * aGeometryInterpolator )
{
    return aCoeff( 0 );
}

TEST_CASE("2D XTK WITH HMR 2D","[XTK_HMR_2D]")
{
    if(par_size()<=1)
    {
        std::string tFieldName = "Cylinder";

         moris::uint tLagrangeMeshIndex = 0;
         moris::uint tBSplineMeshIndex = 0;

         hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

         tParameters.set( "number_of_elements_per_dimension", "5, 2" );
         tParameters.set( "domain_dimensions", "5, 2" );
         tParameters.set( "domain_offset", "-2.5, -1.0" );
         tParameters.set( "domain_sidesets", "2, 4" );
         tParameters.set( "lagrange_output_meshes", "0" );

         tParameters.set( "lagrange_orders", "1" );
         tParameters.set( "lagrange_pattern", "0" );
         tParameters.set( "bspline_orders", "1" );
         tParameters.set( "bspline_pattern", "0" );

         tParameters.set( "lagrange_to_bspline", "0" );

         tParameters.set( "truncate_bsplines", 1 );
         tParameters.set( "refinement_buffer", 2 );
         tParameters.set( "staircase_buffer", 2 );
         tParameters.set( "initial_refinement", 2 );

         tParameters.set( "use_multigrid", 0 );
         tParameters.set( "severity_level", 2 );

         hmr::HMR tHMR( tParameters );

         // initial refinement
         tHMR.perform_initial_refinement( 0 );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

         // create field
         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

         tField->evaluate_scalar_function( LvlSetLin );

//         for( uint k=0; k<2; ++k )
//         {
//             tHMR.flag_surface_elements_on_working_pattern( tField );
//             tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//             tField->evaluate_scalar_function( CircleFunc );
//         }

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

        // output to exodus file ----------------------------------------------------------
        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

         // Declare the fields related to enrichment strategy in output options
         Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh1);

        std::string tMeshOutputFile ="./xtk_exo/xtk_hmr_2d_cut.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        // end output -------------------------------------------------------

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        uint tSpatialDimension = 2;

        // create IWG user defined info
        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 6 );
        tIWGUserDefinedInfo( 0 ).resize( 1 );
        tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 1 ).resize( 1 );
        tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );

        tIWGUserDefinedInfo( 2 ).resize( 1 );
        tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 3 ).resize( 1 );
        tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, tSpatialDimension, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );

        tIWGUserDefinedInfo( 4 ).resize( 1 );
        tIWGUserDefinedInfo( 4 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, tSpatialDimension, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );
        tIWGUserDefinedInfo( 5 ).resize( 1 );
        tIWGUserDefinedInfo( 5 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN, tSpatialDimension, { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::TEMP_NEUMANN },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );

        // create property user defined info
        Cell< Cell< fem::Property_User_Defined_Info > > tPropertyUserDefinedInfo( 6 );
        tPropertyUserDefinedInfo( 0 ).resize( 1 );
        tPropertyUserDefinedInfo( 0 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 1 ).resize( 1 );
        tPropertyUserDefinedInfo( 1 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );

        tPropertyUserDefinedInfo( 2 ).resize( 1 );
        tPropertyUserDefinedInfo( 2 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 3 ).resize( 1 );
        tPropertyUserDefinedInfo( 3 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );

        tPropertyUserDefinedInfo( 4 ).resize( 2 );
        tPropertyUserDefinedInfo( 4 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::CONDUCTIVITY,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 1.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 4 )( 1 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_DIRICHLET,
                                                                             Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 5.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );
        tPropertyUserDefinedInfo( 5 ).resize( 1 );
        tPropertyUserDefinedInfo( 5 )( 0 ) = fem::Property_User_Defined_Info( fem::Property_Type::TEMP_NEUMANN,
                                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                                              {{{ 20.0 }}},
                                                                              tConstValFunction,
                                                                              Cell< fem::PropertyFunc >( 0 ) );


        // create constitutive user defined info
        Cell< Cell< fem::Constitutive_User_Defined_Info > > tConstitutiveUserDefinedInfo( 6 );
        tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 0 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      tSpatialDimension );
        tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 1 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      tSpatialDimension );
        tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 2 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      tSpatialDimension );
        tConstitutiveUserDefinedInfo( 3 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 4 )( 0 ) = fem::Constitutive_User_Defined_Info( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                                      { fem::Property_Type::CONDUCTIVITY },
                                                                                      tSpatialDimension );

        // create a list of active block-sets
//        moris::Cell< moris_index >  tSetList = { 4, 5, 1, 3 };

        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );

        // create a list of active block-sets
         moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                                  tEnrIntegMesh.get_side_set_index(tInterfaceSideSetName),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p0")};

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET };

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager, tBSplineMeshIndex,
                                              tIWGUserDefinedInfo,
                                              tSetList, tSetTypeList,
                                              tPropertyUserDefinedInfo,
                                              tConstitutiveUserDefinedInfo );

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

//        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
//        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
//        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
//        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

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


        // TODO: add gold solution data for this problem

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );


//        moris::print_fancy(tIntegSol,"tIntegSol");

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tIntegSol);


//        Matrix<DDRMat> tFullSol;
//        tTimeSolver.get_full_solution(tFullSol);
//
//        print_fancy(tFullSol,"Full Solution");

        // verify solution
//        CHECK(norm(tSolution11 - tGoldSolution)<1e-08);
        tModel->output_solution( "Circle" );
        tField->put_scalar_values_on_field( tModel->get_mSolHMR() );
        tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_stk_bar_plane_interp_l2_b2.e" );


        // output solution and meshes
//        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_stk_bar_hole_integ.e";
//        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        //    delete tInterpMesh1;
        delete tModel;




        delete tIntegMesh1;
    }
}






}

