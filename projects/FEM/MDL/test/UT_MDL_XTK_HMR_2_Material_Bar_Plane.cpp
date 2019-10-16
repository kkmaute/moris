/*
 * UT_MDL_XTK_HMR_Multi_Material_Bar_Plane.cpp
 *
 *  Created on: Oct 4, 2019
 *      Author: doble
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
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Constitutive_User_Defined_Info.hpp"      //FEM/INT/src

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

moris::real
Plane2MatMDL(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 0.11;
    moris::real mYC = 0.11;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aPoint(0)-mXC) + mNy*(aPoint(1)-mYC));
}

moris::real
Circle2MatMDL(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXCenter = 0.01;
    moris::real mYCenter = 0.01;
    moris::real mRadius = 0.47334;


    return  (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                    + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                    - (mRadius * mRadius);
}


Matrix< DDRMat > tConstValFunction2MatMDL( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                           moris::Cell< fem::Field_Interpolator* > & aDofFI,
                                           moris::Cell< fem::Field_Interpolator* > & aDvFI,
                                           fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

TEST_CASE("XTK HMR 2 Material Bar Intersected By Plane","[XTK_HMR_PLANE_BAR_2D]")
{


    if(par_size() == 1)
    {
        std::string tFieldName = "Geometry";

        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {11}, {4}} );
        tParameters.set_domain_dimensions({ {6}, {2} });
        tParameters.set_domain_offset({ {-3.0}, {-1.0} });
        tParameters.set_bspline_truncation( true );

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_side_sets({{1},{2},{3},{4} });
        tParameters.set_max_refinement_level( 2 );
        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tPlaneField  = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tPlaneField->evaluate_scalar_function( Plane2MatMDL);

        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tPlaneField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
            tPlaneField->evaluate_scalar_function( Plane2MatMDL);
        }

        tPlaneField->evaluate_scalar_function( Plane2MatMDL);
        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./xtk_exo/xtk_hmr_2d_enr_ip2.e" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        xtk::Geom_Field tPlaneFieldAsGeom(tPlaneField);

        moris::Cell<xtk::Geometry*> tGeometryVector = {&tPlaneFieldAsGeom};

        size_t tModelDimension = 2;
        xtk::Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
        xtk::Model tXTKModel(tModelDimension,tInterpMesh.get(),tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        tEnrIntegMesh.print();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create IWG user defined info
        Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 7 );
        tIWGUserDefinedInfo( 0 ).resize( 1 );
        tIWGUserDefinedInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    {fem::Property_Type::TEMP_LOAD },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 1 ).resize( 1 );
        tIWGUserDefinedInfo( 1 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    {fem::Property_Type::TEMP_LOAD },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 2 ).resize( 1 );
        tIWGUserDefinedInfo( 2 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    {fem::Property_Type::TEMP_LOAD },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 3 ).resize( 1 );
        tIWGUserDefinedInfo( 3 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    {fem::Property_Type::TEMP_LOAD },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 4 ).resize( 1 );
        tIWGUserDefinedInfo( 4 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::TEMP_DIRICHLET },
                                                                    { fem::Constitutive_Type::DIFF_LIN_ISO } );
        tIWGUserDefinedInfo( 5 ).resize( 1 );
        tIWGUserDefinedInfo( 5 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    { fem::Property_Type::TEMP_NEUMANN },
                                                                    moris::Cell< fem::Constitutive_Type >( 0 ) );

        tIWGUserDefinedInfo( 6 ).resize( 1 );
        tIWGUserDefinedInfo( 6 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_INTERFACE,
                                                                    { MSI::Dof_Type::TEMP },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    Cell< fem::Property_Type >( 0 ),
                                                                    {fem::Constitutive_Type::DIFF_LIN_ISO },
                                                                    {{ MSI::Dof_Type::TEMP }},
                                                                    Cell< fem::Property_Type >( 0 ),
                                                                    {fem::Constitutive_Type::DIFF_LIN_ISO } );

        // create the property user defined infos
        fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
                                                       Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                       {{{ 1.0 }}},
                                                       tConstValFunction2MatMDL,
                                                       Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tConductivity2( fem::Property_Type::CONDUCTIVITY,
                                                       Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                       {{{ 5.0 }}},
                                                       tConstValFunction2MatMDL,
                                                       Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
                                                        Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                        {{{ 5.0 }}},
                                                        tConstValFunction2MatMDL,
                                                        Cell< fem::PropertyFunc >( 0 ) );
        fem::Property_User_Defined_Info tTempNeumann( fem::Property_Type::TEMP_NEUMANN,
                                                      Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                      {{{ 20.0 }}},
                                                      tConstValFunction2MatMDL,
                                                      Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tTempLoad1( fem::Property_Type::TEMP_LOAD,
                                                              Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                              {{{ 100.0 }}},
                                                              tConstValFunction2MatMDL,
                                                              Cell< fem::PropertyFunc >( 0 ) );

        fem::Property_User_Defined_Info tTempLoad2( fem::Property_Type::TEMP_LOAD,
                                                      Cell< Cell< MSI::Dof_Type > >( 0 ),
                                                      {{{ 100.0 }}},
                                                      tConstValFunction2MatMDL,
                                                      Cell< fem::PropertyFunc >( 0 ) );

        // create property user defined info
        Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 7 );
        tPropertyUserDefinedInfo( 0 ).resize( 1 );
        tPropertyUserDefinedInfo( 0 )( 0 ).resize( 2 );
        tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tConductivity;
        tPropertyUserDefinedInfo( 0 )( 0 )( 1 ) = tTempLoad1;

        tPropertyUserDefinedInfo( 1 ).resize( 1 );
        tPropertyUserDefinedInfo( 1 )( 0 ).resize( 2 );
        tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tConductivity;
        tPropertyUserDefinedInfo( 1 )( 0 )( 1 ) = tTempLoad1;

        tPropertyUserDefinedInfo( 2 ).resize( 1 );
        tPropertyUserDefinedInfo( 2 )( 0 ).resize( 2 );
        tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tConductivity2;
        tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tTempLoad2;

        tPropertyUserDefinedInfo( 3 ).resize( 1 );
        tPropertyUserDefinedInfo( 3 )( 0 ).resize( 2 );
        tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tConductivity2;
        tPropertyUserDefinedInfo( 3 )( 0 )( 1 ) = tTempLoad2;

        tPropertyUserDefinedInfo( 4 ).resize( 1 );
        tPropertyUserDefinedInfo( 4 )( 0 ).resize( 2 );
        tPropertyUserDefinedInfo( 4 )( 0 )( 0 ) = tConductivity2;
        tPropertyUserDefinedInfo( 4 )( 0 )( 1 ) = tTempDirichlet;

        tPropertyUserDefinedInfo( 5 ).resize( 1 );
        tPropertyUserDefinedInfo( 5 )( 0 ).resize( 1 );
        tPropertyUserDefinedInfo( 5 )( 0 )( 0 ) = tTempNeumann;

        tPropertyUserDefinedInfo( 6 ).resize( 2 );
        tPropertyUserDefinedInfo( 6 )( 0 ).resize( 1 );
        tPropertyUserDefinedInfo( 6 )( 0 )( 0 ) = tConductivity;
        tPropertyUserDefinedInfo( 6 )( 1 ).resize( 1 );
        tPropertyUserDefinedInfo( 6 )( 1 )( 0 ) = tConductivity2;

        // create constitutive user defined info
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

        tConstitutiveUserDefinedInfo( 6 ).resize( 2 );
        tConstitutiveUserDefinedInfo( 6 )( 0 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 6 )( 0 )( 0 ) = tDiffLinIso;
        tConstitutiveUserDefinedInfo( 6 )( 1 ).resize( 1 );
        tConstitutiveUserDefinedInfo( 6 )( 1 )( 0 ) = tDiffLinIso;

        // create a list of active block-sets
        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
        std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
        moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1"),
                                                  tEnrIntegMesh.get_side_set_index("SideSet_4_n_p0"),
                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName)};

        moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::BULK,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::SIDESET,
                                                          fem::Element_Type::DOUBLE_SIDESET};


        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              tBSplineMeshIndex,
                                              tSetList,
                                              tSetTypeList,
                                              tIWGUserDefinedInfo,
                                              tPropertyUserDefinedInfo,
                                              tConstitutiveUserDefinedInfo,
                                              0,
                                              false);

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

        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

        // Declare the fields related to enrichment strategy in output options
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
            moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);

            tSTKIntegSol(i) = tIntegSol(tMyIndex);
        }

        // crate field in integration mesh
        moris::moris_index tFieldIndex = tEnrIntegMesh.create_field("Solution",EntityRank::NODE);
        tEnrIntegMesh.add_field_data(tFieldIndex,EntityRank::NODE,tSTKIntegSol);

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldName,EntityRank::NODE,tSTKIntegSol);

        std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_plane_2_mat_integ_2d.e";
        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        //    delete tInterpMesh1;
        delete tModel;
        delete tIntegMesh1;
    }
}
