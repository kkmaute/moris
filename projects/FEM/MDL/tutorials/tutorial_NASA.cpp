/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * tutorial_NASA.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

// dynamik linker function
#include "dlfcn.h"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_Geom_Field.hpp"
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

#include "cl_Library_IO.hpp"

#include "banner.hpp" // COR/src

#include "fn_norm.hpp"

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Logger.hpp" // MRS/IOS/src

#include "fn_properties_and_constitutive_model.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int
main( int    argc,
      char * argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm.initialize(&argc, &argv);

    // Severity level 0 - all outputs
    gLogger.initialize( 2 );

//    moris::print_banner( argc, argv );

    Library_IO tLibrary( argv[ 1 ] );

    // load user defined function
    MORIS_USER_FUNCTION user_input = tLibrary.load_function( "Input_Parameters" );

    Cell< real > tInput;

    user_input( tInput );

    uint tDim               = static_cast<uint>( tInput(0) );
    uint tBsplineOrder      = static_cast<uint>( tInput(1) );
    uint tRefinementLevels  = static_cast<uint>( tInput(2) );
    real tKappa_1           = tInput(3);
    real tKappa_2           = tInput(4);

    real tDirichlet         = tInput(5);
    real tNeumann           = tInput(6);

    real tLoad_1            = tInput(7);
    real tLoad_2            = tInput(8);

    std::string tFieldName = "Geometry";

    moris::uint tLagrangeMeshIndex = 0;
    moris::uint tBSplineMeshIndex = 0;

    moris::Matrix<DDLUMat> tElementsPerDimenson;
    moris::Matrix<DDRMat> tDimensons;
    moris::Matrix<DDRMat> tOffset;
    moris::Matrix<DDUMat> tSideSets;

    Cell<enum Subdivision_Method> tDecompositionMethods(2);

    if ( tDim ==2 )
    {
        tElementsPerDimenson = { {22}, {8} };
        tDimensons = { {6}, {2} };
        tOffset = { {-3.0}, {-1.0} };
        tSideSets = { {1},{2},{3},{4} };

        tDecompositionMethods( 0 ) = Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4;
        tDecompositionMethods( 1 ) = Subdivision_Method::C_TRI3;
    }
    else if ( tDim ==3 )
    {
        tElementsPerDimenson = { {22}, {8}, {8} };
        tDimensons = { {6}, {2}, {2} };
        tOffset = { {-3.0}, {-1.0},{-1.0} };
        tSideSets = {{1},{2},{3},{4},{5},{6}};

        tDecompositionMethods( 0 ) = Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8;
        tDecompositionMethods( 1 ) = Subdivision_Method::C_HIERARCHY_TET4;
    }

    moris::uint tLagrangeOrder = tBsplineOrder;

    moris::hmr::Parameters tParameters;
    tParameters.set_number_of_elements_per_dimension( tElementsPerDimenson );
    tParameters.set_domain_dimensions( tDimensons );
    tParameters.set_domain_offset( tOffset );
    tParameters.set_bspline_truncation( true );
    tParameters.set_output_meshes( {{ {0} }} );
    tParameters.set_lagrange_orders  ( { {tLagrangeOrder} });
    tParameters.set_lagrange_patterns({ {0} });
    tParameters.set_bspline_orders   ( { {tBsplineOrder} } );
    tParameters.set_bspline_patterns ( { {0} } );
    tParameters.set_side_sets( tSideSets );
    tParameters.set_refinement_buffer( 2 );
    tParameters.set_staircase_buffer( 2 );
    tParameters.set_lagrange_to_bspline_mesh( {{ {0} }});

    hmr::HMR tHMR(tParameters);
    Cell<std::shared_ptr< moris::hmr::Field > > tHMRFields;

    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

    tHMRFields.resize(2);

    // create field
    tHMRFields(0) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    tHMRFields(1) = tMesh->create_field( "Geom", tLagrangeMeshIndex );

    tHMRFields(0)->evaluate_scalar_function( Circle4MatMDL );
    tHMRFields(1)->evaluate_scalar_function( Plane4MatMDL1 );

    for( uint k=0; k<tRefinementLevels; ++k )
    {
        tHMR.flag_surface_elements_on_working_pattern( tHMRFields(0) );
        tHMR.flag_surface_elements_on_working_pattern( tHMRFields(1) );

        tHMR.perform_refinement_based_on_working_pattern( 0 );

        tHMRFields(0)->evaluate_scalar_function( Circle4MatMDL );
        tHMRFields(1)->evaluate_scalar_function( Plane4MatMDL1 );
    }

    tHMR.finalize();

    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

    xtk::Geom_Field tCircleFieldAsGeom(tHMRFields(0));
    xtk::Geom_Field tPlaneFieldAsGeom2(tHMRFields(1));
    moris::Cell<xtk::Geometry*> tGeometryVector = {&tCircleFieldAsGeom,&tPlaneFieldAsGeom2};

    size_t tModelDimension = tDim;
    xtk::Phase_Table     tPhaseTable (tGeometryVector.size());
    xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
    xtk::Model           tXTKModel(tModelDimension,tInterpMesh.get(),&tGeometryEngine);
    tXTKModel.mVerbose = false;

    tXTKModel.decompose(tDecompositionMethods);

    tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);

    xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
    xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

    // place the pair in mesh manager
    std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
    tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

    Cell< fem::IWG_User_Defined_Info > tBulkIWG(1);
    Cell< fem::IWG_User_Defined_Info > tDBCIWG(1);
    Cell< fem::IWG_User_Defined_Info > tNBCIWG(1);
    Cell< fem::IWG_User_Defined_Info > tIntIWG(1);

    // create IWG user defined info
    tBulkIWG( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK,
            { MSI::Dof_Type::TEMP },
            {{ MSI::Dof_Type::TEMP }},
            {fem::Property_Type::TEMP_LOAD },
            { fem::Constitutive_Type::DIFF_LIN_ISO } );
    tDBCIWG( 0 )  = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE,
            { MSI::Dof_Type::TEMP },
            {{ MSI::Dof_Type::TEMP }},
            { fem::Property_Type::TEMP_DIRICHLET },
            { fem::Constitutive_Type::DIFF_LIN_ISO } );
     tNBCIWG( 0 )  = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_NEUMANN,
            { MSI::Dof_Type::TEMP },
            {{ MSI::Dof_Type::TEMP }},
            { fem::Property_Type::TEMP_NEUMANN },
            moris::Cell< fem::Constitutive_Type >( 0 ) );
     tIntIWG( 0 )  = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE,
            { MSI::Dof_Type::TEMP },
            {{ MSI::Dof_Type::TEMP }},
            Cell< fem::Property_Type >( 0 ),
            {fem::Constitutive_Type::DIFF_LIN_ISO },
            {{ MSI::Dof_Type::TEMP }},
            Cell< fem::Property_Type >( 0 ),
            {fem::Constitutive_Type::DIFF_LIN_ISO } );

    Cell< Cell< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 14 );

    tIWGUserDefinedInfo( 0 )  = tBulkIWG;
    tIWGUserDefinedInfo( 1 )  = tBulkIWG;
    tIWGUserDefinedInfo( 2 )  = tBulkIWG;
    tIWGUserDefinedInfo( 3 )  = tBulkIWG;
    tIWGUserDefinedInfo( 4 )  = tBulkIWG;
    tIWGUserDefinedInfo( 5 )  = tBulkIWG;
    tIWGUserDefinedInfo( 6 )  = tBulkIWG;
    tIWGUserDefinedInfo( 7 )  = tBulkIWG;
    tIWGUserDefinedInfo( 8 )  = tDBCIWG;
    tIWGUserDefinedInfo( 9 )  = tNBCIWG;
    tIWGUserDefinedInfo( 10 ) = tIntIWG;
    tIWGUserDefinedInfo( 11 ) = tIntIWG;
    tIWGUserDefinedInfo( 12 ) = tIntIWG;
    tIWGUserDefinedInfo( 13 ) = tIntIWG;

    // create the property user defined infos
    fem::Property_User_Defined_Info tConductivity( fem::Property_Type::CONDUCTIVITY,
            Cell< Cell< MSI::Dof_Type > >( 0 ),
            {{{ tKappa_1 }}},
            tConstValFunction2MatMDL,
            Cell< fem::PropertyFunc >( 0 ) );

    fem::Property_User_Defined_Info tConductivity2( fem::Property_Type::CONDUCTIVITY,
            Cell< Cell< MSI::Dof_Type > >( 0 ),
            {{{ tKappa_2 }}},
            tConstValFunction2MatMDL,
            Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tTempDirichlet( fem::Property_Type::TEMP_DIRICHLET,
            Cell< Cell< MSI::Dof_Type > >( 0 ),
            {{{ tDirichlet }}},
            tConstValFunction2MatMDL,
            Cell< fem::PropertyFunc >( 0 ) );
    fem::Property_User_Defined_Info tNeumannFlux( fem::Property_Type::TEMP_NEUMANN,
            Cell< Cell< MSI::Dof_Type > >( 0 ),
            {{{ tNeumann }}},
            tConstValFunction2MatMDL,
            Cell< fem::PropertyFunc >( 0 ) );

    fem::Property_User_Defined_Info tTempLoad1( fem::Property_Type::TEMP_LOAD,
            Cell< Cell< MSI::Dof_Type > >( 0 ),
            {{{ tLoad_1 }}},
            tConstValFunction2MatMDL,
            Cell< fem::PropertyFunc >( 0 ) );

    fem::Property_User_Defined_Info tTempLoad2( fem::Property_Type::TEMP_LOAD,
            Cell< Cell< MSI::Dof_Type > >( 0 ),
            {{{ tLoad_2 }}},
            tConstValFunction2MatMDL,
            Cell< fem::PropertyFunc >( 0 ) );

    // create property user defined info
    Cell< Cell< Cell< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 14 );
    tPropertyUserDefinedInfo(0)  = create_bulk_properties(tConductivity2,tTempLoad2);
    tPropertyUserDefinedInfo(1)  = create_bulk_properties(tConductivity2,tTempLoad2);
    tPropertyUserDefinedInfo(2)  = create_bulk_properties(tConductivity2,tTempLoad2);
    tPropertyUserDefinedInfo(3)  = create_bulk_properties(tConductivity2,tTempLoad2);
    tPropertyUserDefinedInfo(4)  = create_bulk_properties(tConductivity,tTempLoad1);
    tPropertyUserDefinedInfo(5)  = create_bulk_properties(tConductivity,tTempLoad1);
    tPropertyUserDefinedInfo(6)  = create_bulk_properties(tConductivity,tTempLoad1);
    tPropertyUserDefinedInfo(7)  = create_bulk_properties(tConductivity,tTempLoad1);
    tPropertyUserDefinedInfo(8)  = create_dirichlet_properties(tConductivity,tTempDirichlet);
    tPropertyUserDefinedInfo(9)  = create_neumann_properties(tNeumannFlux);
    tPropertyUserDefinedInfo(10) = create_interface_properties(tConductivity2,tConductivity2);
    tPropertyUserDefinedInfo(11) = create_interface_properties(tConductivity2,tConductivity);
    tPropertyUserDefinedInfo(12) = create_interface_properties(tConductivity2,tConductivity);
    tPropertyUserDefinedInfo(13) = create_interface_properties(tConductivity,tConductivity);

    // create constitutive user defined info
    fem::Constitutive_User_Defined_Info tDiffLinIso = create_diff_lin_constitutive_info();
    // create constitutive user defined info
    Cell< Cell< Cell< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 14 );
    tConstitutiveUserDefinedInfo(0) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(1) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(2) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(3) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(4) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(5) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(6) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(7) = create_bulk_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo(8) = create_dbc_diff_lin_constitutive(tDiffLinIso);
    tConstitutiveUserDefinedInfo( 9 ).resize( 1 ); // neumann
    tConstitutiveUserDefinedInfo( 10 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
    tConstitutiveUserDefinedInfo( 11 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
    tConstitutiveUserDefinedInfo( 12 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);
    tConstitutiveUserDefinedInfo( 13 ) = create_interface_diff_lin_constitutive(tDiffLinIso,tDiffLinIso);

    // create a list of active block-sets
    std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
    std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
    std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
    std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);

    moris::Cell< moris_index >  tSetList = {  tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p0"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p0"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p1"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p1"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p2"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p2"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p3"),
            tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p3"),
            tEnrIntegMesh.get_set_index_by_name("SideSet_4_n_p2"),
            tEnrIntegMesh.get_set_index_by_name("SideSet_2_n_p3"),
            tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName01),
            tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName02),
            tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName13),
            tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName23)};

    moris::Cell< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::BULK,
            fem::Element_Type::SIDESET,
            fem::Element_Type::SIDESET,
            fem::Element_Type::DOUBLE_SIDESET,
            fem::Element_Type::DOUBLE_SIDESET,
            fem::Element_Type::DOUBLE_SIDESET,
            fem::Element_Type::DOUBLE_SIDESET,
    };

    // create model
    mdl::Model * tModel = new mdl::Model( tMeshManager, tBSplineMeshIndex, tSetList, tSetTypeList, tIWGUserDefinedInfo, tPropertyUserDefinedInfo, tConstitutiveUserDefinedInfo, 0, false);

    moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 1: create linear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    dla::Solver_Factory  tSolFactory;
    std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

    tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
    tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
    tLinearSolverAlgorithm->set_param("AZ_orthog") = 1;
    tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
    tLinearSolverAlgorithm->set_param("AZ_precond") = AZ_dom_decomp;
    tLinearSolverAlgorithm->set_param("AZ_ilut_fill") = 5.0;
    tLinearSolverAlgorithm->set_param("AZ_max_iter") = 100;
    tLinearSolverAlgorithm->set_param("rel_residual") = 1e-8;

    dla::Linear_Solver tLinSolver;

    tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // STEP 2: create nonlinear solver and algorithm
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    NLA::Nonlinear_Solver_Factory tNonlinFactory;
    std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

    tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

    NLA::Nonlinear_Solver tNonlinearSolver;
    tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 10;
    tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
    tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
    tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;
    tNonlinearSolverAlgorithm->set_param("NLA_rel_residual") = 1e-6;

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
    Matrix<DDRMat> tFullSol;
    tTimeSolver.get_full_solution(tFullSol);

    delete tModel;

    // finalize moris global communication manager
    gMorisComm.finalize();

    return 0;

}

