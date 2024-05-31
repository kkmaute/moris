/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_Transient.cpp
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

#include "cl_FEM_IWG_Factory.hpp"                   //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"                   //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"                 //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src
#include "cl_HMR_Database.hpp"

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
#include "cl_GEN_Geometry_Field_HMR.hpp"

// PRM
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "cl_GEN_Line.hpp"

#include "fn_norm.hpp"

namespace moris
{

    // define free function for properties
    inline void
    tPropConstFunc_MDLTransient( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&       aParameters,
            moris::fem::Field_Interpolator_Manager*              aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    inline void
    tPropTimeFunc_MDLTransient(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tTime  = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );
        aPropMatrix = aParameters( 0 ) * tTime;
    }

    inline bool
    tSolverOutputCriteria_MDLTransient( moris::tsa::Time_Solver* )
    {
        return true;
    }

    TEST_CASE( "MDL Transient", "[MDL_Transient]" )
    {
        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;

            // empty container for B-Spline meshes
            Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

            // create settings object
            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 1 }, { 1 } } );
            tParameters.set_domain_dimensions( 1, 1 );
            tParameters.set_domain_offset( 0.0, 0.0 );
            tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );

            tParameters.set_bspline_truncation( true );
            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );
            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_staircase_buffer( 1 );
            tParameters.set_initial_refinement( { { 2 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );
            tParameters.set_number_aura( true );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            // create the HMR object by passing the settings to the constructor
            moris::hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            tHMR.finalize();

            // construct a mesh manager for the fem
            moris::hmr::Interpolation_Mesh_HMR* tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            moris::hmr::Integration_Mesh_HMR*   tIGMesh = tHMR.create_integration_mesh( 1, 0, tIPMesh );

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( tIPMesh, tIGMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
            tPropConductivity->set_parameters( { { { 1.0 } } } );
            tPropConductivity->set_val_function( tPropConstFunc_MDLTransient );

            std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
            tPropDensity->set_parameters( { { { 1.0 } } } );
            tPropDensity->set_val_function( tPropConstFunc_MDLTransient );

            std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
            tPropHeatCapacity->set_parameters( { { { 1.0 } } } );
            tPropHeatCapacity->set_val_function( tPropConstFunc_MDLTransient );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { 100.0 } } } );
            tPropNeumann->set_val_function( tPropTimeFunc_MDLTransient );

            std::shared_ptr< fem::Property > tPropInitCondition = std::make_shared< fem::Property >();
            tPropInitCondition->set_parameters( { { { 0.0 } } } );
            tPropInitCondition->set_val_function( tPropConstFunc_MDLTransient );

            std::shared_ptr< fem::Property > tPropWeightCurrent = std::make_shared< fem::Property >();
            tPropWeightCurrent->set_parameters( { { { 100.0 } } } );
            tPropWeightCurrent->set_val_function( tPropConstFunc_MDLTransient );

            std::shared_ptr< fem::Property > tPropWeightPrevious = std::make_shared< fem::Property >();
            tPropWeightPrevious->set_parameters( { { { 100.0 } } } );
            tPropWeightPrevious->set_val_function( tPropConstFunc_MDLTransient );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
                    tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffusion->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
            tCMDiffusion->set_space_dim( 2 );
            tCMDiffusion->set_local_properties();

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGDiffusionBulk =
                    tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGDiffusionBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDiffusionBulk->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDiffusionBulk->set_constitutive_model( tCMDiffusion, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDiffusionBulk->set_property( tPropDensity, "Density", mtk::Leader_Follower::LEADER );
            tIWGDiffusionBulk->set_property( tPropHeatCapacity, "HeatCapacity", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGNeumann =
                    tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGTimeContinuity =
                    tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
            tIWGTimeContinuity->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGTimeContinuity->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGTimeContinuity->set_property( tPropWeightCurrent, "WeightCurrent", mtk::Leader_Follower::LEADER );
            tIWGTimeContinuity->set_property( tPropWeightPrevious, "WeightPrevious", mtk::Leader_Follower::LEADER );
            tIWGTimeContinuity->set_property( tPropInitCondition, "InitialCondition", mtk::Leader_Follower::LEADER );

            // define the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP =
                    tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQITEMP->set_output_type_index( 0 );

            // define set info
            Vector< fem::Set_User_Info > tSetInfo( 3 );

            tSetInfo( 0 ).set_mesh_index( 0 );
            tSetInfo( 0 ).set_IWGs( { tIWGDiffusionBulk } );
            tSetInfo( 0 ).set_IQIs( { tIQITEMP } );

            tSetInfo( 1 ).set_mesh_index( 2 );
            tSetInfo( 1 ).set_IWGs( { tIWGNeumann } );

            tSetInfo( 2 ).set_mesh_index( 0 );
            tSetInfo( 2 ).set_time_continuity( true );
            tSetInfo( 2 ).set_IWGs( { tIWGTimeContinuity } );

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_Transient.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy" },
                    { "Temperature" },
                    { vis::Field_Type::NODAL },
                    { vis::Output_Type::TEMP } );
            tModel->set_output_manager( &tOutputData );

            // --------------------------------------------------------------------------------------
            // define linear solver and algorithm
            dla::Solver_Factory tSolFactory;

            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm =
                    tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // --------------------------------------------------------------------------------------
            // define nonlinear solver and algorithm
            NLA::Nonlinear_Solver_Factory tNonlinFactory;

            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm =
                    tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolver;
            tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

            // --------------------------------------------------------------------------------------
            // define time solver and algorithm
            tsa::Time_Solver_Factory tTimeSolverFactory;

            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm =
                    tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
            tTimeSolverAlgorithm->set_param( "TSA_Num_Time_Steps" ) = 100;
            tTimeSolverAlgorithm->set_param( "TSA_Time_Frame" )     = 1.0;

            tsa::Time_Solver tTimeSolver;
            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            tTimeSolver.set_param( "TSA_Initialize_Sol_Vec" )  = "TEMP,0.0";
            tTimeSolver.set_param( "TSA_time_level_per_type" ) = "TEMP,2";

            sol::SOL_Warehouse tSolverWarehouse;
            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tTimeSolver.set_dof_type_list( { { MSI::Dof_Type::TEMP } } );

            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLTransient );

            // --------------------------------------------------------------------------------------
            // solve and check
            tTimeSolver.solve();

            // clean up
            delete tIPMesh;
            delete tIGMesh;
        }

    } /* END_TEST_CASE */

    TEST_CASE( "MDL Transient XFEM", "[MDL_Transient_XFEM]" )
    {
        // Geometry Parameters
        moris::real tDomainLX   = 2.0; /* Length of full domain in x (m) */
        moris::real tDomainLY   = 1.0; /* Length of full domain in y (m) */
        moris::real tPlaneLeft  = 0.0; /* x left plane   (m) */
        moris::real tPlaneRight = 1.0; /* x right plane  (m) */

        // Mesh Setup
        moris::uint tNumX   = 2; /* Number of elements in x*/
        moris::uint tNumY   = 1; /* Number of elements in y*/
        moris::uint tNumRef = 0; /* Number of HMR refinements */
        moris::uint tOrder  = 1; /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint          tLagrangeMeshIndex = 0;
        ParameterList tParameters        = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", std::to_string( tNumX ) + "," + std::to_string( tNumY ) );
        tParameters.set( "domain_dimensions", std::to_string( tDomainLX ) + "," + std::to_string( tDomainLY ) );
        tParameters.set( "domain_offset", std::to_string( -tDomainLX / 2.1 ) + "," + std::to_string( 0 ) );
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
        tParameters.set( "initial_refinement", "2" );
        tParameters.set( "initial_refinement_pattern", "0" );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        // initial refinement
        tHMR.perform_initial_refinement();

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //    for( uint k=0; k<tNumRef; ++k )
        //    {
        //        Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry( 2 );
        //        tGeometry( 0 ) = std::make_shared< moris::gen::Plane >( tPlaneLeft, 0.0, 1.0, 0.0 );
        //        tGeometry( 1 ) = std::make_shared< moris::gen::Plane >( tPlaneRight, 0.0, 1.0, 0.0 );
        //
        //        moris::gen::Phase_Table tPhaseTable (1);
        //        moris::gen::Geometry_Engine tGENGeometryEngine( tGeometry, tPhaseTable, 2 );
        //
        //        moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );
        //
        //        uint tNumIPNodes = tMesh->get_num_nodes();
        //        Matrix<DDRMat> tFieldData( tNumIPNodes,1 );
        //        Matrix<DDRMat> tFieldData0( tNumIPNodes,1 );
        //
        //        tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumIPNodes );
        //        Matrix< DDRMat > tCoords( tNumIPNodes, 2 );
        //        for( uint i = 0; i < tNumIPNodes; i++ )
        //        {
        //            tCoords.set_row( i, tMesh->get_mtk_vertex(i).get_coords() );
        //        }
        //
        //        tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );
        //
        //        for(uint i=0; i<tNumIPNodes; i++)
        //        {
        //            tFieldData( i )  = tGENGeometryEngine.get_entity_phase_val( i, 0 );
        //            tFieldData0( i ) = tGENGeometryEngine.get_entity_phase_val( i, 1 );
        //        }
        //
        //        tHMR.based_on_field_put_elements_on_queue( tFieldData, tLagrangeMeshIndex );
        //        tHMR.based_on_field_put_elements_on_queue( tFieldData0, tLagrangeMeshIndex );
        //
        //        tHMR.perform_refinement_based_on_working_pattern( 0, false );
        //    }
        tHMR.finalize();

        moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        //-----------------------------------------------------------------------------------------------

        Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry0( 2 );
        tGeometry0( 0 ) = std::make_shared< moris::gen::Line >( tPlaneLeft, 0.0, 1.0, 0.0 );
        tGeometry0( 1 ) = std::make_shared< moris::gen::Line >( tPlaneRight, 0.0, 1.0, 0.0 );

        size_t                     tModelDimension = 2;
        moris::gen::Geometry_Engine tGENGeometryEngine0( tGeometry0, tModelDimension );

        // --------------------------------------------------------------------------------------
        xtk::Model tXTKModel( tModelDimension, tInterpolationMesh, &tGENGeometryEngine0 );
        tXTKModel.mVerbose = true;

        // Specify decomposition Method and Cut Mesh ---------------------------------------
        Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
        tXTKModel.decompose( tDecompositionMethods );

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
        //        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // get meshes for FEM
        xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
        tPropConductivity->set_parameters( { { { 1.0 } } } );
        tPropConductivity->set_val_function( tPropConstFunc_MDLTransient );

        std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
        tPropDensity->set_parameters( { { { 1.0 } } } );
        tPropDensity->set_val_function( tPropConstFunc_MDLTransient );

        std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
        tPropHeatCapacity->set_parameters( { { { 1.0 } } } );
        tPropHeatCapacity->set_val_function( tPropConstFunc_MDLTransient );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { { { 100.0 } } } );
        tPropNeumann->set_val_function( tPropTimeFunc_MDLTransient );

        std::shared_ptr< fem::Property > tPropInitCondition = std::make_shared< fem::Property >();
        tPropInitCondition->set_parameters( { { { 0.0 } } } );
        tPropInitCondition->set_val_function( tPropConstFunc_MDLTransient );

        std::shared_ptr< fem::Property > tPropWeightCurrent = std::make_shared< fem::Property >();
        tPropWeightCurrent->set_parameters( { { { 100.0 } } } );
        tPropWeightCurrent->set_val_function( tPropConstFunc_MDLTransient );

        std::shared_ptr< fem::Property > tPropWeightPrevious = std::make_shared< fem::Property >();
        tPropWeightPrevious->set_parameters( { { { 100.0 } } } );
        tPropWeightPrevious->set_val_function( tPropConstFunc_MDLTransient );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffusion->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
        tCMDiffusion->set_space_dim( 2 );
        tCMDiffusion->set_local_properties();

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGDiffusionBulk =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGDiffusionBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDiffusionBulk->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGDiffusionBulk->set_constitutive_model( tCMDiffusion, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDiffusionBulk->set_property( tPropDensity, "Density", mtk::Leader_Follower::LEADER );
        tIWGDiffusionBulk->set_property( tPropHeatCapacity, "HeatCapacity", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGTimeContinuity =
                tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
        tIWGTimeContinuity->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGTimeContinuity->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGTimeContinuity->set_property( tPropWeightCurrent, "WeightCurrent", mtk::Leader_Follower::LEADER );
        tIWGTimeContinuity->set_property( tPropWeightPrevious, "WeightPrevious", mtk::Leader_Follower::LEADER );
        tIWGTimeContinuity->set_property( tPropInitCondition, "InitialCondition", mtk::Leader_Follower::LEADER );

        // define the IQIs
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIQITEMP->set_output_type_index( 0 );

        // define set info
        Vector< fem::Set_User_Info > tSetInfo( 5 );

        tSetInfo( 0 ).set_mesh_set_name( "HMR_dummy_c_p2" );
        tSetInfo( 0 ).set_IWGs( { tIWGDiffusionBulk } );
        tSetInfo( 0 ).set_IQIs( { tIQITEMP } );

        tSetInfo( 1 ).set_mesh_set_name( "HMR_dummy_n_p2" );
        tSetInfo( 1 ).set_IWGs( { tIWGDiffusionBulk } );
        tSetInfo( 1 ).set_IQIs( { tIQITEMP } );

        tSetInfo( 2 ).set_mesh_set_name( "iside_1_b0_2_b1_3" );
        tSetInfo( 2 ).set_IWGs( { tIWGNeumann } );

        tSetInfo( 3 ).set_mesh_set_name( "HMR_dummy_c_p2" );
        tSetInfo( 3 ).set_time_continuity( true );
        tSetInfo( 3 ).set_IWGs( { tIWGTimeContinuity } );

        tSetInfo( 4 ).set_mesh_set_name( "HMR_dummy_n_p2" );
        tSetInfo( 4 ).set_time_continuity( true );
        tSetInfo( 4 ).set_IWGs( { tIWGTimeContinuity } );

        // create model
        mdl::Model* tModel = new mdl::Model( tMeshManager,
                0,
                tSetInfo );

        // --------------------------------------------------------------------------------------
        // define outputs
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                vis::VIS_Mesh_Type::STANDARD,
                "./",
                "UT_MDL_Transient_XFEM.exo",
                "./",
                "temp.exo",
                { "HMR_dummy_n_p2", "HMR_dummy_c_p2" },
                { "Temperature" },
                { vis::Field_Type::NODAL },
                { vis::Output_Type::TEMP } );
        tModel->set_output_manager( &tOutputData );

        // --------------------------------------------------------------------------------------
        // define linear solver and algorithm
        dla::Solver_Factory tSolFactory;

        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm =
                tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

        dla::Linear_Solver tLinSolver;
        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // --------------------------------------------------------------------------------------
        // define nonlinear solver and algorithm
        NLA::Nonlinear_Solver_Factory tNonlinFactory;

        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm =
                tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // --------------------------------------------------------------------------------------
        // define time solver and algorithm
        tsa::Time_Solver_Factory tTimeSolverFactory;

        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm =
                tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
        tTimeSolverAlgorithm->set_param( "TSA_Num_Time_Steps" ) = 100;
        tTimeSolverAlgorithm->set_param( "TSA_Time_Frame" )     = 1.0;

        tsa::Time_Solver tTimeSolver;
        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        tTimeSolver.set_param( "TSA_Initialize_Sol_Vec" )  = "TEMP,0.0";
        tTimeSolver.set_param( "TSA_time_level_per_type" ) = "TEMP,2";

        sol::SOL_Warehouse tSolverWarehouse;
        tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tTimeSolver.set_dof_type_list( { { MSI::Dof_Type::TEMP } } );

        tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLTransient );

        // --------------------------------------------------------------------------------------
        // solve and check
        tTimeSolver.solve();

        // clean up
        //------------------------------------------------------------------------------
        delete tInterpolationMesh;
        delete tModel;

    } /*END_TEST_CASE */
}    // namespace moris
