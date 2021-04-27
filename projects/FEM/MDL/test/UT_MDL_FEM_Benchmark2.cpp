/*
 * UT_MDL_FEM_Benchmark2.cpp
 *
 *  Created on: Feb 04, 2020
 *      Author: noel
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
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
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "fn_norm.hpp"

#include "cl_GEN_Circle.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

#include <functional>

namespace moris
{

    //-------------------------------------------------------------------------------------

    void ConstFunctionVal_MDLFEMBench2
    ( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void AnalyticalTempFunc_MDLFEMBench2
    ( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // get parameters
        real RInner  = aParameters( 0 )( 0 ); // inner radius
        real ROuter  = aParameters( 1 )( 0 ); // outer radius
        real xCenter = aParameters( 2 )( 0 ); // x coord of center
        real yCenter = aParameters( 2 )( 1 ); // y coord of center
        real TInner  = aParameters( 3 )( 0 ); // imposed temperature at inner radius
        real Q       = aParameters( 4 )( 0 ) * 2 * M_PI * ROuter; // heat load (W)
        real kappa   = aParameters( 5 )( 0 ); // conductivity (W/m^2)

        // get x and y coords
        real xCoord = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real yCoord = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // compute radius
        real R = std::sqrt( std::pow( xCoord - xCenter, 2 ) + std::pow( yCoord - yCenter, 2 ) );

        aPropMatrix = { { TInner + ( Q * std::log( R/RInner ) )/( kappa * 2 * M_PI ) } };
    }

    void AnalyticalTemp2MatFunc_MDLFEMBench2
    ( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        // get parameters
        real RInner  = aParameters( 0 )( 0 ); // inner radius
        real RMiddle = aParameters( 1 )( 0 ); // middle radius
        real ROuter  = aParameters( 2 )( 0 ); // outer radius
        real xCenter = aParameters( 3 )( 0 ); // x coord of center
        real yCenter = aParameters( 3 )( 1 ); // y coord of center
        real TInner  = aParameters( 4 )( 0 ); // imposed temperature at inner radius
        real Q       = aParameters( 5 )( 0 ) * 2 * M_PI * ROuter; // heat load (W)
        real kappaA  = aParameters( 6 )( 0 ); // conductivity for phase A (W/m^2)
        real kappaB  = aParameters( 7 )( 0 ); // conductivity for pahse B (W/m^2)

        // get x and y coords
        real xCoord = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real yCoord = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        // compute radius
        real R = std::sqrt( std::pow( xCoord - xCenter, 2 ) + std::pow( yCoord - yCenter, 2 ) );

        Matrix< DDRMat > tT;
        if( R < RMiddle )
        {
            tT = { { TInner + ( Q * ( std::log( R/RInner ) * ( 1/kappaB ) ) ) / ( 2 * M_PI ) } };
        }
        else
        {
            tT = { { TInner + ( Q * ( std::log( RMiddle/RInner ) * ( 1/kappaB ) + std::log( R/RMiddle ) * ( 1/kappaA ) ) ) / ( 2 * M_PI ) } };
        }

        aPropMatrix = tT;
    }

    bool tSolverOutputCriteria_MDLFEMBench2( moris::tsa::Time_Solver * )
    {
        return true;
    }
    //-------------------------------------------------------------------------------------

    TEST_CASE("MDL_FEM_Benchmark_Diffusion_1Mat","[MDL_FEM_Benchmark_Diffusion_1Mat]")
    {
        if(par_size()<=1)
        {
            // Geometry Parameters
            moris::real tDomainLX = 10.0;                   /* Length of full domain in x (m) */
            moris::real tDomainLY = 10.0;                   /* Length of full domain in y (m) */
            Matrix<DDRMat> tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
            moris::real tRInner = 0.55;                     /* Inner circle radius (m) */
            moris::real tROuter = 1.05;                     /* Outer circle radius (m) */

            //Material Parameters
            moris::real tKappaA = 1.0; /* Conductivity material A (W/m^2) */

            // Boundary Conditions
            moris::real tTDirichlet = 5.0; /* Imposed temperature for Dirichlet BC (K) */
            moris::real tDBCGamma = 100.0; /* Penalty for Dirichlet BC */
            moris::real tQ = 20.0;         /* Imposed heat flux for Neumann BC (W/m) */

            // Mesh Setup
            moris::uint tNumX   = 20; /* Number of elements in x*/
            moris::uint tNumY   = 20; /* Number of elements in y*/


            uint tLagrangeMeshIndex = 0;
            std::string tOuterFieldName   = "Outercircle";
            std::string tInnerFieldName   = "Innercircle";
            ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
            tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
            tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) );
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
            tParameters.set( "initial_refinement", "1" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );
            tParameters.set( "use_number_aura", 0 );

            std::shared_ptr<hmr::HMR> tHMR = std::make_shared<hmr::HMR>( tParameters );

            // Initial refinement
            tHMR->perform_initial_refinement();

            // Create geometry engine
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(2);
            tGeometry(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
            tGeometry(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);

            // Perform additional refinement
            //tGENGeometryEngine.perform_refinement(tHMR);

            // Get interpolation mesh
            tHMR->finalize();
            moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR->create_interpolation_mesh(tLagrangeMeshIndex);

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGENGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);

            //-----------------------------------------------------------------------------------------------

            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry0(2);
            tGeometry0(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
            tGeometry0(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);

            size_t tModelDimension = 2;
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters0;
            tGeometryEngineParameters0.mGeometries = tGeometry0;
            moris::ge::Geometry_Engine tGENGeometryEngine0(tInterpolationMesh, tGeometryEngineParameters0);

            // --------------------------------------------------------------------------------------
            xtk::Model tXTKModel(tModelDimension,tInterpolationMesh,&tGENGeometryEngine0);
            tXTKModel.mVerbose = true;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
            //        tXTKModel.construct_face_oriented_ghost_penalization_cells();

            xtk::Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = false;
            tOutputOptions.mAddSideSets = true;
            tOutputOptions.mAddClusters = false;

            // output integration mesh
            moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
            std::string tOutputFile = "./mdl_exo/MDL_FEM_Benchmark_Diffusion_1Mat.exo";
            tIntegMesh1->create_output_mesh(tOutputFile);

            // get meshes for FEM
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

            // create the properties
            // --------------------------------------------------------------------------------------
            std::shared_ptr< fem::Property > tPropKappaA = std::make_shared< fem::Property >();
            tPropKappaA->set_parameters( { {{ tKappaA }} } );
            tPropKappaA->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { {{ tTDirichlet }} } );
            tPropDirichlet->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { tQ } } } );
            tPropNeumann->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropL2Analytic = std::make_shared< fem::Property >();
            tPropL2Analytic->set_parameters( { {{ tRInner }}, {{ tROuter }}, tCenterPoint, {{ tTDirichlet }}, {{ tQ }}, {{ tKappaA }} } );
            tPropL2Analytic->set_val_function( AnalyticalTempFunc_MDLFEMBench2 );

            // create constitutive models
            // --------------------------------------------------------------------------------------
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIsoA = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIsoA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIsoA->set_property( tPropKappaA, "Conductivity" );
            tCMDiffLinIsoA->set_space_dim( 2 );
            tCMDiffLinIsoA->set_local_properties();

            // create stabilization parameters
            // --------------------------------------------------------------------------------------
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { {{ tDBCGamma }} } );
            tSPDirichletNitsche->set_property( tPropKappaA, "Material", mtk::Master_Slave::MASTER );

            // create the IWGs
            // --------------------------------------------------------------------------------------
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulkA = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkA->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulkA->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            // create the IQIs
            // --------------------------------------------------------------------------------------
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            std::shared_ptr< fem::IQI > tIQIL2 = tIQIFactory.create_IQI( fem::IQI_Type::L2_ERROR_ANALYTIC );
            tIQIL2->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQIL2->set_property( tPropL2Analytic, "L2Check", mtk::Master_Slave::MASTER );
            tIQIL2->set_name( "IQI_L2" );

            std::shared_ptr< fem::IQI > tIQITempExact = tIQIFactory.create_IQI( fem::IQI_Type::PROPERTY );
            tIQITempExact->set_property( tPropL2Analytic, "Property", mtk::Master_Slave::MASTER );
            tIQITempExact->set_name( "IQI_Exact" );

            // create set info
            // --------------------------------------------------------------------------------------
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk1.set_IWGs( { tIWGBulkA } );
            tSetBulk1.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk2.set_IWGs( { tIWGBulkA } );
            tSetBulk2.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetDirichlet1;
            tSetDirichlet1.set_mesh_set_name( "iside_b0_1_b1_0" );
            tSetDirichlet1.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann1;
            tSetNeumann1.set_mesh_set_name( "iside_b0_1_b1_3" );
            tSetNeumann1.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
            tSetInfo( 0 )  = tSetBulk1;
            tSetInfo( 1 )  = tSetBulk2;
            tSetInfo( 2 )  = tSetDirichlet1;
            tSetInfo( 3 )  = tSetNeumann1;

            // create model
            // --------------------------------------------------------------------------------------
            mdl::Model * tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0, false );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                    "./",
                    "MDL_FEM_Benchmark_Diffusion_1Mat.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p1", "HMR_dummy_n_p1" },
                    { "TEMP", "L2", "TEMP_EXACT" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_TEMP","IQI_L2","IQI_Exact" } );
            tModel->set_output_manager( &tOutputData );

            // create linear solver and algorithm
            // --------------------------------------------------------------------------------------
            // define dof type for solve
            moris::Cell< enum MSI::Dof_Type > tDofTypesU( 1 );
            tDofTypesU( 0 ) = MSI::Dof_Type::TEMP;

            dla::Solver_Factory  tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
            //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            //        = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            //        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10000;
            //        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
            //        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
            //        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // create nonlinear solver and algorithm
            // --------------------------------------------------------------------------------------
            NLA::Nonlinear_Solver_Factory tNonlinFactory;

            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
            = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU
            //        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolverMain;
            tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
            tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

            // create solver database
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
            tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

            // create time Solver and algorithm
            // --------------------------------------------------------------------------------------
            tsa::Time_Solver_Factory tTimeSolverFactory;

            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
            = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

            tsa::Time_Solver tTimeSolver;
            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_dof_type_list( tDofTypesU );
            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench2 );

            // solve
            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            //         //print full solution for debug
            //         Matrix<DDRMat> tFullSol;
            //         tTimeSolver.get_full_solution(tFullSol);
            //         print(tFullSol,"tFullSol");

            // clean up
            //------------------------------------------------------------------------------
            delete tInterpolationMesh;
            delete tIntegMesh1;
            delete tModel;
        }
    }

    TEST_CASE("MDL_FEM_Benchmark_Diffusion_1Mat_Ghost","[MDL_FEM_Benchmark_Diffusion_1Mat_Ghost]")
    {
        if(par_size()<=1)
        {
            // Geometry Parameters
            moris::real tDomainLX = 10.0;                   /* Length of full domain in x (m) */
            moris::real tDomainLY = 10.0;                   /* Length of full domain in y (m) */
            Matrix<DDRMat> tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
            moris::real tRInner = 0.55;                     /* Inner circle radius (m) */
            moris::real tROuter = 1.05;                     /* Outer circle radius (m) */

            //Material Parameters
            moris::real tKappaA = 1.0; /* Conductivity material A (W/m^2) */

            // Boundary Conditions
            moris::real tTDirichlet = 5.0; /* Imposed temperature for Dirichlet BC (K) */
            moris::real tDBCGamma = 100.0; /* Penalty for Dirichlet BC */
            moris::real tQ = 20.0;         /* Imposed heat flux for Neumann BC (W/m) */
            moris::real tGammaDisplGhost = 0.0001; /* Penalty for displacement based Ghost */

            // Mesh Setup
            moris::uint tNumX   = 20; /* Number of elements in x*/
            moris::uint tNumY   = 20; /* Number of elements in y*/
            moris::uint tLagrangeMeshIndex = 0;

            std::string tOuterFieldName   = "Outercircle";
            std::string tInnerFieldName   = "Innercircle";
            ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
            tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
            tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) );
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
            tParameters.set( "initial_refinement", "1" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );
            tParameters.set( "use_number_aura", 0 );

            std::shared_ptr<hmr::HMR> tHMR = std::make_shared<hmr::HMR>( tParameters );

            // Initial refinement
            tHMR->perform_initial_refinement();

            // Create geometry engine
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(2);
            tGeometry(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
            tGeometry(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);

            // Perform additional refinement
            //tGENGeometryEngine.perform_refinement(tHMR);

            // Get interpolation mesh
            tHMR->finalize();
            moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR->create_interpolation_mesh(tLagrangeMeshIndex);

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGENGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);

            //-----------------------------------------------------------------------------------------------

            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry0(2);
            tGeometry0(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
            tGeometry0(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);

            size_t tModelDimension = 2;
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters0;
            tGeometryEngineParameters0.mGeometries = tGeometry0;
            moris::ge::Geometry_Engine tGENGeometryEngine0(tInterpolationMesh, tGeometryEngineParameters0);

            // --------------------------------------------------------------------------------------
            xtk::Model tXTKModel(tModelDimension,tInterpolationMesh,&tGENGeometryEngine0);
            tXTKModel.mVerbose = true;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
            tXTKModel.construct_face_oriented_ghost_penalization_cells();

            // get meshes for FEM
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();
            xtk::Ghost_Stabilization         & tGhost         = tXTKModel.get_ghost_stabilization( 0 );

            tGhost.visualize_ghost_on_mesh( 0 );
            tGhost.visualize_ghost_on_mesh( 1 );
            tGhost.visualize_ghost_on_mesh( 2 );

            moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
            writer.write_mesh("", "benchmark_enriched_ig_w_ghost.exo", "", "temp.exo");
            writer.close_file();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

            // create the properties
            // --------------------------------------------------------------------------------------
            std::shared_ptr< fem::Property > tPropKappaA = std::make_shared< fem::Property >();
            tPropKappaA->set_parameters( { {{ tKappaA }} } );
            tPropKappaA->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { {{ tTDirichlet }} } );
            tPropDirichlet->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { tQ } } } );
            tPropNeumann->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropL2Analytic = std::make_shared< fem::Property >();
            tPropL2Analytic->set_parameters( { {{ tRInner }}, {{ tROuter }}, tCenterPoint, {{ tTDirichlet }}, {{ tQ }}, {{ tKappaA }} } );
            tPropL2Analytic->set_val_function( AnalyticalTempFunc_MDLFEMBench2 );

            // create constitutive models
            // --------------------------------------------------------------------------------------
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIsoA = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIsoA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIsoA->set_property( tPropKappaA, "Conductivity" );
            tCMDiffLinIsoA->set_space_dim( 2 );
            tCMDiffLinIsoA->set_local_properties();

            // create stabilization parameters
            // --------------------------------------------------------------------------------------
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { {{ tDBCGamma }} } );
            tSPDirichletNitsche->set_property( tPropKappaA, "Material", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPDisplGhost =
                    tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
            tSPDisplGhost->set_parameters( {{{ tGammaDisplGhost }} });
            tSPDisplGhost->set_property( tPropKappaA, "Material", mtk::Master_Slave::MASTER );

            // create the IWGs
            // --------------------------------------------------------------------------------------
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulkA = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkA->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulkA->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGGhost = tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
            tIWGGhost->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGGhost->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGGhost->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );
            tIWGGhost->set_stabilization_parameter( tSPDisplGhost, "GhostSP" );

            // create the IQIs
            // --------------------------------------------------------------------------------------
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            std::shared_ptr< fem::IQI > tIQIL2 = tIQIFactory.create_IQI( fem::IQI_Type::L2_ERROR_ANALYTIC );
            tIQIL2->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQIL2->set_property( tPropL2Analytic, "L2Check", mtk::Master_Slave::MASTER );
            tIQIL2->set_name( "IQI_L2" );

            std::shared_ptr< fem::IQI > tIQITempExact = tIQIFactory.create_IQI( fem::IQI_Type::PROPERTY );
            tIQITempExact->set_property( tPropL2Analytic, "Property", mtk::Master_Slave::MASTER );
            tIQITempExact->set_name( "IQI_Exact" );

            // create set info
            // --------------------------------------------------------------------------------------
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulk1.set_IWGs( { tIWGBulkA } );
            tSetBulk1.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk2.set_IWGs( { tIWGBulkA } );
            tSetBulk2.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetDirichlet1;
            tSetDirichlet1.set_mesh_set_name( "iside_b0_1_b1_0" );
            tSetDirichlet1.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann1;
            tSetNeumann1.set_mesh_set_name( "iside_b0_1_b1_3" );
            tSetNeumann1.set_IWGs( { tIWGNeumann } );

            fem::Set_User_Info tSetDisplGhost;
            tSetDisplGhost.set_mesh_set_name( tGhost.get_ghost_dbl_side_set_name(1) );
            tSetDisplGhost.set_IWGs( { tIWGGhost } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 5 );
            tSetInfo( 0 )  = tSetBulk1;
            tSetInfo( 1 )  = tSetBulk2;
            tSetInfo( 2 )  = tSetDirichlet1;
            tSetInfo( 3 )  = tSetNeumann1;
            tSetInfo( 4 )  = tSetDisplGhost;

            // create model
            // --------------------------------------------------------------------------------------
            mdl::Model * tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0, false );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                    "./",
                    "MDL_FEM_Benchmark_Diffusion_1Mat_Ghost.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p1", "HMR_dummy_n_p1" },
                    { "TEMP", "L2", "TEMP_EXACT" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_TEMP","IQI_L2","IQI_Exact" } );
            tModel->set_output_manager( &tOutputData );

            // create linear solver and algorithm
            // --------------------------------------------------------------------------------------
            // define dof type for solve
            moris::Cell< enum MSI::Dof_Type > tDofTypesU( 1 );
            tDofTypesU( 0 ) = MSI::Dof_Type::TEMP;

            dla::Solver_Factory  tSolFactory;
            //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            //        = tSolFactory.create_solver( SolverType::AMESOS_IMPL );
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            //        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10000;
            //        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
            //        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
            //        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "";
            tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
            tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;
            tLinearSolverAlgorithm->set_param("AZ_precond") = AZ_none;
            tLinearSolverAlgorithm->set_param("AZ_kspace") = 500;

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // create nonlinear solver and algorithm
            // --------------------------------------------------------------------------------------
            NLA::Nonlinear_Solver_Factory tNonlinFactory;

            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
            = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU
            //        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolverMain;
            tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
            tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

            // create solver database
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
            tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

            // create time Solver and algorithm
            // --------------------------------------------------------------------------------------
            tsa::Time_Solver_Factory tTimeSolverFactory;

            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
            = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

            tsa::Time_Solver tTimeSolver;
            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_dof_type_list( tDofTypesU );
            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench2 );

            // solve
            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            //         //print full solution for debug
            //         Matrix<DDRMat> tFullSol;
            //         tTimeSolver.get_full_solution(tFullSol);
            //         print(tFullSol,"tFullSol");

            // clean up
            //------------------------------------------------------------------------------
            delete tModel;
            delete tInterpolationMesh;
        }
    }


    TEST_CASE("FEM Benchmark 2 - 2Mat","[MDL_FEM_Benchmark2_2Mat]")
    {
        if(par_size()<=1)
        {
            // define problem parameters
            //------------------------------------------------------------------------------
            // geometry Parameters
            moris::real tDomainLX  = 10.0;                   /* Length of full domain in x (m) */
            moris::real tDomainLY  = 10.0;                   /* Length of full domain in y (m) */
            moris::real tRInner    = 0.55;                   /* Inner circle radius (m) */
            moris::real tRMiddle   = 0.77;                   /* Middle circle radius (m) */
            moris::real tROuter    = 1.05;                   /* Outer circle radius (m) */
            Matrix<DDRMat> tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/

            // material Parameters
            moris::real tKappaA = 1.0; /* Conductivity material A (W/m^2) */
            moris::real tKappaB = 100.0; /* Conductivity material B (W/m^2) */

            // boundary Conditions
            moris::real tTDirichlet = 5.0;   /* Imposed temperature for Dirichlet BC (K) */
            moris::real tDBCGamma   = 100.0; /* Penalty for Dirichlet BC */
            moris::real tQ          = 20.0;  /* Imposed heat flux for Neumann BC (W/m) */
            moris::real tIBCGamma   = 100.0;   /* Penalty for Interface BC */

            // mesh Setup
            moris::uint tNumX   = 20; /* Number of elements in x*/
            moris::uint tNumY   = 20; /* Number of elements in y*/

            // define hmr parameters
            //------------------------------------------------------------------------------
            uint tLagrangeMeshIndex = 0;
            std::string tOuterFieldName  = "OuterCircle";
            std::string tMiddleFieldName = "MiddleCircle";
            std::string tInnerFieldName  = "InnerCircle";
            ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
            tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
            tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) );
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
            tParameters.set( "initial_refinement", "1" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );
            tParameters.set( "use_number_aura", 0 );

            std::shared_ptr<hmr::HMR> tHMR = std::make_shared<hmr::HMR>( tParameters );

            // Initial refinement
            tHMR->perform_initial_refinement();

            // Create geometry engine
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(3);
            tGeometry(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
            tGeometry(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRMiddle);
            tGeometry(2) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);

            // Perform additional refinement
            //tGENGeometryEngine.perform_refinement(tHMR);

            // Get interpolation mesh
            tHMR->finalize();
            moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR->create_interpolation_mesh(tLagrangeMeshIndex);

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGENGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);

            // --------------------------------------------------------------------------------------
            xtk::Model tXTKModel( 2, tInterpolationMesh, &tGENGeometryEngine );
            tXTKModel.mVerbose = true;

            // specify decomposition method and cut mesh
            Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );
            //        tXTKModel.construct_face_oriented_ghost_penalization_cells();

            xtk::Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = false;
            tOutputOptions.mAddSideSets = true;
            tOutputOptions.mAddClusters = false;

            // output integration mesh
            moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh( tOutputOptions );
            std::string tOutputFile = "./mdl_exo/FEM_Bench2_2Mat.exo";
            tIntegMesh1->create_output_mesh( tOutputFile );

            // get meshes for FEM
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

            // create the properties
            //------------------------------------------------------------------------------
            std::shared_ptr< fem::Property > tPropKappaA = std::make_shared< fem::Property >();
            tPropKappaA->set_parameters( { {{ tKappaA }} } );
            tPropKappaA->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropKappaB = std::make_shared< fem::Property >();
            tPropKappaB->set_parameters( { {{ tKappaB }} } );
            tPropKappaB->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { {{ tTDirichlet }} } );
            tPropDirichlet->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { tQ } } } );
            tPropNeumann->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropL2Analytic = std::make_shared< fem::Property >();
            tPropL2Analytic->set_parameters( { {{ tRInner }}, {{ tRMiddle }}, {{ tROuter }}, tCenterPoint, {{ tTDirichlet }}, {{ tQ }}, {{ tKappaA }}, {{ tKappaB }} } );
            tPropL2Analytic->set_val_function( AnalyticalTemp2MatFunc_MDLFEMBench2 );

            // create constitutive models
            //------------------------------------------------------------------------------
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIsoA = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIsoA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIsoA->set_property( tPropKappaA, "Conductivity" );
            tCMDiffLinIsoA->set_space_dim( 2 );
            tCMDiffLinIsoA->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIsoB = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIsoB->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIsoB->set_property( tPropKappaB, "Conductivity" );
            tCMDiffLinIsoB->set_space_dim( 2 );
            tCMDiffLinIsoB->set_local_properties();

            // create stabilization parameters
            //------------------------------------------------------------------------------
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { {{ tDBCGamma }} } );
            tSPDirichletNitsche->set_property( tPropKappaA, "Material", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { {{ tIBCGamma }} } );
            tSPNitscheInterface->set_property( tPropKappaB, "Material", mtk::Master_Slave::MASTER );
            tSPNitscheInterface->set_property( tPropKappaA, "Material", mtk::Master_Slave::SLAVE );

            // create the IWGs
            //------------------------------------------------------------------------------
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulkA = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkA->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulkA->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGBulkB = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkB->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkB->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulkB->set_constitutive_model( tCMDiffLinIsoB, "Diffusion", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIsoB, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
            tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );
            tIWGInterface->set_constitutive_model( tCMDiffLinIsoB, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGInterface->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::SLAVE );
            tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

            // create the IQIs
            //------------------------------------------------------------------------------
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_Temp" );

            std::shared_ptr< fem::IQI > tIQIL2 = tIQIFactory.create_IQI( fem::IQI_Type::L2_ERROR_ANALYTIC );
            tIQIL2->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQIL2->set_property( tPropL2Analytic, "L2Check", mtk::Master_Slave::MASTER );
            tIQIL2->set_name( "IQI_L2" );

            std::shared_ptr< fem::IQI > tIQITempExact = tIQIFactory.create_IQI( fem::IQI_Type::PROPERTY );
            tIQITempExact->set_property( tPropL2Analytic, "Property", mtk::Master_Slave::MASTER );
            tIQITempExact->set_name( "IQI_Exact" );

            // create set info
            //------------------------------------------------------------------------------
            fem::Set_User_Info tSetBulkB1;
            tSetBulkB1.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulkB1.set_IWGs( { tIWGBulkB } );
            tSetBulkB1.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulkB2;
            tSetBulkB2.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulkB2.set_IWGs( { tIWGBulkB } );
            tSetBulkB2.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulkA1;
            tSetBulkA1.set_mesh_set_name( "HMR_dummy_c_p3" );
            tSetBulkA1.set_IWGs( { tIWGBulkA } );
            tSetBulkA1.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulkA2;
            tSetBulkA2.set_mesh_set_name( "HMR_dummy_n_p3" );
            tSetBulkA2.set_IWGs( { tIWGBulkA } );
            tSetBulkA2.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetDirichletB1;
            tSetDirichletB1.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name( 2, 1, 0 ) );
            tSetDirichletB1.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumannA1;
            tSetNeumannA1.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name( 0, 3, 7 ) );
            tSetNeumannA1.set_IWGs( { tIWGNeumann } );

            fem::Set_User_Info tSetInterfaceBA1;
            tSetInterfaceBA1.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 1, 3 ) );
            tSetInterfaceBA1.set_IWGs( { tIWGInterface } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 7 );
            tSetInfo( 0 )  = tSetBulkB1;
            tSetInfo( 1 )  = tSetBulkB2;
            tSetInfo( 2 )  = tSetBulkA1;
            tSetInfo( 3 )  = tSetBulkA2;
            tSetInfo( 4 )  = tSetDirichletB1;
            tSetInfo( 5 )  = tSetNeumannA1;
            tSetInfo( 6 )  = tSetInterfaceBA1;

            // create model
            //------------------------------------------------------------------------------
            mdl::Model * tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0, false );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                    "./",
                    "UT_MDL_FEM_Bench2_Output_2Mat.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p1", "HMR_dummy_n_p1", "HMR_dummy_c_p3", "HMR_dummy_n_p3" },
                    { "TEMP", "L2", "TEMP_EXACT" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_Temp","IQI_L2","IQI_Exact" } );
            tModel->set_output_manager( &tOutputData );

            // create linear solver and algorithm
            // --------------------------------------------------------------------------------------
            // define dof type for solve
            moris::Cell< enum MSI::Dof_Type > tSolveDofTypes( 1 );
            tSolveDofTypes( 0 ) = MSI::Dof_Type::TEMP;

            dla::Solver_Factory  tSolFactory;

            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
            //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            //        = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            //        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10000;
            //        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
            //        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
            //        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // create nonlinear solver and algorithm
            // --------------------------------------------------------------------------------------
            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolverMain;
            tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
            tNonlinearSolverMain.set_dof_type_list( tSolveDofTypes );

            // Create solver database
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
            tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

            // create time Solver and algorithm
            //------------------------------------------------------------------------------
            tsa::Time_Solver_Factory tTimeSolverFactory;

            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
            = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

            tsa::Time_Solver tTimeSolver;
            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_dof_type_list( tSolveDofTypes );
            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench2 );

            // solve
            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            //        // print full solution for debug
            //        Matrix<DDRMat> tFullSol;
            //        tTimeSolver.get_full_solution( tFullSol );
            //        print( tFullSol, "tFullSol" );

            // clean up
            //------------------------------------------------------------------------------
            delete tIntegMesh1;
            delete tModel;
            delete tInterpolationMesh;
        }
    }

    TEST_CASE("FEM Benchmark Diffusion Inclusion - 2Mat","[MDL_FEM_Benchmark_Diffusion_Inclusion]")
    {
        if(par_size()<=1)
        {
            // define problem parameters
            //------------------------------------------------------------------------------
            // geometry Parameters
            moris::real tDomainLX  = 10.0;                   /* Length of full domain in x (m) */
            moris::real tDomainLY  = 10.0;                   /* Length of full domain in y (m) */
            moris::real tRInner    = 0.55;                   /* Inner circle radius (m) */
            moris::real tRMiddle   = 0.77;                   /* Middle circle radius (m) */
            moris::real tROuter    = 1.05;                   /* Outer circle radius (m) */
            Matrix<DDRMat> tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/

            // material Parameters
            moris::real tKappaA = 1.0; /* Conductivity material A (W/m^2) */
            moris::real tKappaB = 100.0; /* Conductivity material B (W/m^2) */

            // boundary Conditions
            moris::real tTDirichlet = 5.0;   /* Imposed temperature for Dirichlet BC (K) */
            moris::real tDBCGamma   = 100.0; /* Penalty for Dirichlet BC */
            moris::real tQ          = 20.0;  /* Imposed heat flux for Neumann BC (W/m) */
            moris::real tIBCGamma   = 100.0;   /* Penalty for Interface BC */

            // mesh Setup
            moris::uint tNumX   = 20; /* Number of elements in x*/
            moris::uint tNumY   = 20; /* Number of elements in y*/

            // define hmr parameters
            //------------------------------------------------------------------------------
            uint tLagrangeMeshIndex = 0;
            std::string tOuterFieldName  = "OuterCircle";
            std::string tMiddleFieldName = "MiddleCircle";
            std::string tInnerFieldName  = "InnerCircle";
            ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
            tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
            tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) );
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
            tParameters.set( "initial_refinement", "1" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );
            tParameters.set( "use_number_aura", 0 );

            std::shared_ptr<hmr::HMR> tHMR = std::make_shared<hmr::HMR>( tParameters );

            // Initial refinement
            tHMR->perform_initial_refinement();

            // Create geometry engine
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(3);
            tGeometry(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
            tGeometry(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRMiddle);
            tGeometry(2) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);

            // Perform additional refinement
            //tGENGeometryEngine.perform_refinement(tHMR);

            // Get interpolation mesh
            tHMR->finalize();
            moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR->create_interpolation_mesh(tLagrangeMeshIndex);

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGENGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);

            // --------------------------------------------------------------------------------------
            xtk::Model tXTKModel( 2, tInterpolationMesh, &tGENGeometryEngine );
            tXTKModel.mVerbose = true;

            // specify decomposition method and cut mesh
            Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );
            //        tXTKModel.construct_face_oriented_ghost_penalization_cells();

            xtk::Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = false;
            tOutputOptions.mAddSideSets = true;
            tOutputOptions.mAddClusters = false;

            // output integration mesh
            moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh( tOutputOptions );
            std::string tOutputFile = "./mdl_exo/FEM_Bench2_2Mat.exo";
            tIntegMesh1->create_output_mesh( tOutputFile );

            // get meshes for FEM
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

            // create the properties
            //------------------------------------------------------------------------------
            std::shared_ptr< fem::Property > tPropKappaA = std::make_shared< fem::Property >();
            tPropKappaA->set_parameters( { {{ tKappaA }} } );
            tPropKappaA->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropKappaB = std::make_shared< fem::Property >();
            tPropKappaB->set_parameters( { {{ tKappaB }} } );
            tPropKappaB->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { {{ tTDirichlet }} } );
            tPropDirichlet->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { tQ } } } );
            tPropNeumann->set_val_function( ConstFunctionVal_MDLFEMBench2 );

            std::shared_ptr< fem::Property > tPropL2Analytic = std::make_shared< fem::Property >();
            tPropL2Analytic->set_parameters( { {{ tRInner }}, {{ tRMiddle }}, {{ tROuter }}, tCenterPoint, {{ tTDirichlet }}, {{ tQ }}, {{ tKappaA }}, {{ tKappaB }} } );
            tPropL2Analytic->set_val_function( AnalyticalTemp2MatFunc_MDLFEMBench2 );

            // create constitutive models
            //------------------------------------------------------------------------------
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIsoA = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIsoA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIsoA->set_property( tPropKappaA, "Conductivity" );
            tCMDiffLinIsoA->set_space_dim( 2 );
            tCMDiffLinIsoA->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIsoB = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIsoB->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMDiffLinIsoB->set_property( tPropKappaB, "Conductivity" );
            tCMDiffLinIsoB->set_space_dim( 2 );
            tCMDiffLinIsoB->set_local_properties();

            // create stabilization parameters
            //------------------------------------------------------------------------------
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                    tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { {{ tDBCGamma }} } );
            tSPDirichletNitsche->set_property( tPropKappaA, "Material", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                    tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
            tSPNitscheInterface->set_parameters( { {{ tIBCGamma }} } );
            tSPNitscheInterface->set_property( tPropKappaB, "Material", mtk::Master_Slave::MASTER );
            tSPNitscheInterface->set_property( tPropKappaA, "Material", mtk::Master_Slave::SLAVE );

            // create the IWGs
            //------------------------------------------------------------------------------
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulkA = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkA->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkA->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulkA->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGBulkB = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkB->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkB->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGBulkB->set_constitutive_model( tCMDiffLinIsoB, "Diffusion", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIsoB, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
            tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );
            tIWGInterface->set_constitutive_model( tCMDiffLinIsoB, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGInterface->set_constitutive_model( tCMDiffLinIsoA, "Diffusion", mtk::Master_Slave::SLAVE );
            tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

            // create the IQIs
            //------------------------------------------------------------------------------
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_Temp" );

            std::shared_ptr< fem::IQI > tIQIL2 = tIQIFactory.create_IQI( fem::IQI_Type::L2_ERROR_ANALYTIC );
            tIQIL2->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQIL2->set_property( tPropL2Analytic, "L2Check", mtk::Master_Slave::MASTER );
            tIQIL2->set_name( "IQI_L2" );

            std::shared_ptr< fem::IQI > tIQITempExact = tIQIFactory.create_IQI( fem::IQI_Type::PROPERTY );
            tIQITempExact->set_property( tPropL2Analytic, "Property", mtk::Master_Slave::MASTER );
            tIQITempExact->set_name( "IQI_Exact" );

            // create set info
            //------------------------------------------------------------------------------
            fem::Set_User_Info tSetBulkB1;
            tSetBulkB1.set_mesh_set_name( "HMR_dummy_c_p1" );
            tSetBulkB1.set_IWGs( { tIWGBulkB } );
            tSetBulkB1.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulkB2;
            tSetBulkB2.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulkB2.set_IWGs( { tIWGBulkB } );
            tSetBulkB2.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulkA1;
            tSetBulkA1.set_mesh_set_name( "HMR_dummy_c_p3" );
            tSetBulkA1.set_IWGs( { tIWGBulkA } );
            tSetBulkA1.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetBulkA2;
            tSetBulkA2.set_mesh_set_name( "HMR_dummy_n_p3" );
            tSetBulkA2.set_IWGs( { tIWGBulkA } );
            tSetBulkA2.set_IQIs( { tIQITEMP, tIQIL2, tIQITempExact } );

            fem::Set_User_Info tSetDirichletB1;
            tSetDirichletB1.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name( 2, 1, 0 ) );
            tSetDirichletB1.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumannA1;
            tSetNeumannA1.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name( 0, 3, 7 ) );
            tSetNeumannA1.set_IWGs( { tIWGNeumann } );

            fem::Set_User_Info tSetInterfaceBA1;
            tSetInterfaceBA1.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 1, 3 ) );
            tSetInterfaceBA1.set_IWGs( { tIWGInterface } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 7 );
            tSetInfo( 0 )  = tSetBulkB1;
            tSetInfo( 1 )  = tSetBulkB2;
            tSetInfo( 2 )  = tSetBulkA1;
            tSetInfo( 3 )  = tSetBulkA2;
            tSetInfo( 4 )  = tSetDirichletB1;
            tSetInfo( 5 )  = tSetNeumannA1;
            tSetInfo( 6 )  = tSetInterfaceBA1;

            // create model
            //------------------------------------------------------------------------------
            mdl::Model * tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0, false );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                    "./",
                    "UT_MDL_FEM_Bench2_Output_2Mat.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p1", "HMR_dummy_n_p1", "HMR_dummy_c_p3", "HMR_dummy_n_p3" },
                    { "TEMP", "L2", "TEMP_EXACT" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_Temp","IQI_L2","IQI_Exact" } );
            tModel->set_output_manager( &tOutputData );

            // create linear solver and algorithm
            // --------------------------------------------------------------------------------------
            // define dof type for solve
            moris::Cell< enum MSI::Dof_Type > tSolveDofTypes( 1 );
            tSolveDofTypes( 0 ) = MSI::Dof_Type::TEMP;

            dla::Solver_Factory  tSolFactory;

            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
            //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
            //        = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            //        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
            //        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10000;
            //        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
            //        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
            //        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 10;
            //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";

            dla::Linear_Solver tLinSolver;
            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // create nonlinear solver and algorithm
            // --------------------------------------------------------------------------------------
            NLA::Nonlinear_Solver_Factory tNonlinFactory;
            std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

            tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 3;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

            tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
            //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

            NLA::Nonlinear_Solver tNonlinearSolverMain;
            tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
            tNonlinearSolverMain.set_dof_type_list( tSolveDofTypes );

            // Create solver database
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
            tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

            // create time Solver and algorithm
            //------------------------------------------------------------------------------
            tsa::Time_Solver_Factory tTimeSolverFactory;

            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
            = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

            tsa::Time_Solver tTimeSolver;
            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_dof_type_list( tSolveDofTypes );
            tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFEMBench2 );

            // solve
            //------------------------------------------------------------------------------
            tTimeSolver.solve();

            //        // print full solution for debug
            //        Matrix<DDRMat> tFullSol;
            //        tTimeSolver.get_full_solution( tFullSol );
            //        print( tFullSol, "tFullSol" );

            // clean up
            //------------------------------------------------------------------------------
            delete tIntegMesh1;
            delete tModel;
            delete tInterpolationMesh;
        }
    }
}
