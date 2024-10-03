/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_Gyroid.cpp
 *
 */

#include "catch.hpp"
#include "cl_Star.hpp"
#include "cl_Circle.hpp"
#include "cl_Plane.hpp"

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

#include "cl_FEM_IWG_Factory.hpp"      //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"      //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"       //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"       //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"    //FEM/INT/src

#include "cl_MDL_Model.hpp"
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
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

#include "cl_GE_Geometry_Library.hpp"
#include "cl_GEN_Geom_Field_HMR.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

#include "fn_norm.hpp"

namespace moris
{

    void tPropValConstFunc( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&       aParameters,
            moris::fem::Field_Interpolator_Manager*         aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }
    void tPropValFunc( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&  aParameters,
            moris::fem::Field_Interpolator_Manager*    aFIManager )
    {
        aPropMatrix = aParameters( 0 ) + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
    }
    void tPropDerFunc( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&  aParameters,
            moris::fem::Field_Interpolator_Manager*    aFIManager )
    {
        aPropMatrix = aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
    }

    bool tSolverOutputCriteria_UTGyroid( moris::tsa::Time_Solver* )
    {
        return true;
    }

    int user_defined_refinement( hmr::Element* aElement,
            const Vector< Matrix< DDRMat > >&  aElementLocalValues,
            ParameterList&                     aParameters )
    {
        int aDoRefine = -1;

        // abs variable field threshold
        real lsth = 0.0;

        // abs variable field bandwidth (absolute)
        real lsbwabs = 0.4;

        // maximum refinement level
        uint maxlevel = 4;

        // min refinement level
        uint minlevel = 0;

        // max refinement level along interface
        uint maxifcref = 4;

        // max refinement level in volume
        uint maxvolref = 1;

        // current refinement level of element
        uint curlevel = aElement->get_level();

        // refinement strategy
        if ( aElementLocalValues( 0 ).max() >= lsth - lsbwabs )
        {
            // for volume refinement
            if ( aElementLocalValues( 0 ).min() >= lsth + lsbwabs )
            {
                if ( curlevel < maxvolref && curlevel < maxlevel )
                {
                    aDoRefine = 1;    // refine
                }
                else if ( curlevel == maxvolref || curlevel == minlevel )
                {
                    aDoRefine = 0;    // keep
                }
                else
                {
                    aDoRefine = -1;    // coarsen
                }
            }
            // for interface refinement
            else
            {
                if ( curlevel < maxifcref && curlevel < maxlevel )
                {
                    aDoRefine = 1;    // refine
                }
                else if ( curlevel == maxifcref || curlevel == minlevel )
                {
                    aDoRefine = 0;    // keep
                }
                else
                {
                    aDoRefine = -1;    // coarsen
                }
            }
        }
        else
        {
            if ( curlevel < minlevel )
            {
                aDoRefine = 1;    // refine
            }
            else if ( curlevel == minlevel )
            {
                aDoRefine = 0;    // keep
            }
        }

        return aDoRefine;
    }

    TEST_CASE( "MDL Gyroid", "[MDL_Gyroid]" )
    {

        //    if(par_size() == 1)
        //    {
        //	std::cout<<"I am proc: "<<par_rank()<<std::endl;

        uint tLagrangeMeshIndex = 0;
        // empty container for B-Spline meshes
        Vector< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

        // create settings object
        moris::hmr::Parameters tParameters;

        // Dummy parameter list
        ParameterList tParam = prm::create_hmr_parameter_list();

        tParameters.set_number_of_elements_per_dimension( { { 10 }, { 10 }, { 10 } } );
        tParameters.set_domain_dimensions( 5, 5, 5 );
        tParameters.set_domain_offset( 0.1, 0.1, 0.1 );
        tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

        tParameters.set_bspline_truncation( true );

        tParameters.set_lagrange_orders( { { 1 } } );
        tParameters.set_lagrange_patterns( { { 0 } } );

        tParameters.set_bspline_orders( { { 1 } } );
        tParameters.set_bspline_patterns( { { 0 } } );

        tParameters.set_output_meshes( { { 0 } } );
        tParameters.set_lagrange_input_mesh( { { 0 } } );

        tParameters.set_staircase_buffer( 1 );

        tParameters.set_initial_refinement( { { 2 } } );
        tParameters.set_initial_refinement_patterns( { { 0 } } );

        tParameters.set_number_aura( true );


        // create the HMR object by passing the settings to the constructor
        moris::hmr::HMR tHMR( tParameters );

        // std::shared_ptr< Database >
        //        auto tDatabase = tHMR.get_database();

        tHMR.perform_initial_refinement();

        std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );    // HMR Lagrange mesh
        //==============================
        std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex );

        tField->evaluate_scalar_function( moris::gen::getDistanceToGyroidsMassive );

        Vector< std::shared_ptr< moris::hmr::Field > > tFields( 1, tField );

        for ( uint k = 0; k < 1; ++k )
        {
            //            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.user_defined_flagging( user_defined_refinement, tFields, tParam, 0 );
            tHMR.perform_refinement_based_on_working_pattern( 0, true );
            tField->evaluate_scalar_function( moris::gen::getDistanceToGyroidsMassive );
        }
        tHMR.finalize();

        //        tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_initial_3x3x3.vtk");

        //==============================
        tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

        hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        moris::gen::GEN_Geom_Field_HMR tFieldAsGeom( tField );

        Vector< moris::gen::GEN_Geometry* > tGeometryVector = { &tFieldAsGeom };

        size_t                      tModelDimension = 3;
        moris::gen::Geometry_Engine tGeometryEngine( tGeometryVector, tModelDimension );

        //        moris::gen::Geometry_Engine tGeometryEngine;

        xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
        tXTKModel.mVerbose = false;

        // Specify decomposition Method and Cut Mesh ---------------------------------------
        Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
        tXTKModel.decompose( tDecompositionMethods );

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

        xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        //==============================

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
        tPropConductivity1->set_parameters( { { { 1.0 } } } );
        tPropConductivity1->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
        tPropConductivity2->set_parameters( { { { 5.0 } } } );
        tPropConductivity2->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { { { 5.0 } } } );
        tPropDirichlet->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { { { 20.0 } } } );
        tPropNeumann->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { { { 100.0 } } } );
        tPropTempLoad1->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
        tPropTempLoad2->set_parameters( { { { 100.0 } } } );
        tPropTempLoad2->set_val_function( tPropValConstFunc );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
        tCMDiffLinIso1->set_space_dim( 3 );
        tCMDiffLinIso1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
        tCMDiffLinIso2->set_space_dim( 3 );
        tCMDiffLinIso2->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory                                 tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
        tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::FOLLOWER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk1 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGBulk2 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk2->set_property( tPropTempLoad2, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGInterface =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
        tIWGInterface->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        // create the IQIs
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQITEMP->set_output_type( vis::Output_Type::TEMP );
        tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIQITEMP->set_output_type_index( 0 );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
        tSetBulk1.set_IWGs( { tIWGBulk2 } );
        tSetBulk1.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
        tSetBulk2.set_IWGs( { tIWGBulk2 } );
        tSetBulk2.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk3;
        tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
        tSetBulk3.set_IWGs( { tIWGBulk1 } );
        tSetBulk3.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk4.set_IWGs( { tIWGBulk1 } );
        tSetBulk4.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_set_name( "SideSet_2_n_p1" );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );

        fem::Set_User_Info tSetInterface1;
        tSetInterface1.set_mesh_set_name( tEnrIntegMesh.get_set_index_by_name( tDblInterfaceSideSetName ) );
        tSetInterface1.set_IWGs( { tIWGInterface } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 7 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetBulk3;
        tSetInfo( 3 ) = tSetBulk4;
        tSetInfo( 4 ) = tSetDirichlet;
        tSetInfo( 5 ) = tSetNeumann;
        tSetInfo( 6 ) = tSetInterface1;

        // create model
        mdl::Model* tModel = new mdl::Model( tMeshManager,
                0,
                tSetInfo );

        // --------------------------------------------------------------------------------------
        // Define outputs
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                vis::VIS_Mesh_Type::STANDARD,
                "UT_Gyroid_Output_01.exo",
                "./",
                "temp.exo",
                { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                { 0, 1, 2, 3 },
                { "Temperature" },
                { vis::Field_Type::NODAL },
                { vis::Output_Type::TEMP } );

        tModel->set_output_manager( &tOutputData );

        dla::Solver_Factory                             tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

        //       tLinearSolverParameterList.set( "AZ_diagnostics", AZ_all );
        //       tLinearSolverParameterList.set( "AZ_output", AZ_all );
        //       tLinearSolverParameterList.set( "AZ_solver", AZ_gmres_condnum );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory               tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory                      tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        //       tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        //       tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
        tNonlinearSolver.set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tTimeSolver.set_dof_type_list( { { MSI::Dof_Type::TEMP } } );

        tTimeSolver.set_output( 0, tSolverOutputCriteria_UTGyroid );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tTimeSolver.solve();
        Matrix< DDRMat > tFullSol;
        tTimeSolver.get_full_solution( tFullSol );

        delete tInterpMesh;
        //    }
    } /* END_TEST_CASE */

}    // namespace moris
