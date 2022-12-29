/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_DiffusionElement.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "typedefs.hpp"

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
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_Matrix.hpp"    //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"    // ALG/src

#include "cl_FEM_NodeProxy.hpp"          //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"       //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"          //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"    //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"        //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"      //FEM/INT/src

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

#include "fn_norm.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"

namespace moris
{
    inline moris::real
    LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat >& aPoint )
    {

        real mXn = 0;
        real mYn = 0;
        real mZn = 1.0;
        real mXc = 1.011;
        real mYc = 1.011;
        real mZc = 1.411;
        return mXn * ( aPoint( 0 ) - mXc ) + mYn * ( aPoint( 1 ) - mYc ) + mZn * ( aPoint( 2 ) - mZc );
    }

    inline moris::real
    LevelSetSphereFunction( const moris::Matrix< moris::DDRMat >& aPoint )
    {

        real mR  = 0.77;
        real mXc = 0.0;
        real mYc = 0.0;
        real mZc = 0.0;
        return ( aPoint( 0 ) - mXc ) * ( aPoint( 0 ) - mXc )
             + ( aPoint( 1 ) - mYc ) * ( aPoint( 1 ) - mYc )
             + ( aPoint( 2 ) - mZc ) * ( aPoint( 2 ) - mZc )
             - ( mR * mR );
    }

    inline moris::real
    LevelSetSphereCylinder( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::Matrix< moris::DDRMat > aCenter = { { 0.0 }, { 0.0 }, { 0.0 } };
        moris::Matrix< moris::DDRMat > aAxis   = { { 0.0 }, { 1.0 }, { 0.0 } };
        moris::real                    aRad    = 0.77;
        moris::real                    aLength = 5;

        MORIS_ASSERT( aCenter.numel() == 3, "Centers need to have length 3" );
        MORIS_ASSERT( aAxis.numel() == 3, "axis need to have length 3" );

        Cell< moris::real > relativePosition = { ( aPoint( 0 ) - aCenter( 0 ) ), ( aPoint( 1 ) - aCenter( 1 ) ), ( aPoint( 2 ) - aCenter( 2 ) ) };
        moris::real         lsFromLeft       = ( relativePosition( 0 ) * ( -aAxis( 0 ) ) + relativePosition( 1 ) * ( -aAxis( 1 ) ) + relativePosition( 2 ) * ( -aAxis( 2 ) ) ) - aLength / 2.0;
        moris::real         lsFromRight      = ( relativePosition( 0 ) * ( aAxis( 0 ) ) + relativePosition( 1 ) * ( aAxis( 1 ) ) + relativePosition( 2 ) * ( aAxis( 2 ) ) ) - aLength / 2.0;

        moris::real         axialCrd  = ( relativePosition( 0 ) * ( aAxis( 0 ) ) + relativePosition( 1 ) * ( aAxis( 1 ) ) + relativePosition( 2 ) * ( aAxis( 2 ) ) );
        Cell< moris::real > radDir    = { ( relativePosition( 0 ) - aAxis( 0 ) * axialCrd ), ( relativePosition( 1 ) - aAxis( 1 ) * axialCrd ), ( relativePosition( 2 ) - aAxis( 2 ) * axialCrd ) };
        moris::real         radDist   = std::pow( radDir( 0 ) * radDir( 0 ) + radDir( 1 ) * radDir( 1 ) + radDir( 2 ) * radDir( 2 ), 0.5 );
        moris::real         lsFromRad = radDist - aRad;

        return -std::max( std::max( lsFromLeft, lsFromRight ), lsFromRad );
    }

    inline moris::real
    LevelSetSphereCylinderGeometry( const moris::Matrix< moris::DDRMat >& aCoordinates, const moris::Cell< moris::real* >& aParameters )
    {
        return LevelSetSphereCylinder( aCoordinates );
    }

    inline moris::real
    LevelSetFunction_star( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real tPhi = std::atan2( aPoint( 0 ), aPoint( 2 ) );

        moris::real tLevelSetVaue = 0.5 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aPoint( 0 ), 2 ) + std::pow( aPoint( 2 ), 2 ) );

        return tLevelSetVaue;
    }

    inline void
    tConstValFunction_MDL_XTK_HMR(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    inline bool
    tSolverOutputCriteria( moris::tsa::Time_Solver* )
    {
        return true;
    }

    TEST_CASE( "HMR Interpolation STK Cut Diffusion Model Lag Order 2", "[XTK_HMR_STK_DIFF]" )
    {
        if ( par_size() == 1 )
        {
            std::string tFieldName = "Circle";

            moris::uint tLagrangeMeshIndex = 0;
            moris::uint tBSplineMeshIndex  = 0;

            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 4 }, { 4 } } );
            tParameters.set_domain_dimensions( { { 1 }, { 1 }, { 2 } } );
            tParameters.set_domain_offset( { { 0.0 }, { 0.0 }, { 0.0 } } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_side_sets( { { 5 }, { 6 } } );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_union_pattern( 2 );
            tParameters.set_working_pattern( 3 );

            tParameters.set_refinement_buffer( 2 );
            tParameters.set_staircase_buffer( 2 );

            Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LevelSetPlaneFunction );

            // FIXME what is the following test about
            if ( false )
            {
                for ( uint k = 0; k < 1; ++k )
                {
                    tHMR.flag_surface_elements_on_working_pattern( tField );
                    tHMR.perform_refinement_based_on_working_pattern( 0 );

                    tField->evaluate_scalar_function( LevelSetPlaneFunction );
                }
            }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_bar_hole_interp_l1_b1.e" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            moris::Cell< std::shared_ptr< moris::ge::Geometry > > tGeometryVector( 1 );
            tGeometryVector( 0 ) = std::make_shared< moris::ge::Plane >( 1.011, 1.011, 1.411, 0.0, 0.0, 1.0 );

            // Tell the geometry engine about the discrete field mesh and how to interpret phases
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::ge::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
            size_t                          tModelDimension       = 3;
            Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            xtk::Model                      tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Do the cutting
            tXTKModel.decompose( tDecompositionMethods );

            // Perform the enrichment
            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );
            //        tXTKModel.construct_face_oriented_ghost_penalization_cells();

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
            tPropConductivity->set_parameters( { { { 1.0 } } } );
            tPropConductivity->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 5.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { 20.0 } } } );
            tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
            tPropTempLoad->set_parameters( { { { 0.0 } } } );
            tPropTempLoad->set_val_function( tConstValFunction_MDL_XTK_HMR );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );    // FIXME through the factory?
            tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
            tCMDiffLinIso->set_space_dim( 3 );
            tCMDiffLinIso->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Master_Slave::MASTER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk } );

            fem::Set_User_Info tSetDirichlet;
            std::string        tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
            tSetDirichlet.set_mesh_set_name( tInterfaceSideSetName );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_1_n_p0" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet;
            tSetInfo( 3 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    tBSplineMeshIndex,
                    tSetInfo );

            moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            dla::Solver_Factory                             tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param( "AZ_diagnostics" ) = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_output" )      = AZ_none;

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
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
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( tDofTypes1 );
            tTimeSolver.set_dof_type_list( tDofTypes1 );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tTimeSolver.solve();

            // TODO: add gold solution data for this problem
            // verify solution
            Matrix< DDRMat > tFullSol;
            tTimeSolver.get_full_solution( tFullSol );

            // clean up
            delete tModel;
            delete tInterpMesh;
        }
    }

    TEST_CASE( "HMR Interpolation XTK Cut Diffusion Model Lag Order 2", "[XTK_HMR_DIFF]" )
    {
        if ( par_size() == 1 )
        {
            std::string tFieldName = "Cylinder";

            moris::uint tLagrangeMeshIndex = 0;
            moris::uint tBSplineMeshIndex  = 0;

            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 4 }, { 4 }, { 4 } } );
            tParameters.set_domain_dimensions( { { 1 }, { 1 }, { 2 } } );
            tParameters.set_domain_offset( { { 0.0 }, { 0.0 }, { 0.0 } } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_side_sets( { { 5 }, { 6 } } );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_union_pattern( 2 );
            tParameters.set_working_pattern( 3 );

            tParameters.set_refinement_buffer( 2 );
            tParameters.set_staircase_buffer( 2 );

            tParameters.set_number_aura( true );

            Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LevelSetPlaneFunction );

            for ( uint k = 0; k < 1; ++k )
            {
                //            tHMR.finalize();
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LevelSetPlaneFunction );
            }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./mdl_exo/xtk_hmr_bar_plane_interp_l2_b2.e" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            moris::Cell< std::shared_ptr< moris::ge::Geometry > > tGeometryVector( 1 );
            tGeometryVector( 0 ) = std::make_shared< moris::ge::Plane >( 1.011, 1.011, 1.411, 0.0, 0.0, 1.0 );

            // Tell the geometry engine about the discrete field mesh and how to interpret phases
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::ge::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
            size_t                          tModelDimension       = 3;
            Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            xtk::Model                      tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Do the cutting
            tXTKModel.decompose( tDecompositionMethods );

            // Perform the enrichment
            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
            tPropConductivity->set_parameters( { { { 1.0 } } } );
            tPropConductivity->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 5.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { 20.0 } } } );
            tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropTempLoad = std::make_shared< fem::Property >();
            tPropTempLoad->set_parameters( { { { 0.0 } } } );
            tPropTempLoad->set_val_function( tConstValFunction_MDL_XTK_HMR );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );    // FIXME through the factory?
            tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
            tCMDiffLinIso->set_space_dim( 3 );
            tCMDiffLinIso->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity, "Material", mtk::Master_Slave::MASTER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGBulk->set_property( tPropTempLoad, "Load", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            // define the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_Temp" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk } );
            tSetBulk1.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk } );
            tSetBulk2.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetDirichlet;
            std::string        tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
            tSetDirichlet.set_mesh_set_name( tInterfaceSideSetName );
            tSetDirichlet.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "SideSet_1_n_p0" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet;
            tSetInfo( 3 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    tBSplineMeshIndex,
                    tSetInfo );

            // --------------------------------------------------------------------------------------
            // Define outputs

            vis::Output_Manager tOutputData;

            tOutputData.set_outputs( 0,
                    //                                         VIS_Mesh_Type::STANDARD,
                    vis::VIS_Mesh_Type::OVERLAPPING_INTERFACE,
                    "./",
                    "XTK_HMR_DIFF.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_c_p0", "HMR_dummy_n_p0" },
                    { "Temperature" },
                    { vis::Field_Type::NODAL },
                    { "IQI_Temp" } );

            tModel->set_output_manager( &tOutputData );

            // --------------------------------------------------------------------------------------

            moris::Cell< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 1: create linear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            dla::Solver_Factory                             tSolFactory;
            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );

            tLinearSolverAlgorithm->set_param( "AZ_diagnostics" ) = AZ_none;
            tLinearSolverAlgorithm->set_param( "AZ_output" )      = AZ_none;

            dla::Linear_Solver tLinSolver;

            tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 2: create nonlinear solver and algorithm
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            NLA::Nonlinear_Solver_Factory               tNonlinFactory;
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
            tsa::Time_Solver_Factory                      tTimeSolverFactory;
            std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

            tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

            tsa::Time_Solver tTimeSolver;

            tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

            sol::SOL_Warehouse tSolverWarehouse;

            tSolverWarehouse.set_solver_interface( tModel->get_solver_interface() );

            tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
            tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

            tNonlinearSolver.set_dof_type_list( tDofTypes1 );
            tTimeSolver.set_dof_type_list( tDofTypes1 );

            tTimeSolver.set_output( 0, tSolverOutputCriteria );

            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            // STEP 4: Solve and check
            // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            tTimeSolver.solve();

            // TODO: add gold solution data for this problem
            // verify solution

            // clean up
            delete tModel;
            delete tInterpMesh;
        }
    }

    TEST_CASE( "HMR Interpolation XTK Cut Diffusion Model Multigrid", "[XTK_HMR_DIFF_MULTIGRID]" )
    {
        if ( par_size() == 1 )
        {
            gLogger.set_severity_level( 0 );

            std::string tFieldName = "Cylinder";

            // start timer
            tic tTimer_HMR;

            moris::uint tLagrangeMeshIndex = 0;

            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 }, { 4 } } );
            tParameters.set_domain_dimensions( { { 2 }, { 2 }, { 4 } } );
            tParameters.set_domain_offset( { { -1.0 }, { -1.0 }, { -2.0 } } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_side_sets( { { 5 }, { 6 } } );

            tParameters.set_multigrid( true );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_union_pattern( 2 );
            tParameters.set_working_pattern( 3 );

            tParameters.set_refinement_buffer( 2 );
            tParameters.set_staircase_buffer( 2 );

            Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LevelSetSphereCylinder );

            for ( uint k = 0; k < 2; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LevelSetSphereCylinder );
            }

            tHMR.finalize();

            // stop timer
            real tElapsedTime = tTimer_HMR.toc< moris::chronos::milliseconds >().wall;

            MORIS_LOG_INFO( " HMR took %5.3f seconds.", (double)tElapsedTime / 1000 );

            //        tHMR.save_to_exodus( "./mdl_exo/xtk_hmr_bar_hole_interp_l1_b1.e" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            // start timer
            tic tTimer_XTK;

            moris::Cell< std::shared_ptr< moris::ge::Geometry > > tGeometryVector( 1 );
            tGeometryVector( 0 ) = std::make_shared< moris::ge::User_Defined_Geometry >( Matrix< DDRMat >( 0, 0 ), &( LevelSetSphereCylinderGeometry ) );

            // Tell the geometry engine about the discrete field mesh and how to interpret phases
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::ge::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
            size_t                          tModelDimension       = 3;
            Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            xtk::Model                      tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Do the cutting
            tXTKModel.decompose( tDecompositionMethods );

            // Perform the enrichment
            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );

            tXTKModel.construct_multigrid();

            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
            tPropConductivity1->set_parameters( { { { 1.0 } } } );
            tPropConductivity1->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
            tPropConductivity2->set_parameters( { { { 1.0 } } } );
            tPropConductivity2->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 5.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { 20.0 } } } );
            tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
            tPropTempLoad1->set_parameters( { { { 0.0 } } } );
            tPropTempLoad1->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
            tPropTempLoad2->set_parameters( { { { 0.0 } } } );
            tPropTempLoad2->set_val_function( tConstValFunction_MDL_XTK_HMR );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 3 );
            tCMDiffLinIso1->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Master_Slave::MASTER );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            // define the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Master_Slave::MASTER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_Temp" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );
            tSetBulk1.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk1 } );
            tSetBulk2.set_IQIs( { tIQITEMP } );

            fem::Set_User_Info tSetDirichlet1;
            tSetDirichlet1.set_mesh_set_name( "SideSet_1_n_p0" );
            tSetDirichlet1.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetDirichlet2;
            tSetDirichlet2.set_mesh_set_name( "SideSet_2_n_p0" );
            tSetDirichlet2.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "iside_b0_0_b1_1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 5 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet1;
            tSetInfo( 3 ) = tSetDirichlet2;
            tSetInfo( 4 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0,
                    true );

            // --------------------------------------------------------------------------------------
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_Multigrid.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy" },
                    { "Temperature" },
                    { vis::Field_Type::NODAL },
                    { "IQI_Temp" } );
            tModel->set_output_manager( &tOutputData );

            //-------------------------------------------------------
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

            moris::Cell< moris::Cell< moris::ParameterList > > tParameterlist( 7 );
            for ( uint Ik = 0; Ik < 7; Ik++ )
            {
                tParameterlist( Ik ).resize( 1 );
            }

            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC );
            tParameterlist( 0 )( 0 ).set( "KSPType", std::string( "fgmres" ) );
            tParameterlist( 0 )( 0 ).set( "PCType", std::string( "mg" ) );
            tParameterlist( 0 )( 0 ).set( "ILUFill", 3 );
            tParameterlist( 0 )( 0 ).set( "ILUTol", 1e-6 );

            tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Crteria", "" );

            tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
            tParameterlist( 6 )( 0 ).set( "SOL_TPL_Type", static_cast< uint >( sol::MapType::Petsc ) );

            tSolverWarehouse.set_parameterlist( tParameterlist );

            tSolverWarehouse.initialize();

            tsa::Time_Solver* tTimeSolver = tSolverWarehouse.get_main_time_solver();

            tTimeSolver->solve();

            moris::Matrix< DDRMat > tSolution;
            tTimeSolver->get_full_solution( tSolution );

            CHECK( equal_to( tSolution( 0, 0 ), 4.999702007294005e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 1, 0 ), 5.000060374039293e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 2, 0 ), 4.999702006950493e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 3, 0 ), 5.000060374198882e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 4, 0 ), 1.712285753707913e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 5, 0 ), 1.707413784393956e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 6, 0 ), 1.712285753672601e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 7, 0 ), 1.707413784389796e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 8, 0 ), 5.000177244537252e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 382, 0 ), 5.149916424143507e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 461, 0 ), 1.712285752727456e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 505, 0 ), 2.947091742301102e+01, 1.0e+08 ) );

            // clean up
            delete tModel;
            delete tInterpMesh;
        }
    }

    TEST_CASE( "XTK Diffusion Multigrid", "[XTK_DIFF_MULTIGRID]" )
    {
        if ( par_size() == 1 )
        {
            gLogger.set_severity_level( 0 );

            std::string tFieldName = "Cylinder";

            // start timer
            tic tTimer_HMR;

            moris::uint tLagrangeMeshIndex = 0;

            moris::hmr::Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 }, { 4 } } );
            tParameters.set_domain_dimensions( { { 2 }, { 2 }, { 4 } } );
            tParameters.set_domain_offset( { { -1.0 }, { -1.0 }, { -2.0 } } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_side_sets( { { 5 }, { 6 } } );

            tParameters.set_multigrid( true );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_union_pattern( 2 );
            tParameters.set_working_pattern( 3 );

            tParameters.set_refinement_buffer( 2 );
            tParameters.set_staircase_buffer( 2 );

            Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            hmr::HMR tHMR( tParameters );

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LevelSetSphereCylinder );

            for ( uint k = 0; k < 2; ++k )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LevelSetSphereCylinder );
            }

            tHMR.finalize();

            // stop timer
            real tElapsedTime = tTimer_HMR.toc< moris::chronos::milliseconds >().wall;

            MORIS_LOG_INFO( " HMR took %5.3f seconds.", (double)tElapsedTime / 1000 );

            //        tHMR.save_to_exodus( "./mdl_exo/xtk_hmr_bar_hole_interp_l1_b1.e" );

            hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            // start timer
            tic tTimer_XTK;

            moris::Cell< std::shared_ptr< moris::ge::Geometry > > tGeometryVector( 1 );
            tGeometryVector( 0 ) = std::make_shared< moris::ge::User_Defined_Geometry >( Matrix< DDRMat >( 0, 0 ), LevelSetSphereCylinderGeometry );

            // Tell the geometry engine about the discrete field mesh and how to interpret phases
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::ge::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
            size_t                          tModelDimension       = 3;
            Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            xtk::Model                      tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            // Do the cutting
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );

            tXTKModel.construct_multigrid();

            // get meshes
            xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            //------------------------------------------------------------------------------
            // create the properties
            std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
            tPropConductivity1->set_parameters( { { { 1.0 } } } );
            tPropConductivity1->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
            tPropConductivity2->set_parameters( { { { 1.0 } } } );
            tPropConductivity2->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
            tPropDirichlet->set_parameters( { { { 5.0 } } } );
            tPropDirichlet->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
            tPropNeumann->set_parameters( { { { 20.0 } } } );
            tPropNeumann->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
            tPropTempLoad1->set_parameters( { { { 0.0 } } } );
            tPropTempLoad1->set_val_function( tConstValFunction_MDL_XTK_HMR );

            std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
            tPropTempLoad2->set_parameters( { { { 0.0 } } } );
            tPropTempLoad2->set_val_function( tConstValFunction_MDL_XTK_HMR );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 3 );
            tCMDiffLinIso1->set_local_properties();

            // define stabilization parameters
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPReciprocalVolume = tSPFactory.create_SP( fem::Stabilization_Type::RECIPROCAL_TOTAL_VOLUME );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
            tIWGDirichlet->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Master_Slave::MASTER );
            tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

            std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
            tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

            // define the IQIs
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIVolFraction = tIQIFactory.create_IQI( fem::IQI_Type::VOLUME_FRACTION );
            tIQIVolFraction->set_stabilization_parameter( tSPReciprocalVolume, "Reciprocal_total_vol" );

            // define set info
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_n_p0" );
            tSetBulk1.set_IWGs( { tIWGBulk1 } );
            tSetBulk1.set_IQIs( { tIQIVolFraction } );

            fem::Set_User_Info tSetBulk2;
            tSetBulk2.set_mesh_set_name( "HMR_dummy_c_p0" );
            tSetBulk2.set_IWGs( { tIWGBulk1 } );
            tSetBulk2.set_IQIs( { tIQIVolFraction } );

            fem::Set_User_Info tSetDirichlet1;
            tSetDirichlet1.set_mesh_set_name( "SideSet_1_n_p0" );
            tSetDirichlet1.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetDirichlet2;
            tSetDirichlet2.set_mesh_set_name( "SideSet_2_n_p0" );
            tSetDirichlet2.set_IWGs( { tIWGDirichlet } );

            fem::Set_User_Info tSetNeumann;
            tSetNeumann.set_mesh_set_name( "iside_b0_0_b1_1" );
            tSetNeumann.set_IWGs( { tIWGNeumann } );

            // create a cell of set info
            moris::Cell< fem::Set_User_Info > tSetInfo( 5 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetBulk2;
            tSetInfo( 2 ) = tSetDirichlet1;
            tSetInfo( 3 ) = tSetDirichlet2;
            tSetInfo( 4 ) = tSetNeumann;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0,
                    true );

            // --------------------------------------------------------------------------------------
            // define outputs
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,
                    "./",
                    "UT_MDL_Multigrid.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy" },
                    { "Temperature" },
                    { vis::Field_Type::NODAL },
                    { "IQI_Temp" } );
            tModel->set_output_manager( &tOutputData );

            //-------------------------------------------------------
            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

            moris::Cell< moris::Cell< moris::ParameterList > > tParameterlist( 7 );
            for ( uint Ik = 0; Ik < 7; Ik++ )
            {
                tParameterlist( Ik ).resize( 1 );
            }

            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC );
            tParameterlist( 0 )( 0 ).set( "KSPType", std::string( "fgmres" ) );
            tParameterlist( 0 )( 0 ).set( "PCType", std::string( "mg" ) );
            tParameterlist( 0 )( 0 ).set( "ILUFill", 0 );
            tParameterlist( 0 )( 0 ).set( "ILUTol", 1e-6 );

            tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "TEMP" );

            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "TEMP" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "" );
            tParameterlist( 5 )( 0 ).set( "TSA_Output_Crteria", "" );

            tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
            tParameterlist( 6 )( 0 ).set( "SOL_TPL_Type", static_cast< uint >( sol::MapType::Petsc ) );

            tSolverWarehouse.set_parameterlist( tParameterlist );

            tSolverWarehouse.initialize();

            tsa::Time_Solver* tTimeSolver = tSolverWarehouse.get_main_time_solver();

            tTimeSolver->solve();

            moris::Matrix< DDRMat > tSolution;
            tTimeSolver->get_full_solution( tSolution );

            CHECK( equal_to( tSolution( 0, 0 ), 4.999702007294010e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 1, 0 ), 5.000060374039292e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 2, 0 ), 4.999702006950497e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 3, 0 ), 5.000060374198884e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 4, 0 ), 1.712285753707912e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 5, 0 ), 1.707413784393955e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 6, 0 ), 1.712285753672600e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 7, 0 ), 1.707413784389795e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 8, 0 ), 5.000177244537259e+00, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 382, 0 ), 5.149916424143507e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 461, 0 ), 1.712285752727455e+01, 1.0e+08 ) );
            CHECK( equal_to( tSolution( 505, 0 ), 2.947091742301103e+01, 1.0e+08 ) );

            delete tInterpMesh;
        }
    }

}    // namespace moris
