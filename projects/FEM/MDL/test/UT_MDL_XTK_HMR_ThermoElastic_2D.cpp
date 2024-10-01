/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_ThermoElastic_2D.cpp
 *
 */

#include "catch.hpp"
#include "paths.hpp"
#include <memory>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "moris_typedefs.hpp"

#include "HDF5_Tools.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster_Input.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Side_Cluster_Input.hpp"

#include "cl_Matrix.hpp"    //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"    // ALG/src

#include "cl_FEM_NodeProxy.hpp"                     //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"                  //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                     //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"               //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"                   //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"                   //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"                    //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"                 //FEM/INT/src
#include "cl_FEM_Property.hpp"                      //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"    //FEM/INT/src
// FEM/VIS/src
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"

#include "cl_MDL_Model.hpp"

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

#include "cl_GEN_Line.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris
{

    inline moris::real
    LvlSetCircle_2D( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        return std::sqrt( aPoint( 0 ) * aPoint( 0 ) + aPoint( 1 ) * aPoint( 1 ) ) - 0.2505;
    }

    inline moris::real
    LvlSetPlane( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        return aPoint( 0 ) + 500;
    }

    inline moris::real
    LevelSetFunction_star1( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real tPhi          = std::atan2( aPoint( 0 ), aPoint( 1 ) );
        moris::real tLevelSetVaue = 0.501 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aPoint( 0 ), 2 ) + std::pow( aPoint( 1 ), 2 ) );
        return -tLevelSetVaue;
    }

    inline moris::real
    Plane4MatMDL1( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real mXC = 0.1;
        moris::real mYC = 0.1;
        moris::real mNx = 1.0;
        moris::real mNy = 0.0;
        return ( mNx * ( aPoint( 0 ) - mXC ) + mNy * ( aPoint( 1 ) - mYC ) );
    }

    inline moris::real
    Circle4MatMDL( const moris::Matrix< moris::DDRMat >& aPoint )
    {
        moris::real mXCenter = 0.01;
        moris::real mYCenter = 0.01;
        moris::real mRadius  = 0.47334;

        return ( aPoint( 0 ) - mXCenter ) * ( aPoint( 0 ) - mXCenter )
             + ( aPoint( 1 ) - mYCenter ) * ( aPoint( 1 ) - mYCenter )
             - ( mRadius * mRadius );
    }

    inline void
    tConstValFunction(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    inline bool
    tSolverOutputCriteria_thermolast( moris::tsa::Time_Solver* )
    {
        return true;
    }

    TEST_CASE( "2D XTK WITH HMR ThermoElastic 2D", "[XTK_HMR_thermoelastic_2D]" )
    {
        if ( par_size() <= 1 )
        {
            uint        tLagrangeMeshIndex = 0;
            std::string tFieldName         = "Cylinder";

            Parameter_List tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::string( "2, 1" ) );
            tParameters.set( "domain_dimensions", "2, 2" );
            tParameters.set( "domain_offset", "-1.0, -1.0" );
            tParameters.set( "lagrange_output_meshes", std::string( "0" ) );

            tParameters.set( "lagrange_orders", std::string( "1" ) );
            tParameters.set( "lagrange_pattern", std::string( "0" ) );
            tParameters.set( "bspline_orders", std::string( "1" ) );
            tParameters.set( "bspline_pattern", std::string( "0" ) );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 1 );
            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            //// create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LvlSetPlane );
            //
            // for( uint k=0; k<2; ++k )
            // {
            // tHMR.flag_surface_elements_on_working_pattern( tField );
            // tHMR.perform_refinement_based_on_working_pattern( 0 );

            // tField->evaluate_scalar_function( LvlSetCircle_2D );
            // }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto                                              tPlane          = std::make_shared< moris::gen::Line >( -500.0, 0.0, 1.0, 0.0 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

            size_t tModelDimension = 2;
            //------------------------------------------------------------------------------
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( tModelDimension, tInterpolationMesh, &tGeometryEngine );

            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

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
            tPropConductivity1->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
            tPropConductivity2->set_parameters( { { { 5.0 } } } );
            tPropConductivity2->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichletU = std::make_shared< fem::Property >();
            tPropDirichletU->set_parameters( { { { 0.0 }, { 0.0 } } } );
            tPropDirichletU->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichletTEMP = std::make_shared< fem::Property >();
            tPropDirichletTEMP->set_parameters( { { { 3.0 } } } );
            tPropDirichletTEMP->set_val_function( tConstValFunction );

            //        std::shared_ptr< fem::Property > tPropNeumannTEMP = std::make_shared< fem::Property >();
            //        tPropNeumann->set_parameters( { {{ 20.0 }} } );
            //        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

            std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
            tPropEMod1->set_parameters( { { { 1.0 } } } );
            tPropEMod1->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
            tPropEMod2->set_parameters( { { { 1.0 } } } );
            tPropEMod2->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropPoisson = std::make_shared< fem::Property >();
            tPropPoisson->set_parameters( { { { 0.0 } } } );
            tPropPoisson->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tCTE = std::make_shared< fem::Property >();
            tCTE->set_parameters( { { { 1.0 } } } );
            tCTE->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tTRef = std::make_shared< fem::Property >();
            tTRef->set_parameters( { { { 1.0 } } } );
            tTRef->set_val_function( tConstValFunction );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY }, { MSI::Dof_Type::TEMP } },
                    { "Displacement", "Temperature" } );
            tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropPoisson, "PoissonRatio" );
            tCMStrucLinIso1->set_property( tCTE, "CTE" );
            tCMStrucLinIso1->set_property( tTRef, "ReferenceTemperature" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso1->set_space_dim( 2 );
            tCMStrucLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY }, { MSI::Dof_Type::TEMP } },
                    { "Displacement", "Temperature" } );
            tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
            tCMStrucLinIso2->set_property( tPropPoisson, "PoissonRatio" );
            tCMStrucLinIso2->set_property( tCTE, "CTE" );
            tCMStrucLinIso2->set_property( tTRef, "ReferenceTemperature" );
            tCMStrucLinIso2->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso2->set_space_dim( 2 );
            tCMStrucLinIso2->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 2 );
            tCMDiffLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
            tCMDiffLinIso2->set_space_dim( 2 );
            tCMDiffLinIso2->set_local_properties();

            //----------------------------------------------------------------------------------------------------------
            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheU = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitscheU->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitscheU->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheTEMP = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitscheTEMP->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitscheTEMP->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );

            //----------------------------------------------------------------------------------------------------------
            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulkU_1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulkU_1->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulkU_2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulkU_2->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichletU = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichletU->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichletU->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichletU->set_stabilization_parameter( tSPDirichletNitscheU, "DirichletNitsche" );
            tIWGDirichletU->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichletU->set_property( tPropDirichletU, "Dirichlet", mtk::Leader_Follower::LEADER );

            //        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            //        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            //        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
            //        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulkTEMP_1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkTEMP_1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            //        tIWGBulkTEMP_1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulkTEMP_2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkTEMP_2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            //        tIWGBulkTEMP_2->set_property( tPropTempLoad2, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichletTEMP = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichletTEMP->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichletTEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichletTEMP->set_stabilization_parameter( tSPDirichletNitscheTEMP, "DirichletNitsche" );
            tIWGDirichletTEMP->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDirichletTEMP->set_property( tPropDirichletTEMP, "Dirichlet", mtk::Leader_Follower::LEADER );

            // create the IQIs
            // --------------------------------------------------------------------------------------
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUX->set_quantity_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
            tIQIUX->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } }, mtk::Leader_Follower::LEADER );
            tIQIUX->set_output_type_index( 0 );
            tIQIUX->set_name( "IQI_UX" );

            std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUY->set_quantity_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
            tIQIUY->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } }, mtk::Leader_Follower::LEADER );
            tIQIUY->set_output_type_index( 1 );
            tIQIUY->set_name( "IQI_UY" );

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            //----------------------------------------------------------------------------------------------------------
            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk1.set_IWGs( { tIWGBulkU_1, tIWGBulkTEMP_1 } );
            tSetBulk1.set_IQIs( { tIQIUX, tIQIUY, tIQITEMP } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
            tSetDirichlet.set_IWGs( { tIWGDirichletU, tIWGDirichletTEMP } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 2 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetDirichlet;

            // create model
            mdl::Model* tModel = new mdl::Model( tMeshManager,
                    0,
                    tSetInfo,
                    0,
                    false );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,    // STANDARD_WITH_OVERLAP
                    "./",
                    "MDL_Thermoelastic_Coupled_Test_2D_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_n_p1" },
                    { "UX", "UY", "TEMP" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_UX", "IQI_UY", "IQI_TEMP" } );
            tModel->set_output_manager( &tOutputData );

            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

            Vector< Vector< moris::Parameter_List > > tParameterlist( 8 );
            for ( uint Ik = 0; Ik < 8; Ik++ )
            {
                tParameterlist( Ik ).resize( 1 );
            }

            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );
            tParameterlist( 0 )( 0 ).set( "AZ_diagnostics", AZ_none );
            tParameterlist( 0 )( 0 ).set( "AZ_output", AZ_none );
            tParameterlist( 0 )( 0 ).set( "AZ_max_iter", 10000 );
            tParameterlist( 0 )( 0 ).set( "AZ_solver", AZ_gmres );
            tParameterlist( 0 )( 0 ).set( "AZ_subdomain_solve", AZ_ilu );
            tParameterlist( 0 )( 0 ).set( "AZ_graph_fill", 10 );
            tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );

            tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY;TEMP" );

            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY;TEMP" );

            tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

            tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML );
            tParameterlist( 7 )( 0 ).set( "ml_prec_type", "SA" );

            tSolverWarehouse.set_parameterlist( tParameterlist );

            tSolverWarehouse.initialize();

            tsa::Time_Solver* tTimeSolver = tSolverWarehouse.get_main_time_solver();

            tTimeSolver->set_output( 0, tSolverOutputCriteria_thermolast );
            tTimeSolver->solve();

            Matrix< DDRMat > tFullSolution;
            Matrix< DDRMat > tGoldSolution;
            tTimeSolver->get_full_solution( tFullSolution );

            std::string tMorisRoot    = moris::get_base_moris_dir();
            std::string tHdf5FilePath = tMorisRoot + "/projects/FEM/MDL/test/data/Thermoelastic_test_2d.hdf5";

            //        //------------------------------------------------------------------------------
            //        //    write solution ( uncomment this if you want to recreate solution files )
            //        //------------------------------------------------------------------------------
            //
            //        // create file
            //        hid_t tFileID = create_hdf5_file( tHdf5FilePath );
            //
            //        // error handler
            //        herr_t tStatus = 0;
            //
            //        // save data
            //        save_matrix_to_hdf5_file( tFileID, "Gold Solution", tFullSolution, tStatus );
            //
            //        // close file
            //        close_hdf5_file( tFileID );

            //------------------------------------------------------------------------------
            //    check solution
            //------------------------------------------------------------------------------

            // open file
            hid_t tFileID = open_hdf5_file( tHdf5FilePath );

            // error handler
            herr_t tStatus = 0;

            // read solution from file
            load_matrix_from_hdf5_file( tFileID, "Gold Solution", tGoldSolution, tStatus );

            // close file
            close_hdf5_file( tFileID );

            // verify solution
            bool tSolutionCheck = true;
            for ( uint i = 0; i < tFullSolution.numel(); i++ )
            {
                tSolutionCheck = tSolutionCheck && ( tFullSolution( i ) - tGoldSolution( i ) < 1e-03 );
                if ( !tSolutionCheck )
                {
                    std::cout << "tFullSolution( i ) " << tFullSolution( i ) << " tGoldSolution( i ) " << tGoldSolution( i ) << '\n';
                }
            }
            CHECK( tSolutionCheck );

            delete tModel;
            delete tInterpolationMesh;
        }
    }

    // TEST_CASE("2D XTK WITH HMR ThermoElastic 2D Input","[XTK_HMR_thermoelastic_2D_Input]")
    //{
    //     if(par_size()<=1)
    //     {
    //         uint tLagrangeMeshIndex = 0;
    //         std::string tFieldName = "Cylinder";
    //
    //         ParameterList tParameters = prm::create_hmr_parameter_list();
    //
    //         tParameters.set( "number_of_elements_per_dimension", std::string( "2, 1"));
    //         tParameters.set( "domain_dimensions", "2, 2" );
    //         tParameters.set( "domain_offset", "-1.0, -1.0" );
    //         tParameters.set( "lagrange_output_meshes",std::string( "0") );
    //
    //         tParameters.set( "lagrange_orders", std::string("1" ));
    //         tParameters.set( "lagrange_pattern", std::string("0" ));
    //         tParameters.set( "bspline_orders", std::string("1" ));
    //         tParameters.set( "bspline_pattern", std::string("0" ));
    //
    //         tParameters.set( "lagrange_to_bspline", "0" );
    //
    //         tParameters.set( "truncate_bsplines", 1 );
    //         tParameters.set( "refinement_buffer", 3 );
    //         tParameters.set( "staircase_buffer", 3 );
    //         tParameters.set( "initial_refinement", 0 );
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
    //         //// create field
    //         std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
    //         tField->evaluate_scalar_function( LvlSetPlane );
    //         //
    //         // for( uint k=0; k<2; ++k )
    //         // {
    //             // tHMR.flag_surface_elements_on_working_pattern( tField );
    //             // tHMR.perform_refinement_based_on_working_pattern( 0 );
    //
    //             // tField->evaluate_scalar_function( LvlSetCircle_2D );
    //         // }
    //
    //         tHMR.finalize();
    //
    //          tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );
    //
    //         std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);
    //
    //         moris::gen::Level_Set_Geometry_Field_HMR tPlaneFieldAsGeom(tField);
    //
    //         Vector<moris::gen::GEN_Geometry*> tGeometryVector = {&tPlaneFieldAsGeom};
    //
    //         size_t tModelDimension = 2;
    //         moris::gen::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    //         moris::gen::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tModelDimension);
    //
    //         xtk::Model tXTKModel(tModelDimension, tInterpolationMesh.get(), &tGeometryEngine);
    //
    //         tXTKModel.mVerbose = false;
    //
    //         //Specify decomposition Method and Cut Mesh ---------------------------------------
    //         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
    //         tXTKModel.decompose(tDecompositionMethods);
    //
    //         tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);
    //
    //         // get meshes
    //         xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
    //         xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
    //
    //         // place the pair in mesh manager
    //         mtk::Mesh_Manager tMeshManager;
    //         tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
    //
    //         std::string tInputFilePath = moris::get_base_moris_dir();
    //         tInputFilePath = tInputFilePath + "projects/FEM/MDL/test/data/Input_test.so";
    //
    //         // create a pointer to library for input reading
    //         std::shared_ptr< Library_IO > tLibrary = std::make_shared< Library_IO >( tInputFilePath );
    //
    //         // create model
    //         mdl::Model * tModel = new mdl::Model( tLibrary,
    //                                               &tMeshManager,
    //                                               0 );
    //
    //         // solve
    //         tModel->solve();
    //
    //         // clean up
    //         delete tModel;
    //     }
    // }

    TEST_CASE( "2D XTK WITH HMR ThermoElastic 2D Staggered", "[XTK_HMR_thermoelastic_2D_staggered]" )
    {
        if ( par_size() <= 1 )
        {
            uint        tLagrangeMeshIndex = 0;
            std::string tFieldName         = "Cylinder";

            Parameter_List tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "2, 1" );
            tParameters.set( "domain_dimensions", "2, 2" );
            tParameters.set( "domain_offset", "-1.0, -1.0" );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "1" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement();

            std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            //// create field
            std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

            tField->evaluate_scalar_function( LvlSetPlane );
            //
            // for( uint k=0; k<2; ++k )
            // {
            // tHMR.flag_surface_elements_on_working_pattern( tField );
            // tHMR.perform_refinement_based_on_working_pattern( 0 );

            // tField->evaluate_scalar_function( LvlSetCircle_2D );
            // }

            tHMR.finalize();

            tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            auto                                              tPlane          = std::make_shared< moris::gen::Line >( -500.0, 0.0, 1.0, 0.0 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

            size_t tModelDimension = 2;

            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometryVector;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            xtk::Model tXTKModel( tModelDimension, tInterpolationMesh, &tGeometryEngine );

            tXTKModel.mVerbose = false;

            // Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

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
            tPropConductivity1->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
            tPropConductivity2->set_parameters( { { { 5.0 } } } );
            tPropConductivity2->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichletU = std::make_shared< fem::Property >();
            tPropDirichletU->set_parameters( { { { 0.0 }, { 0.0 } } } );
            tPropDirichletU->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropDirichletTEMP = std::make_shared< fem::Property >();
            tPropDirichletTEMP->set_parameters( { { { 3.0 } } } );
            tPropDirichletTEMP->set_val_function( tConstValFunction );

            //        std::shared_ptr< fem::Property > tPropNeumannTEMP = std::make_shared< fem::Property >();
            //        tPropNeumann->set_parameters( { {{ 20.0 }} } );
            //        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

            std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
            tPropEMod1->set_parameters( { { { 1.0 } } } );
            tPropEMod1->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
            tPropEMod2->set_parameters( { { { 1.0 } } } );
            tPropEMod2->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tPropPoisson = std::make_shared< fem::Property >();
            tPropPoisson->set_parameters( { { { 0.0 } } } );
            tPropPoisson->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tCTE = std::make_shared< fem::Property >();
            tCTE->set_parameters( { { { 1.0 } } } );
            tCTE->set_val_function( tConstValFunction );

            std::shared_ptr< fem::Property > tTRef = std::make_shared< fem::Property >();
            tTRef->set_parameters( { { { 1.0 } } } );
            tTRef->set_val_function( tConstValFunction );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY }, { MSI::Dof_Type::TEMP } },
                    { "Displacement", "Temperature" } );
            tCMStrucLinIso1->set_property( tPropEMod1, "YoungsModulus" );
            tCMStrucLinIso1->set_property( tPropPoisson, "PoissonRatio" );
            tCMStrucLinIso1->set_property( tCTE, "CTE" );
            tCMStrucLinIso1->set_property( tTRef, "ReferenceTemperature" );
            tCMStrucLinIso1->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso1->set_space_dim( 2 );
            tCMStrucLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            tCMStrucLinIso2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY }, { MSI::Dof_Type::TEMP } },
                    { "Displacement", "Temperature" } );
            tCMStrucLinIso2->set_property( tPropEMod2, "YoungsModulus" );
            tCMStrucLinIso2->set_property( tPropPoisson, "PoissonRatio" );
            tCMStrucLinIso2->set_property( tCTE, "CTE" );
            tCMStrucLinIso2->set_property( tTRef, "ReferenceTemperature" );
            tCMStrucLinIso2->set_model_type( fem::Model_Type::PLANE_STRESS );
            tCMStrucLinIso2->set_space_dim( 2 );
            tCMStrucLinIso2->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
            tCMDiffLinIso1->set_space_dim( 2 );
            tCMDiffLinIso1->set_local_properties();

            std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
            tCMDiffLinIso2->set_space_dim( 2 );
            tCMDiffLinIso2->set_local_properties();

            //----------------------------------------------------------------------------------------------------------

            // define stabilization parameters
            fem::SP_Factory                                 tSPFactory;
            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheU = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitscheU->set_parameters( { { { 100.0 } } } );
            tSPDirichletNitscheU->set_property( tPropEMod1, "Material", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitscheTEMP = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
            tSPDirichletNitscheTEMP->set_parameters( { { { 1.0 } } } );
            tSPDirichletNitscheTEMP->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );

            //----------------------------------------------------------------------------------------------------------
            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWGBulkU_1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulkU_1->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_1->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_1->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulkU_2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
            tIWGBulkU_2->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_2->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGBulkU_2->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichletU = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichletU->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichletU->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            tIWGDirichletU->set_stabilization_parameter( tSPDirichletNitscheU, "DirichletNitsche" );
            tIWGDirichletU->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIWGDirichletU->set_property( tPropDirichletU, "Dirichlet", mtk::Leader_Follower::LEADER );

            //         std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            //         tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
            //         tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
            //         tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulkTEMP_1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkTEMP_1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
            //         tIWGBulkTEMP_1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGBulkTEMP_2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGBulkTEMP_2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGBulkTEMP_2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            //         tIWGBulkTEMP_2->set_property( tPropTempLoad2, "Load", mtk::Leader_Follower::LEADER );

            std::shared_ptr< fem::IWG > tIWGDirichletTEMP = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
            tIWGDirichletTEMP->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichletTEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
            tIWGDirichletTEMP->set_stabilization_parameter( tSPDirichletNitscheTEMP, "DirichletNitsche" );
            tIWGDirichletTEMP->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWGDirichletTEMP->set_property( tPropDirichletTEMP, "Dirichlet", mtk::Leader_Follower::LEADER );

            // create the IQIs
            // --------------------------------------------------------------------------------------
            fem::IQI_Factory tIQIFactory;

            std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUX->set_quantity_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
            tIQIUX->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } }, mtk::Leader_Follower::LEADER );
            tIQIUX->set_output_type_index( 0 );
            tIQIUX->set_name( "IQI_UX" );

            std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQIUY->set_quantity_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
            tIQIUY->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } }, mtk::Leader_Follower::LEADER );
            tIQIUY->set_output_type_index( 1 );
            tIQIUY->set_name( "IQI_UY" );

            std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
            tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
            tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
            tIQITEMP->set_output_type_index( 0 );
            tIQITEMP->set_name( "IQI_TEMP" );

            //----------------------------------------------------------------------------------------------------------

            fem::Set_User_Info tSetBulk1;
            tSetBulk1.set_mesh_set_name( "HMR_dummy_n_p1" );
            tSetBulk1.set_IWGs( { tIWGBulkU_1, tIWGBulkTEMP_1 } );
            tSetBulk1.set_IQIs( { tIQIUX, tIQIUY, tIQITEMP } );

            fem::Set_User_Info tSetDirichlet;
            tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p1" );
            tSetDirichlet.set_IWGs( { tIWGDirichletU, tIWGDirichletTEMP } );

            // create a cell of set info
            Vector< fem::Set_User_Info > tSetInfo( 2 );
            tSetInfo( 0 ) = tSetBulk1;
            tSetInfo( 1 ) = tSetDirichlet;

            // create model
            mdl::Model* tModel = new mdl::Model(
                    tMeshManager,
                    0,
                    tSetInfo,
                    0,
                    false );

            // define outputs
            // --------------------------------------------------------------------------------------
            vis::Output_Manager tOutputData;
            tOutputData.set_outputs( 0,
                    vis::VIS_Mesh_Type::STANDARD,    // STANDARD_WITH_OVERLAP
                    "./",
                    "MDL_Thermoelastic_Staggered_Test_2D_Output.exo",
                    "./",
                    "temp.exo",
                    { "HMR_dummy_n_p1" },
                    { "UX", "UY", "TEMP" },
                    { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                    { "IQI_UX", "IQI_UY", "IQI_TEMP" } );
            tModel->set_output_manager( &tOutputData );

            sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

            Vector< Vector< moris::Parameter_List > > tParameterlist( 8 );

            tParameterlist( 0 ).resize( 3 );
            tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );
            tParameterlist( 0 )( 0 ).set( "AZ_diagnostics", AZ_none );
            tParameterlist( 0 )( 0 ).set( "AZ_output", AZ_none );
            tParameterlist( 0 )( 0 ).set( "AZ_max_iter", 10000 );
            tParameterlist( 0 )( 0 ).set( "AZ_solver", AZ_gmres );
            tParameterlist( 0 )( 0 ).set( "AZ_subdomain_solve", AZ_ilu );
            tParameterlist( 0 )( 0 ).set( "AZ_graph_fill", 10 );
            tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );

            tParameterlist( 0 )( 1 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AZTEC_IMPL );
            tParameterlist( 0 )( 1 ).set( "preconditioners", "1" );

            tParameterlist( 0 )( 2 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );

            tParameterlist( 1 ).resize( 3 );
            tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );
            tParameterlist( 1 )( 1 ) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 1 )( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
            tParameterlist( 1 )( 2 ) = moris::prm::create_linear_solver_parameter_list();
            tParameterlist( 1 )( 2 ).set( "DLA_Linear_solver_algorithms", "2" );

            tParameterlist( 2 ).resize( 3 );
            tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
            tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );
            tParameterlist( 2 )( 0 ).set( "NLA_Linear_solver", 0 );
            tParameterlist( 2 )( 1 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
            tParameterlist( 2 )( 1 ).set( "NLA_combined_res_jac_assembly", false );
            tParameterlist( 2 )( 1 ).set( "NLA_Linear_solver", 1 );
            tParameterlist( 2 )( 2 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
            tParameterlist( 2 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
            tParameterlist( 2 )( 2 ).set( "NLA_combined_res_jac_assembly", false );
            tParameterlist( 2 )( 2 ).set( "NLA_Linear_solver", 2 );

            tParameterlist( 3 ).resize( 3 );
            tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 0 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
            tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );
            tParameterlist( 3 )( 0 ).set( "NLA_Nonlinear_solver_algorithms", "0" );
            tParameterlist( 3 )( 1 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 1 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NEWTON_SOLVER ) );
            tParameterlist( 3 )( 1 ).set( "NLA_DofTypes", "TEMP" );
            tParameterlist( 3 )( 1 ).set( "NLA_Nonlinear_solver_algorithms", "1" );
            tParameterlist( 3 )( 2 ) = moris::prm::create_nonlinear_solver_parameter_list();
            tParameterlist( 3 )( 2 ).set( "NLA_Solver_Implementation", static_cast< uint >( moris::NLA::NonlinearSolverType::NLBGS_SOLVER ) );
            tParameterlist( 3 )( 2 ).set( "NLA_Sub_Nonlinear_Solver", "1,0" );
            tParameterlist( 3 )( 2 ).set( "NLA_DofTypes", "UX,UY;TEMP" );
            tParameterlist( 3 )( 2 ).set( "NLA_Nonlinear_solver_algorithms", "2" );

            tParameterlist( 4 ).resize( 1 );
            tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
            tParameterlist( 4 )( 0 ).set( "TSA_Nonlinear_Solver", 2 );

            tParameterlist( 5 ).resize( 1 );
            tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
            tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY;TEMP" );

            tParameterlist( 6 ).resize( 1 );
            tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();

            tParameterlist( 7 ).resize( 2 );
            tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML );
            tParameterlist( 7 )( 0 ).set( "ml_prec_type", "SA" );

            tParameterlist( 7 )( 1 ) = moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK );
            tParameterlist( 7 )( 1 ).set( "ifpack_prec_type", "ILU" );

            tSolverWarehouse.set_parameterlist( tParameterlist );
            tSolverWarehouse.initialize();

            tsa::Time_Solver* tTimeSolver = tSolverWarehouse.get_main_time_solver();
            tTimeSolver->set_output( 0, tSolverOutputCriteria_thermolast );
            tTimeSolver->solve();

            //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //        // STEP 1: create linear solver and algorithm
            //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //        Vector< enum MSI::Dof_Type > tDofTypesT( 1 );            tDofTypesT( 0 ) = MSI::Dof_Type::TEMP;
            //        Vector< enum MSI::Dof_Type > tDofTypesU( 2 );            tDofTypesU( 0 ) = MSI::Dof_Type::UX;              tDofTypesU( 1 ) = MSI::Dof_Type::UY;
            //
            //        dla::Solver_Factory  tSolFactory;
            //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
            //
            //        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
            //        tLinearSolverParameterList.set( "AZ_output", AZ_none );
            //        tLinearSolverParameterList.set( "AZ_max_iter", 10000 );
            //        tLinearSolverParameterList.set( "AZ_solver", AZ_gmres );
            //        tLinearSolverParameterList.set( "AZ_subdomain_solve", AZ_ilu );
            //        tLinearSolverParameterList.set( "AZ_graph_fill", 10 );
            ////        tLinearSolverParameterList.set( "ml_prec_type", "SA" );
            //
            //        dla::Linear_Solver tLinSolver;
            //        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
            //
            //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //        // STEP 2: create nonlinear solver and algorithm
            //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //        NLA::Nonlinear_Solver_Factory tNonlinFactory;
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMain = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NLBGS_SOLVER );
            //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythic = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            ////        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
            //
            //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_iter")   = 3;
            //        tNonlinearSolverAlgorithmMain->set_param("NLA_max_iter")   = 1;
            //        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
            //        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
            //        //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;
            //
            //        tNonlinearSolverAlgorithmMonolythic->set_linear_solver( &tLinSolver );
            ////        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );
            //
            //        NLA::Nonlinear_Solver tNonlinearSolverMain;
            //        NLA::Nonlinear_Solver tNonlinearSolverMonolythicT;
            //        NLA::Nonlinear_Solver tNonlinearSolverMonolythicU;
            //        tNonlinearSolverMain       .set_nonlinear_algorithm( tNonlinearSolverAlgorithmMain, 0 );
            //        tNonlinearSolverMonolythicT.set_nonlinear_algorithm( tNonlinearSolverAlgorithmMonolythic, 0 );
            //        tNonlinearSolverMonolythicU.set_nonlinear_algorithm( tNonlinearSolverAlgorithmMonolythic, 0 );
            //
            //        tNonlinearSolverMain       .set_dof_type_list( tDofTypesU );
            //        tNonlinearSolverMain       .set_dof_type_list( tDofTypesT );
            //        tNonlinearSolverMonolythicT.set_dof_type_list( tDofTypesT );
            //        tNonlinearSolverMonolythicU.set_dof_type_list( tDofTypesU );
            //
            //        tNonlinearSolverMonolythicU.set_secondary_dof_type_list(tDofTypesT);
            //
            //        tNonlinearSolverMain.set_sub_nonlinear_solver( &tNonlinearSolverMonolythicT );
            //        tNonlinearSolverMain.set_sub_nonlinear_solver( &tNonlinearSolverMonolythicU );
            //
            //        // Create solver database
            //        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
            //
            //        tNonlinearSolverMain       .set_solver_warehouse( &tSolverWarehouse );
            //        tNonlinearSolverMonolythicT.set_solver_warehouse( &tSolverWarehouse );
            //        tNonlinearSolverMonolythicU.set_solver_warehouse( &tSolverWarehouse );
            //
            //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //        // STEP 3: create time Solver and algorithm
            //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            //        tsa::Time_Solver_Factory tTimeSolverFactory;
            //        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
            //
            //        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );
            //
            //        tsa::Time_Solver tTimeSolver;
            //        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
            //        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
            //
            //        tTimeSolver.set_dof_type_list( tDofTypesU );
            //        tTimeSolver.set_dof_type_list( tDofTypesT );
            //
            //        //------------------------------------------------------------------------------
            //        tTimeSolver.solve();

            Matrix< DDRMat > tFullSolution;
            Matrix< DDRMat > tGoldSolution;
            tTimeSolver->get_full_solution( tFullSolution );

            std::string tMorisRoot    = moris::get_base_moris_dir();
            std::string tHdf5FilePath = tMorisRoot + "/projects/FEM/MDL/test/data/Thermoelastic_test_2d.hdf5";

            //------------------------------------------------------------------------------
            //    check solution
            //------------------------------------------------------------------------------

            // open file
            hid_t tFileID = open_hdf5_file( tHdf5FilePath );

            // error handler
            herr_t tStatus = 0;

            // read solution from file
            load_matrix_from_hdf5_file( tFileID, "Gold Solution", tGoldSolution, tStatus );

            // close file
            close_hdf5_file( tFileID );

            // verify solution

            bool tSolutionCheck = true;
            for ( uint i = 0; i < tFullSolution.numel(); i++ )
            {
                tSolutionCheck = tSolutionCheck && ( tFullSolution( i ) - tGoldSolution( i ) < 1e-03 );
                if ( !tSolutionCheck )
                {
                    std::cout << "tFullSolution( i ) " << tFullSolution( i ) << " tGoldSolution( i ) " << tGoldSolution( i ) << '\n';
                }
            }

            CHECK( tSolutionCheck );

            delete tModel;
            delete tInterpolationMesh;
        }
    }

    // TEST_CASE("2D XTK WITH HMR Struc Interface 3D","[XTK_HMR_Struc_Interface_3D]")
    //{
    //     if(par_size()<=1)
    //     {
    //         uint tLagrangeMeshIndex = 0;
    //         std::string tFieldName = "Cylinder";
    //
    //         uint tSpatialDimension = 3;
    //
    //          hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
    //
    //          tParameters.set( "number_of_elements_per_dimension", "22, 8, 2");
    //          tParameters.set( "domain_dimensions", "6, 2, 1" );
    //          tParameters.set( "domain_offset", "-3.0, -1.0, -0.5" );
    //          tParameters.set( "lagrange_output_meshes", "0" );
    //
    //          tParameters.set( "lagrange_orders", "1" );
    //          tParameters.set( "lagrange_pattern", "0" );
    //          tParameters.set( "bspline_orders", "1" );
    //          tParameters.set( "bspline_pattern", "0" );
    //
    //          tParameters.set( "lagrange_to_bspline", "0" );
    //
    //          tParameters.set( "truncate_bsplines", 1 );
    //          tParameters.set( "refinement_buffer", 3 );
    //          tParameters.set( "staircase_buffer", 3 );
    //          tParameters.set( "initial_refinement", 0 );
    //
    //          tParameters.set( "use_multigrid", 0 );
    //          tParameters.set( "severity_level", 2 );
    //
    //          hmr::HMR tHMR( tParameters );
    //
    //         //initial refinement
    //          tHMR.perform_initial_refinement( 0 );
    //
    //          std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
    //
    //        //  create field
    //        Vector<std::shared_ptr< moris::hmr::Field > > tHMRFields;
    //        tHMRFields.resize(2);
    //
    //        // create field
    //        tHMRFields(0) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    //        tHMRFields(1) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    //
    //        tHMRFields(0)->evaluate_scalar_function( LevelSetFunction_star1 );
    //        tHMRFields(1)->evaluate_scalar_function( Plane4MatMDL1 );
    //
    //        for( uint k=0; k<2; ++k )
    //        {
    //            tHMR.flag_surface_elements_on_working_pattern( tHMRFields(0) );
    //            tHMR.flag_surface_elements_on_working_pattern( tHMRFields(1) );
    //
    //            tHMR.perform_refinement_based_on_working_pattern( 0 );
    //
    //            tHMRFields(0)->evaluate_scalar_function( LevelSetFunction_star1 );
    //            tHMRFields(1)->evaluate_scalar_function( Plane4MatMDL1 );
    //        }
    //
    //       tHMR.finalize();
    //
    ////      tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_3d.e" );
    //
    //      std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
    //
    //      xtk::Geom_Field tCircleFieldAsGeom(tHMRFields(0));
    //      xtk::Geom_Field tPlaneFieldAsGeom2(tHMRFields(1));
    //      Vector<xtk::Geometry*> tGeometryVector = {&tCircleFieldAsGeom,&tPlaneFieldAsGeom2};
    //
    //      xtk::Phase_Table     tPhaseTable (tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2);
    //      xtk::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tSpatialDimension);
    //      xtk::Model           tXTKModel(tSpatialDimension,tInterpMesh.get(),&tGeometryEngine);
    //      tXTKModel.mVerbose = false;
    //
    //      Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    //      tXTKModel.decompose(tDecompositionMethods);
    //      tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);
    //
    //      // get meshes
    //      xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
    //      xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
    //
    //      // place the pair in mesh manager
    //      mtk::Mesh_Manager tMeshManager;
    //      tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
    //
    //      Vector< fem::IWG_User_Defined_Info > tBulkIWG(1);
    //      Vector< fem::IWG_User_Defined_Info > tDBCIWG(1);
    //      Vector< fem::IWG_User_Defined_Info > tNBCIWG(1);
    //      Vector< fem::IWG_User_Defined_Info > tIntIWG(1);
    //
    //        // create IWG user defined info
    //      tBulkIWG( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_BULK,
    //                                                  { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
    //                                                  {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
    //                                                  Vector< fem::Property_Type >( 0 ),
    //                                                  { fem::Constitutive_Type::STRUC_LIN_ISO } );
    //
    //      tDBCIWG( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE,
    //                                                 { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
    //                                                 {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
    //                                                 { fem::Property_Type::STRUC_DIRICHLET },
    //                                                 { fem::Constitutive_Type::STRUC_LIN_ISO } );
    //
    //      tNBCIWG( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_NEUMANN,
    //                                                 { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
    //                                                 {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },},
    //                                                 { fem::Property_Type::STRUC_NEUMANN },
    //                                                 Vector< fem::Constitutive_Type >( 0 ) );
    //      tIntIWG( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE,
    //                                                 { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ },
    //                                                 {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
    //                                                 Vector< fem::Property_Type >( 0 ),
    //                                                 {fem::Constitutive_Type::STRUC_LIN_ISO },
    //                                                 {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
    //                                                 Vector< fem::Property_Type >( 0 ),
    //                                                 {fem::Constitutive_Type::STRUC_LIN_ISO } );
    //
    //        Vector< Vector< fem::IWG_User_Defined_Info > > tIWGUserDefinedInfo( 14 );
    //
    //        tIWGUserDefinedInfo( 0 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 1 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 2 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 3 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 4 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 5 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 6 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 7 )  = tBulkIWG;
    //        tIWGUserDefinedInfo( 8 )  = tDBCIWG;
    //        tIWGUserDefinedInfo( 9 )  = tNBCIWG;
    //        tIWGUserDefinedInfo( 10 ) = tIntIWG;
    //        tIWGUserDefinedInfo( 11 ) = tIntIWG;
    //        tIWGUserDefinedInfo( 12 ) = tIntIWG;
    //        tIWGUserDefinedInfo( 13 ) = tIntIWG;
    //
    //        // create property user defined info
    //        fem::Property_User_Defined_Info tYoungs_Modulus( fem::Property_Type::YOUNGS_MODULUS,
    //                Vector< Vector< MSI::Dof_Type > >( 0 ),
    //                {{{ 1.0 }}},
    //                tConstValFunction,
    //                Vector< fem::PropertyFunc >( 0 ) );
    //        fem::Property_User_Defined_Info tYoungs_Modulus2( fem::Property_Type::YOUNGS_MODULUS,
    //                Vector< Vector< MSI::Dof_Type > >( 0 ),
    //                {{{ 1.0 }}},
    //                tConstValFunction,
    //                Vector< fem::PropertyFunc >( 0 ) );
    //        fem::Property_User_Defined_Info tPoissons_Ratio( fem::Property_Type::POISSONS_RATIO,
    //                Vector< Vector< MSI::Dof_Type > >( 0 ),
    //                {{{ 0.0 }}},
    //                tConstValFunction,
    //                Vector< fem::PropertyFunc >( 0 ) );
    //        fem::Property_User_Defined_Info tStrucDirichlet( fem::Property_Type::STRUC_DIRICHLET,
    //                Vector< Vector< MSI::Dof_Type > >( 0 ),
    //                {{{ 0.0, 0.0, 0.0 }}},
    //                tConstValFunction,
    //                Vector< fem::PropertyFunc >( 0 ) );
    //        fem::Property_User_Defined_Info tStrucNeumann( fem::Property_Type::STRUC_NEUMANN,
    //                Vector< Vector< MSI::Dof_Type > >( 0 ),
    //                {{{ 1.0, 0.0, 0.0 }}},
    //                tConstValFunction,
    //                Vector< fem::PropertyFunc >( 0 ) );
    //
    //        // create property user defined info
    //        Vector< Vector< Vector< fem::Property_User_Defined_Info > > > tPropertyUserDefinedInfo( 14 );
    //        tPropertyUserDefinedInfo( 0 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 0 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 0 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 0 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 1 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 1 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 1 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 1 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 2 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 2 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 2 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 2 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 3 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 3 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 3 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 3 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 4 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 4 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 4 )( 0 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 4 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 5 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 5 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 5 )( 0 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 5 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 6 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 6 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 6 )( 0 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 6 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 7 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 7 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 7 )( 0 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 7 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 8 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 8 )( 0 ).resize( 3 );
    //        tPropertyUserDefinedInfo( 8 )( 0 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 8 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 8 )( 0 )( 2 ) = tStrucDirichlet;
    //        tPropertyUserDefinedInfo( 9 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 9 )( 0 ).resize( 1 );
    //        tPropertyUserDefinedInfo( 9 )( 0 )( 0 ) = tStrucNeumann;
    //        tPropertyUserDefinedInfo( 10 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 10 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 10 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 10 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 10 )( 1 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 10 )( 1 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 10 )( 1 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 11 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 11 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 11 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 11 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 11 )( 1 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 11 )( 1 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 11 )( 1 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 12 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 12 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 12 )( 0 )( 0 ) = tYoungs_Modulus2;
    //        tPropertyUserDefinedInfo( 12 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 12 )( 1 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 12 )( 1 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 12 )( 1 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 13 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 13 )( 0 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 13 )( 0 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 13 )( 0 )( 1 ) = tPoissons_Ratio;
    //        tPropertyUserDefinedInfo( 13 )( 1 ).resize( 2 );
    //        tPropertyUserDefinedInfo( 13 )( 1 )( 0 ) = tYoungs_Modulus;
    //        tPropertyUserDefinedInfo( 13 )( 1 )( 1 ) = tPoissons_Ratio;
    //
    //        fem::Constitutive_User_Defined_Info tStrucLinIso( fem::Constitutive_Type::STRUC_LIN_ISO,
    //                                                          {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},
    //                                                          { fem::Property_Type::YOUNGS_MODULUS, fem::Property_Type::POISSONS_RATIO } );
    //
    //        // create constitutive user defined info
    //        Vector< Vector< Vector< fem::Constitutive_User_Defined_Info > > > tConstitutiveUserDefinedInfo( 14 );
    //        tConstitutiveUserDefinedInfo( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 0 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 0 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 1 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 1 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 1 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 2 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 2 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 2 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 3 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 3 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 3 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 4 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 4 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 4 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 5 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 5 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 5 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 6 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 6 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 6 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 7 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 7 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 7 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 8 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 8 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 8 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 9 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 10 ).resize( 2 );
    //        tConstitutiveUserDefinedInfo( 10 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 10 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 10 )( 1 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 10 )( 1 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 11 ).resize( 2 );
    //        tConstitutiveUserDefinedInfo( 11 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 11 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 11 )( 1 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 11 )( 1 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 12 ).resize( 2 );
    //        tConstitutiveUserDefinedInfo( 12 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 12 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 12 )( 1 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 12 )( 1 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 13 ).resize( 2 );
    //        tConstitutiveUserDefinedInfo( 13 )( 0 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 13 )( 0 )( 0 ) = tStrucLinIso;
    //        tConstitutiveUserDefinedInfo( 13 )( 1 ).resize( 1 );
    //        tConstitutiveUserDefinedInfo( 13 )( 1 )( 0 ) = tStrucLinIso;
    //
    //        // create a list of active block-sets
    //        // create a list of active block-sets
    //        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
    //        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
    //        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
    //        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);
    //
    //        Vector< moris_index >  tSetList = { tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p0"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p0"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p1"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p1"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p2"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p2"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p3"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p3"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("SideSet_4_n_p2"),
    //                                                 tEnrIntegMesh.get_set_index_by_name("SideSet_2_n_p3"),
    //                                                 tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName01),
    //                                                 tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName02),
    //                                                 tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName13),
    //                                                 tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName23)};
    //
    //        Vector< fem::Element_Type > tSetTypeList = { fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::BULK,
    //                                                          fem::Element_Type::SIDESET,
    //                                                          fem::Element_Type::SIDESET,
    //                                                          fem::Element_Type::DOUBLE_SIDESET,
    //                                                          fem::Element_Type::DOUBLE_SIDESET,
    //                                                          fem::Element_Type::DOUBLE_SIDESET,
    //                                                          fem::Element_Type::DOUBLE_SIDESET };
    //
    //        uint tBSplineMeshIndex = 0;
    //        // create model
    //        mdl::Model * tModel = new mdl::Model( &tMeshManager,
    //                tBSplineMeshIndex,
    //                tSetList, tSetTypeList,
    //                tIWGUserDefinedInfo,
    //                tPropertyUserDefinedInfo,
    //                tConstitutiveUserDefinedInfo );
    //
    //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //        // STEP 1: create linear solver and algorithm
    //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //
    //        dla::Solver_Factory  tSolFactory;
    //        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AZTEC_IMPL );
    //        //            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::PETSC );
    //
    //        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
    //        tLinearSolverParameterList.set( "AZ_output", AZ_none );
    //        tLinearSolverParameterList.set( "AZ_max_iter", 5000 );
    //        tLinearSolverParameterList.set( "AZ_solver", AZ_gmres_condnum );
    //        tLinearSolverParameterList.set( "AZ_subdomain_solve", AZ_ilu );
    //        tLinearSolverParameterList.set( "AZ_ilut_fill", 10.0 );
    //        tLinearSolverParameterList.set( "rel_residual", 1e-8 );
    //
    ////            tLinearSolverParameterList.set( "ml_prec_type", "SA" );
    //
    //        dla::Linear_Solver tLinSolver;
    //
    //        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
    //
    //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //        // STEP 2: create nonlinear solver and algorithm
    //        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Epetra );
    //        //    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );
    //
    //        NLA::Nonlinear_Solver_Factory tNonlinFactory;
    //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
    //
    //        tNonlinearSolverParameterList.set( "NLA_max_iter", 2 );
    //        //        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
    //        //        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
    //        //        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
    //
    //        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
    //
    //        NLA::Nonlinear_Solver tNonlinearSolver;
    //
    //        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
    //
    //        tNonlinearSolver.solve( tNonlinearProblem );
    //
    //        std::cout<<" Solution Vector "<<std::endl;
    //        tNonlinearProblem->get_full_vector()->print();
    //
    //
    //        delete tModel;
    //    }
    //}
}    // namespace moris
