/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_Sensitivity_Test.cpp
 *
 */

#include "catch.hpp"

#include "cl_MDL_Model.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "typedefs.hpp"

#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"              //FEM/INT/src

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_TSA_Time_Solver.hpp"
#include "cl_SOL_Warehouse.hpp"

#include "fn_norm.hpp"

#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Property.hpp"

#include "cl_VIS_Output_Manager.hpp"

#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

real plane_evaluate_field_value(const moris::Matrix< DDRMat >    & aCoordinates,
                                 const moris::Cell< moris::real > & aParameters)
{
    moris::real mXC = 0.11;
    moris::real mYC = 0.11;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return (mNx*(aCoordinates(0)-mXC) + mNy*(aCoordinates(1)-mYC));
}

void evaluate_sensitivity(const moris::Matrix< DDRMat >    & aCoordinates,
                                const moris::Cell< moris::real > & aParameters,
                                Matrix<DDRMat>& aSensitivities)
{
    MORIS_ERROR( false, "sensitivities not implemented");
}

void tConstValFunction2MatMDL
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValDvFunction_FDTest
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( ge::PDV_Type::DENSITY )->val();
}

void tFIDerDvFunction_FDTest
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( ge::PDV_Type::DENSITY )->N();
}

real density_function_1(const Matrix<DDRMat>& aCoordinates, const Cell<real>& aPropertyVariables)
{
    return 1.0;
}

real density_function_2(const Matrix<DDRMat>& aCoordinates, const Cell<real>& aPropertyVariables)
{
    return 2.0;
}

TEST_CASE("Sensitivity test","[Sensitivity test]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Geometry";

        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        ParameterList tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "5,2" );
        tParameters.set( "domain_dimensions", "6,2" );
        tParameters.set( "domain_offset", "-3.0,-1.0" );
        tParameters.set( "domain_sidesets", "1,2,3,4" );
        tParameters.set( "lagrange_output_meshes", "0" );

        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_to_bspline", "0" );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );
        tParameters.set( "initial_refinement", "0" );
        tParameters.set( "initial_refinement_pattern", "0" );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement();

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//        for( uint k=0; k<tNumRef; ++k )
//        {
//            Cell<std::shared_ptr<moris::ge::Level_Set_Geometry>> tGeomVec(2);
//            tGeomVec(0) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tROuter);
//            tGeomVec(1) = std::make_shared<moris::ge::Circle>(tCenterPoint(0), tCenterPoint(1), tRInner);
//
//            moris::ge::Phase_Table     tPhaseTable( tGeomVec.size() );
//            moris::ge::Geometry_Engine tGENGeometryEngine( tGeomVec, tPhaseTable, 2 );
//
//            moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );
//
//            uint tNumIPNodes = tMesh->get_num_nodes();
//            Matrix<DDRMat> tFieldData( tNumIPNodes,1 );
//            Matrix<DDRMat> tFieldData0( tNumIPNodes,1 );
//
//            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumIPNodes );
//            Matrix< DDRMat > tCoords( tNumIPNodes, 2 );
//            for( uint i = 0; i < tNumIPNodes; i++ )
//            {
//                tCoords.set_row( i, tMesh->get_mtk_vertex(i).get_coords() );
//            }
//
//            tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );
//
//            for(uint i=0; i<tNumIPNodes; i++)
//            {
//                tFieldData( i )  = tGENGeometryEngine.get_entity_phase_val( i, 0 );
//                tFieldData0( i ) = tGENGeometryEngine.get_entity_phase_val( i, 1 );
//            }
//
//            tHMR.based_on_field_put_elements_on_queue( tFieldData, tLagrangeMeshIndex );
//            tHMR.based_on_field_put_elements_on_queue( tFieldData0, tLagrangeMeshIndex );
//
//            tHMR.perform_refinement_based_on_working_pattern( 0, false );
//        }
        tHMR.finalize();

//        tHMR.save_to_exodus( 0, "./xtk_exo/hmr_ip_mesh_for_sensitivity_test.e" );

        hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::Cell< std::shared_ptr<moris::ge::Level_Set_Geometry> > tGeometryVector(1);

        Matrix<DDRMat> tADVs( 0, 0 );
        auto tField = std::make_shared<moris::ge::User_Defined_Field>(tADVs,
                                                                     Matrix<DDUMat> (0, 0),
                                                                     Matrix<DDUMat> (0, 0),
                                                                     Matrix<DDRMat> (0, 0),
                                                                            &plane_evaluate_field_value,
                                                                            &evaluate_sensitivity);
        tGeometryVector( 0 ) = std::make_shared< Level_Set_Field >( tField );

        size_t tModelDimension = 2;
        moris::ge::Geometry_Engine tGeometryEngine(tGeometryVector, tModelDimension);
        xtk::Model tXTKModel(tModelDimension,tInterpMesh,&tGeometryEngine);
        tXTKModel.mVerbose = false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);
        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        // Write mesh
        moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
        //writer.write_mesh("","./mdl_exo/xtk_ig_mesh_for_sensitivity_test.e", "", "temp.exo");

        // Write the fields
//        writer.set_time(0.0);
//        writer.close_file();

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property > ();
        tPropConductivity1->set_parameters( { {{ 1.0 }} } );
        tPropConductivity1->set_dv_type_list( {{ ge::PDV_Type::DENSITY }} );
        tPropConductivity1->set_val_function( tFIValDvFunction_FDTest );
        tPropConductivity1->set_dv_derivative_functions( { tFIDerDvFunction_FDTest } );

        std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property > ();
        tPropConductivity2->set_parameters( { {{ 1.0 }} } );
        tPropConductivity2->set_dv_type_list( {{ ge::PDV_Type::DENSITY }} );
        tPropConductivity2->set_val_function( tFIValDvFunction_FDTest );
        tPropConductivity2->set_dv_derivative_functions( { tFIDerDvFunction_FDTest } );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { {{ 100.0 }} } );
        tPropTempLoad1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
        tPropTempLoad2->set_parameters( { {{ 100.0 }} } );
        tPropTempLoad2->set_val_function( tConstValFunction2MatMDL );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
        tCMDiffLinIso1->set_space_dim( 2 );
        tCMDiffLinIso1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
        tCMDiffLinIso2->set_space_dim( 2 );
        tCMDiffLinIso2->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
        tSPNitscheInterface->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::FOLLOWER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk1 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk1->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk1->set_property( tPropTempLoad1, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGBulk2 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGBulk2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGBulk2->set_property( tPropTempLoad2, "Load", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGInterface =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Leader_Follower::FOLLOWER );
        tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
        tIWGInterface->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        Cell<std::string> tRequestedIQINames(2);
        tRequestedIQINames(0) = "IQI_StrainEnergy1";
        tRequestedIQINames(1) = "IQI_StrainEnergy2";

        std::shared_ptr< fem::IQI > tIQI_SE1 = tIQIFactory.create_IQI( fem::IQI_Type::STRAIN_ENERGY );
        tIQI_SE1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
        tIQI_SE1->set_constitutive_model( tCMDiffLinIso1, "Elast", mtk::Leader_Follower::LEADER );
        tIQI_SE1->set_name(tRequestedIQINames(0));

        std::shared_ptr< fem::IQI > tIQI_SE2 = tIQIFactory.create_IQI( fem::IQI_Type::STRAIN_ENERGY );
        tIQI_SE2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
        tIQI_SE2->set_constitutive_model( tCMDiffLinIso2, "Elast", mtk::Leader_Follower::LEADER );
        tIQI_SE2->set_name(tRequestedIQINames(1));

        std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIQITEMP->set_output_type_index( 0 );
        tIQITEMP->set_name("IQI_TEMP");

        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name(0,0,1);
        std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
        tSetBulk1.set_IWGs( { tIWGBulk1 } );
        tSetBulk1.set_IQIs( { tIQITEMP, tIQI_SE1 } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
        tSetBulk2.set_IWGs( { tIWGBulk1 } );
        tSetBulk2.set_IQIs( { tIQITEMP, tIQI_SE1 } );

        fem::Set_User_Info tSetBulk3;
        tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
        tSetBulk3.set_IWGs( { tIWGBulk2 } );
        tSetBulk3.set_IQIs( { tIQITEMP, tIQI_SE2 } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk4.set_IWGs( { tIWGBulk2 } );
        tSetBulk4.set_IQIs( { tIQITEMP, tIQI_SE2 } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_2_n_p1" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_set_name( "SideSet_4_n_p0" );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetInterface;
        tSetInterface.set_mesh_set_name( tDblInterfaceSideSetName );
        tSetInterface.set_IWGs( { tIWGInterface } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 7 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetBulk3;
        tSetInfo( 3 ) = tSetBulk4;
        tSetInfo( 4 ) = tSetDirichlet;
        tSetInfo( 5 ) = tSetNeumann;
        tSetInfo( 6 ) = tSetInterface;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               0,
                                               tSetInfo,
                                               0, false );

        // output manager
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                                 vis::VIS_Mesh_Type::STANDARD, //STANDARD_WITH_OVERLAP
                                 "./",
                                 "MDL_Sensitivity_Test.exo",
                                 "./",
                                  "temp.exo",
                                 { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_n_p0", "HMR_dummy_n_p1" },
                                 {  "TEMP" },
                                 {  vis::Field_Type::NODAL },
                                 {  vis::Output_Type::TEMP } );
        tModel->set_output_manager( &tOutputData );

        // PDVs and design variable interface
        tGeometryEngine.register_mesh(&tMeshManager);

        Cell<Cell<Cell<ge::PDV_Type>>> tPdvTypes(15);
        Cell<Cell<Cell<ge::PDV_Type>>> tIGPdvTypes(100);
        for (uint tBulkSetIndex = 0; tBulkSetIndex < 4; tBulkSetIndex++)
        {
            tPdvTypes(tBulkSetIndex).resize(1);
            tPdvTypes(tBulkSetIndex)(0).resize(1);
            tPdvTypes(tBulkSetIndex)(0)(0) = ge::PDV_Type::DENSITY;
        }
        tGeometryEngine.create_ip_pdv_hosts(tPdvTypes);

        Cell< ge::PDV_Type> tRequestedType( 1, ge::PDV_Type::DENSITY );
        reinterpret_cast< ge::Pdv_Host_Manager* >(tGeometryEngine.get_design_variable_interface())->
                create_ig_pdv_hosts(0,Cell< Matrix< DDSMat >>(100), tIGPdvTypes);
        moris::Cell< ge::PDV_Type > tMatPdvTypes = { ge::PDV_Type::DENSITY };
        reinterpret_cast< ge::Pdv_Host_Manager* >(tGeometryEngine.get_design_variable_interface())->
                set_ip_requested_dv_types( tMatPdvTypes );

        Matrix<DDRMat> tAdvs(1, 1, 0.0);
        auto tDensityField1 = std::make_shared<ge::User_Defined_Field>(
                tAdvs,
                Matrix<DDUMat>(0, 0),
                Matrix<DDUMat>(0, 0),
                Matrix<DDRMat>(0, 0),
                Cell<std::shared_ptr<ge::Property>>(0),
                &density_function_1,
                &evaluate_sensitivity);
        auto tDensityProperty1 = std::make_shared< ge::Property >( tDensityField1 );
        auto tDensityField2 = std::make_shared<ge::User_Defined_Field>(
                tAdvs,
                Matrix<DDUMat>(0, 0),
                Matrix<DDUMat>(0, 0),
                Matrix<DDRMat>(0, 0),
                Cell<std::shared_ptr<ge::Property>>(0),
                &density_function_2,
                &evaluate_sensitivity);
        auto tDensityProperty2 = std::make_shared< ge::Property >( tDensityField2 );
        tGeometryEngine.assign_ip_hosts_by_set_index(0, tDensityProperty1, ge::PDV_Type::DENSITY);
        tGeometryEngine.assign_ip_hosts_by_set_index(1, tDensityProperty1, ge::PDV_Type::DENSITY);
        tGeometryEngine.assign_ip_hosts_by_set_index(2, tDensityProperty2, ge::PDV_Type::DENSITY);
        tGeometryEngine.assign_ip_hosts_by_set_index(3, tDensityProperty2, ge::PDV_Type::DENSITY);

        tModel->set_design_variable_interface(tGeometryEngine.get_design_variable_interface());

//        Cell<Matrix<DDRMat>> tDvValues;
//        Matrix<IndexMat> tNodeIndices(70, 1);
//        for (uint tNodeIndex = 0; tNodeIndex < 70; tNodeIndex++)
//        {
//            tNodeIndices(tNodeIndex) = tNodeIndex;
//        }
//        tGeometryEngine.get_design_variable_interface()->get_ip_pdv_value(tNodeIndices, tPdvTypes(0)(0), tDvValues);
//        moris::print(tDvValues, "vals");
//        moris::print(tGeometryEngine.get_design_variable_interface()->get_my_local_global_map(), "map");

        // Solver
        std::shared_ptr< sol::SOL_Warehouse > tSolverWarehouse
        = std::make_shared<sol::SOL_Warehouse>( tModel->get_solver_interface() );

        moris::Cell< moris::Cell< moris::ParameterList > > tParameterlist( 7 );
        for( uint Ik = 0; Ik < 7; Ik ++)
        {
        tParameterlist( Ik ).resize(1);
        }

        tParameterlist( 0 )(0) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL );
        tParameterlist( 0 )( 0 ).set( "ifpack_prec_type", "ILU");

        tParameterlist( 1 )(0) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 2 )(0) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )(0).set( "NLA_max_iter", 2 );
        tParameterlist( 3 )(0) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )(0).set("NLA_DofTypes"      , "TEMP" );

        tParameterlist( 4 )(0) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 5 )(0) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )(0).set("TSA_DofTypes"      , "TEMP" );

        tParameterlist( 6 )(0) = moris::prm::create_solver_warehouse_parameterlist();

        tSolverWarehouse->set_parameterlist( tParameterlist );

        tSolverWarehouse->initialize();

        tsa::Time_Solver * tTimeSolver = tSolverWarehouse->get_main_time_solver();

        tTimeSolver->set_output( 0, nullptr );

        tModel->set_solver_warehouse_hack( tSolverWarehouse );

        tModel->perform();

        //--------------------------------------------------------------------------
        //                          Sensitivity analysis
        //--------------------------------------------------------------------------

        tGeometryEngine.communicate_requested_IQIs(tRequestedIQINames);

        tModel->perform_sensitivity_analysis();

        delete tModel;
        delete tInterpMesh;

    }
}

