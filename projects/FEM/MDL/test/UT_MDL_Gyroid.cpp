/*
 * UT_MDL_Gyroid.cpp
 *
 *  Created on: Oct 23, 2019
 *      Author: noel
 */

#include "catch.hpp"
#include "cl_Star.hpp"
#include "cl_Circle.hpp"
#include "cl_Plane.hpp"

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

#include "cl_FEM_IWG_Factory.hpp"             //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"           //FEM/INT/src

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

#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geom_Field.hpp"
#include "cl_GE_Geometry_Library.hpp"

#include "fn_norm.hpp"


namespace moris
{

// define free function for properties
 Matrix< DDRMat > tPropValConstFunc( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                     moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                     moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                     fem::Geometry_Interpolator              * aGeometryInterpolator )
 {
     return aParameters( 0 );
 }
 Matrix< DDRMat > tPropValFunc( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                fem::Geometry_Interpolator              * aGeometryInterpolator )
 {
     return aParameters( 0 ) + aParameters( 1 ) * aDofFieldInterpolator( 0 )->val();
 }
 Matrix< DDRMat > tPropDerFunc( moris::Cell< Matrix< DDRMat > >         & aParameters,
                                moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                fem::Geometry_Interpolator              * aGeometryInterpolator )
 {
     return aParameters( 1 ) * aDofFieldInterpolator( 0 )->N();
 }

TEST_CASE("MDL Gyroid","[MDL_Gyroid]")
{

    if(par_size() == 1)
    {
        uint tLagrangeMeshIndex = 0;
        // empty container for B-Spline meshes
        moris::Cell< moris::hmr::BSpline_Mesh_Base* > tBSplineMeshes;

        // create settings object
        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {10}, {10}, {10} } );
        tParameters.set_domain_dimensions( 5, 5, 5 );
        tParameters.set_domain_offset( 0.1, 0.1, 0.1 );
        tParameters.set_side_sets({ {1}, {2}, {3}, {4}, {5}, {6} });

        tParameters.set_bspline_truncation( true );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns( { {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_staircase_buffer( 1 );

        tParameters.set_initial_refinement( 1 );

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        // create the HMR object by passing the settings to the constructor
        moris::hmr::HMR tHMR( tParameters );

        // std::shared_ptr< Database >
        auto tDatabase = tHMR.get_database();

        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh01 = tHMR.create_mesh( tLagrangeMeshIndex );   // HMR Lagrange mesh
        //==============================
        std::shared_ptr< hmr::Field > tField = tMesh01->create_field( "gyroid", tLagrangeMeshIndex);

        tField->evaluate_scalar_function( moris::ge::getDistanceToGyroidsMassive );

        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( moris::ge::getDistanceToGyroidsMassive );
        }
            tDatabase->get_background_mesh()->save_to_vtk("Bachgroundmesh_initial_3x3x3.vtk");
        //==============================
        tDatabase->update_bspline_meshes();
        tDatabase->update_lagrange_meshes();

        // calculate T-Matrices etc
        tField->evaluate_scalar_function( moris::ge::getDistanceToGyroidsMassive );
        tDatabase->finalize();
//==============================
        tHMR.save_to_exodus( 0, "gyroid_general_geomEng.g" );

        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::ge::GEN_Geom_Field tFieldAsGeom(tField);

        moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};

        size_t tModelDimension = 3;
        moris::ge::GEN_Phase_Table  tPhaseTable( tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::GEN_Geometry_Engine  tGeometryEngine( tGeometryVector,tPhaseTable,tModelDimension );
        xtk::Model                  tXTKModel( tModelDimension,tInterpMesh.get(),tGeometryEngine );
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        //++++++++++++++++++++++++++++++++
//        xtk::Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

        // Declare the fields related to enrichment strategy in output options
//        Cell<std::string> tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
        //++++++++++++++++++++++++++++++++

               // output solution and meshes
                       xtk::Output_Options tOutputOptions1;
                       tOutputOptions1.mAddNodeSets = false;
                       tOutputOptions1.mAddSideSets = true;
                       tOutputOptions1.mAddClusters = false;

//                       tOutputOptions1.mRealElementExternalFieldNames = tEnrichmentFieldNames;

                       //++++++++++++++++++++++++++++++++
                       // add solution field to integration mesh
                       std::string tIntegSolFieldName1 = "solution";
                       tOutputOptions1.mRealNodeExternalFieldNames = {tIntegSolFieldName1};
                       //++++++++++++++++++++++++++++++++

                       moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);

//                       ++++++++++++++++++++++++++++++++
//                       tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tIntegMesh11);
//                       ++++++++++++++++++++++++++++++++



                       std::string tMeshOutputFile1 = "./output_general_geomEng.e";
                       tIntegMesh11->create_output_mesh(tMeshOutputFile1);
                       tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

       xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
       xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//==============================

       // place the pair in mesh manager
       mtk::Mesh_Manager tMeshManager;
       tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);

       //------------------------------------------------------------------------------
//               // create the properties
//               std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
//               tPropEMod1->set_parameters( { {{ 1.0 }} } );
//               tPropEMod1->set_val_function( tPropValConstFunc );
//
//               std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
//               tPropEMod2->set_parameters( { {{ 1.0 }} } );
//               tPropEMod2->set_val_function( tPropValConstFunc );
//
//               std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
//               tPropNu->set_parameters( { {{ 0.0 }} } );
//               tPropNu->set_val_function( tPropValConstFunc );
//
//               std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
//               tPropDirichlet->set_parameters( { {{ 0.0 }, { 0.0 }, { 0.0 }} } );
//               tPropDirichlet->set_val_function( tPropValConstFunc );
//
//               std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
//               tPropNeumann->set_parameters( {{{ 1.0 } , { 0.0 }, { 0.0 }}} );
//               tPropNeumann->set_val_function( tPropValConstFunc );
//
//               // define constitutive models
//               fem::CM_Factory tCMFactory;
//
//               std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
//               tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tCMStrucLinIso1->set_properties( { tPropEMod1, tPropNu } );
//               tCMStrucLinIso1->set_space_dim( 3 );
//
//               std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
//               tCMStrucLinIso2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tCMStrucLinIso2->set_properties( { tPropEMod2, tPropNu } );
//               tCMStrucLinIso2->set_space_dim( 3 );
//
//               // define the IWGs
//               fem::IWG_Factory tIWGFactory;
//
//               std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
//               tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//               tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tIWGBulk1->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::MASTER );
//
//               std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
//               tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//               tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tIWGBulk2->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );
//
//               std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
//               tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//               tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tIWGDirichlet->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::MASTER );
//               tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );
//
//               std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
//               tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//               tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );
//
//               std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
//               tIWGInterface->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
//               tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},mtk::Master_Slave::SLAVE );
//               tIWGInterface->set_constitutive_models( { tCMStrucLinIso2 } );
//               tIWGInterface->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::SLAVE );

       // create the properties
               std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
               tPropConductivity1->set_parameters( { {{ 1.0 }} } );
               tPropConductivity1->set_val_function( tPropValConstFunc );

               std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
               tPropConductivity2->set_parameters( { {{ 5.0 }} } );
               tPropConductivity2->set_val_function( tPropValConstFunc );

               std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
               tPropDirichlet->set_parameters( { {{ 5.0 }} } );
               tPropDirichlet->set_val_function( tPropValConstFunc );

               std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
               tPropNeumann->set_parameters( { {{ 20.0 }} } );
               tPropNeumann->set_val_function( tPropValConstFunc );

               std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
               tPropTempLoad1->set_parameters( { {{ 100.0 }} } );
               tPropTempLoad1->set_val_function( tPropValConstFunc );

               std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
               tPropTempLoad2->set_parameters( { {{ 100.0 }} } );
               tPropTempLoad2->set_val_function( tPropValConstFunc );

               // define constitutive models
               fem::CM_Factory tCMFactory;

               std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
               tCMDiffLinIso1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tCMDiffLinIso1->set_properties( { tPropConductivity1 } );
               tCMDiffLinIso1->set_space_dim( 3 );

               std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
               tCMDiffLinIso2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tCMDiffLinIso2->set_properties( { tPropConductivity2 } );
               tCMDiffLinIso2->set_space_dim( 3 );

               // define stabilization parameters
               fem::SP_Factory tSPFactory;
               std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
               tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
               tSPDirichletNitsche->set_properties( { tPropConductivity2 }, mtk::Master_Slave::MASTER );

               std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
               tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
               tSPNitscheInterface->set_properties( { tPropConductivity1 }, mtk::Master_Slave::MASTER );
               tSPNitscheInterface->set_properties( { tPropConductivity2 }, mtk::Master_Slave::SLAVE );

               std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
               tSPMasterWeightInterface->set_properties( { tPropConductivity1 }, mtk::Master_Slave::MASTER );
               tSPMasterWeightInterface->set_properties( { tPropConductivity2 }, mtk::Master_Slave::SLAVE );

               std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
               tSPSlaveWeightInterface->set_properties( { tPropConductivity1 }, mtk::Master_Slave::MASTER );
               tSPSlaveWeightInterface->set_properties( { tPropConductivity2 }, mtk::Master_Slave::SLAVE );

               // define the IWGs
               fem::IWG_Factory tIWGFactory;

               std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
               tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
               tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tIWGBulk1->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
               tIWGBulk1->set_properties( { tPropTempLoad1 }, mtk::Master_Slave::MASTER );

               std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
               tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
               tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tIWGBulk2->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
               tIWGBulk2->set_properties( { tPropTempLoad2 }, mtk::Master_Slave::MASTER );

               std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
               tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
               tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tIWGDirichlet->set_stabilization_parameters( { tSPDirichletNitsche } );
               tIWGDirichlet->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::MASTER );
               tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

               std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
               tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
               tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

               std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
               tIWGInterface->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
               tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::TEMP }},mtk::Master_Slave::SLAVE );
               tIWGInterface->set_stabilization_parameters( { tSPNitscheInterface, tSPMasterWeightInterface, tSPSlaveWeightInterface } );
               tIWGInterface->set_constitutive_models( { tCMDiffLinIso1 }, mtk::Master_Slave::MASTER );
               tIWGInterface->set_constitutive_models( { tCMDiffLinIso2 }, mtk::Master_Slave::SLAVE );


               // define set info
               fem::Set_User_Info tSetBulk1;
               tSetBulk1.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0") );
               tSetBulk1.set_set_type( fem::Element_Type::BULK );
               tSetBulk1.set_IWGs( { tIWGBulk2 } );

               fem::Set_User_Info tSetBulk2;
               tSetBulk2.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0") );
               tSetBulk2.set_set_type( fem::Element_Type::BULK );
               tSetBulk2.set_IWGs( { tIWGBulk2 } );

               fem::Set_User_Info tSetBulk3;
               tSetBulk3.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1") );
               tSetBulk3.set_set_type( fem::Element_Type::BULK );
               tSetBulk3.set_IWGs( { tIWGBulk1 } );

               fem::Set_User_Info tSetBulk4;
               tSetBulk4.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1") );
               tSetBulk4.set_set_type( fem::Element_Type::BULK );
               tSetBulk4.set_IWGs( { tIWGBulk1 } );

               fem::Set_User_Info tSetDirichlet;
               tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1") );
               tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
               tSetDirichlet.set_IWGs( { tIWGDirichlet } );

               fem::Set_User_Info tSetNeumann;
               tSetNeumann.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1") );
               tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
               tSetNeumann.set_IWGs( { tIWGNeumann } );

               std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);

               fem::Set_User_Info tSetInterface1;
               tSetInterface1.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index( tDblInterfaceSideSetName ) );
               tSetInterface1.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
               tSetInterface1.set_IWGs( { tIWGInterface } );

               // create a cell of set info
               moris::Cell< fem::Set_User_Info > tSetInfo( 7 );
               tSetInfo( 0 ) = tSetBulk1;
               tSetInfo( 1 ) = tSetBulk2;
               tSetInfo( 2 ) = tSetBulk3;
               tSetInfo( 3 ) = tSetBulk4;
               tSetInfo( 4 ) = tSetDirichlet;
               tSetInfo( 5 ) = tSetNeumann;
               tSetInfo( 6 ) = tSetInterface1;

               // create model
               mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                                      0,
                                                      tSetInfo );

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

//               tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
//               tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
               tNonlinearSolver.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
               tTimeSolver.set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );

               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               // STEP 4: Solve and check
               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

               tTimeSolver.solve();


               Matrix<DDRMat> tFullSol;
               tTimeSolver.get_full_solution(tFullSol);

//               print(tFullSol,"tFullSol");


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


               std::string tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_plane_hole_3d_l" + std::to_string(1) + "_b"+std::to_string(0)+".e";
               tIntegMesh1->create_output_mesh(tMeshOutputFile);

       delete tIntegMesh1;

    }


}/* END_TEST_CASE */

}/* END_MORIS_NAMESPACE */




