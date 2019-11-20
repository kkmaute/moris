/*
 * UT_MDL_XTK_HMR_2D.cpp
 *
 *  Created on: Sep 18, 2019
 *      Author: schmidt
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

#include "cl_FEM_NodeProxy.hpp"                //FEM/INT/src
#include "cl_FEM_ElementProxy.hpp"             //FEM/INT/src
#include "cl_FEM_Node_Base.hpp"                //FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"          //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

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

#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geometry.hpp"
#include "../projects/GEN/src/ripped/geometry/cl_GEN_Circle.hpp"
#include "../projects/GEN/src/ripped/geometry/cl_GEN_Geom_Field.hpp"

#include "fn_norm.hpp"


namespace moris
{

Matrix< DDRMat >
exactTempFunc(moris::Cell< Matrix< DDRMat > >         & aCoeff,
        moris::Cell< fem::Field_Interpolator* > & aFieldInterpolator,
        fem::Geometry_Interpolator             * aGeometryInterpolator )
        {
    Matrix< DDRMat > tCoord = aGeometryInterpolator->valx();
    real xcoord = tCoord(0);
    real ycoord = tCoord(1);

    real rad = std::pow (  std::pow( xcoord - 0, 2.0)
    + std::pow( ycoord - 0, 2.0), 0.5);

    return {{(1.0/3.0)*(1.0/rad-0.501)}};
        }

moris::real LvlSetLin(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tOffset = 200;

    return    aPoint(0) - 0.317 * aPoint(1) - tOffset;
}

moris::real LvlSetCircle_2D(const moris::Matrix< moris::DDRMat > & aPoint )
{
    return    std::sqrt( aPoint( 0 ) * aPoint( 0 ) + aPoint( 1 ) * aPoint( 1 ) ) - 0.2505;
}

Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                    moris::Cell< fem::Field_Interpolator* > & aDofFieldInterpolator,
                                    moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                    fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aCoeff( 0 );
}

moris::real
LevelSetFunction_star1( const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real tPhi = std::atan2( aPoint( 0 ), aPoint( 1 ) );

    moris::real tLevelSetVaue = 0.501 + 0.1 * std::sin( 5 * tPhi ) - std::sqrt( std::pow( aPoint( 0 ), 2 ) + std::pow( aPoint( 1 ), 2 ) );

    return -tLevelSetVaue;
}

moris::real
Plane4MatMDL1(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXC = 0.1;
    moris::real mYC = 0.1;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return ( mNx*( aPoint(0)-mXC ) + mNy*( aPoint(1)-mYC ) );
}

moris::real
Circle4MatMDL(const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXCenter = 0.01;
    moris::real mYCenter = 0.01;
    moris::real mRadius = 0.47334;


    return  (aPoint(0) - mXCenter) * (aPoint(0) - mXCenter)
                    + (aPoint(1) - mYCenter) * (aPoint(1) - mYCenter)
                    - (mRadius * mRadius);
}

TEST_CASE("2D XTK WITH HMR Struc Interface 2D","[XTK_HMR_Struc_Interface_2D]")
{
    if(par_size()<=1)
    {
        uint tLagrangeMeshIndex = 0;
        std::string tFieldName = "Cylinder";

         hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

         tParameters.set( "number_of_elements_per_dimension", "20, 20");
         tParameters.set( "domain_dimensions", "2, 2" );
         tParameters.set( "domain_offset", "-1.0, -1.0" );
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
         tParameters.set( "initial_refinement", 0 );

         tParameters.set( "use_multigrid", 0 );
         tParameters.set( "severity_level", 2 );

         hmr::HMR tHMR( tParameters );

        //initial refinement
         tHMR.perform_initial_refinement( 0 );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

       //  create field
       std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

       tField->evaluate_scalar_function( LvlSetCircle_2D );

       for( uint k=0; k<2; ++k )
       {
           tHMR.flag_surface_elements_on_working_pattern( tField );
           tHMR.perform_refinement_based_on_working_pattern( 0 );

           tField->evaluate_scalar_function( LvlSetCircle_2D );
       }

         tHMR.finalize();

       tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

          std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        //----------------------------------------------------------------------------------------------
//       std::string tPrefix = std::getenv("MORISROOT");
//       std::string tMeshFileName = tPrefix + "build/3x1.g";
//       std::cout<<"Mesh input name = "<<tMeshFileName<<std::endl;
//
//       moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
//       std::string tFieldName1 = "Temp_Field";
//       tNodeField1.set_field_name( tFieldName1 );
//       tNodeField1.set_field_entity_rank( EntityRank::NODE );
//
//       // Initialize field information container
//       moris::mtk::MtkFieldsInfo tFieldsInfo;
//
//       // Place the node field into the field info container
//       add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
//
//       // Declare some supplementary fields
//       mtk::MtkMeshData tMeshData;
//       tMeshData.FieldsInfo = &tFieldsInfo;
//
//
//       // construct the mesh data
//       mtk::Interpolation_Mesh* tInterpolationMesh = mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tMeshData );
        //-----------------------------------------------------------------------------------------------

        moris::ge::Circle tCircle( 0.2501, 0.0, 0.0 );
//        moris::Matrix<moris::DDRMat> tCenters = {{ 0.151,0.1 }};
//        moris::Matrix<moris::DDRMat> tNormals = {{ 1.0,0.0 }};
//        xtk::Plane<2> tPlane(tCenters,tNormals);

        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tCircle,tPhaseTable, 2);

         xtk::Model tXTKModel(2, tInterpolationMesh.get(), tGeometryEngine);
//       xtk::Model tXTKModel(2, tInterpolationMesh, tGeometryEngine);

        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

        // get meshes
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        tEnrInterpMesh.print_enriched_cells();
        tEnrIntegMesh.print_double_side_sets(2);

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //-----------------------------------------------------------------------------
//        // output solution and meshes
//                xtk::Output_Options tOutputOptions1;
//                tOutputOptions1.mAddNodeSets = false;
//                tOutputOptions1.mAddSideSets = false;
//                tOutputOptions1.mAddClusters = false;
//
//
//                moris::mtk::Integration_Mesh* tIntegMesh11 = tXTKModel.get_output_mesh(tOutputOptions1);
//
//                for(moris::uint i = 0; i < tIntegMesh11->get_num_entities(EntityRank::NODE); i++)
//                {
//                    moris::moris_id tID = tIntegMesh11->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//                }
//
//
//                std::string tMeshOutputFile1 = "./mdl_exo/stk_xtk_inv_ilu_quad_bspline1.e";
//                tIntegMesh11->create_output_mesh(tMeshOutputFile1);
        //-----------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
        tPropEMod1->set_parameters( { {{ 1.0 }} } );
        tPropEMod1->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
        tPropEMod2->set_parameters( { {{ 10.0 }} } );
        tPropEMod2->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
        tPropNu->set_parameters( { {{ 0.0 }} } );
        tPropNu->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 0.0 }, { 0.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( {{{ 1.0 } , { 0.0 }}} );
        tPropNeumann->set_val_function( tConstValFunction );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tCMStrucLinIso1->set_properties( { tPropEMod1, tPropNu } );
        tCMStrucLinIso1->set_space_dim( 2 );

        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMStrucLinIso2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tCMStrucLinIso2->set_properties( { tPropEMod2, tPropNu } );
        tCMStrucLinIso2->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_properties( { tPropEMod1 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
        tSPNitscheInterface->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
        tSPNitscheInterface->set_properties( { tPropEMod1 }, mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
        tSPMasterWeightInterface->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
        tSPMasterWeightInterface->set_properties( { tPropEMod1 }, mtk::Master_Slave::SLAVE );

        std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
        tSPSlaveWeightInterface->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
        tSPSlaveWeightInterface->set_properties( { tPropEMod1 }, mtk::Master_Slave::SLAVE );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGBulk1->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGBulk2->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGDirichlet->set_stabilization_parameters( { tSPDirichletNitsche } );
        tIWGDirichlet->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::MASTER );
        tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGInterface = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
        tIWGInterface->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }},mtk::Master_Slave::SLAVE );
        tIWGInterface->set_stabilization_parameters( { tSPNitscheInterface, tSPMasterWeightInterface, tSPSlaveWeightInterface } );
        tIWGInterface->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );
        tIWGInterface->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::SLAVE );

        // create a list of active block-sets
        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );
        std::string tDblInterfaceSideSetName = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);

        //        // create a list of active block-sets
        ////       moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("block_1_c_p0"),
        ////                                                tEnrIntegMesh.get_block_set_index("block_1_n_p0"),
        ////                                                tEnrIntegMesh.get_block_set_index("block_1_c_p1"),
        ////                                                tEnrIntegMesh.get_block_set_index("block_1_n_p1"),
        ////                                                tEnrIntegMesh.get_side_set_index("surface_4_n_p0"),
        ////                                                tEnrIntegMesh.get_side_set_index("surface_2_n_p1"),
        ////                                                tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName)};
        //
        //         moris::Cell< moris_index >  tSetList = { tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p0"),
        //                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p0"),
        //                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1"),
        //                                                  tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1"),
        //                                                  tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1"),
        //                                                  tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1"),
        //                                                  tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName)};

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

        fem::Set_User_Info tSetInterface;
        tSetInterface.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName) );
        tSetInterface.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
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
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
        //            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 1000;
        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 5;

        //    tLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Epetra );
        //    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 2;
        //        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
        //        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        //        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        tNonlinearSolver.solve( tNonlinearProblem );

        std::cout<<" Solution Vector "<<std::endl;
        tNonlinearProblem->get_full_vector()->print();

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = false;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldNameUX = "UX";
        std::string tIntegSolFieldNameUY = "UY";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldNameUX, tIntegSolFieldNameUY};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
        Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );


        Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tSTKIntegSolUX(i) = tIntegSolUX(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
            tSTKIntegSolUY(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
        }

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);

        std::string tMeshOutputFile = "./mdl_exo/stk_xtk_inv_ilu_quad_bspline.e";

        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tIntegMesh1;

        delete tModel;
    }
}

TEST_CASE("2D XTK WITH HMR Struc 2D","[XTK_HMR_Struc_2D]")
{
    if(par_size()<=1)
    {
        uint tLagrangeMeshIndex = 0;
        std::string tFieldName = "Cylinder";

        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "20, 20");
        tParameters.set( "domain_dimensions", "2, 2" );
        tParameters.set( "domain_offset", "-1.0, -1.0" );
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
        tParameters.set( "initial_refinement", 0 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );

        hmr::HMR tHMR( tParameters );

        // initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //// create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( LvlSetCircle_2D );
        //
        for( uint k=0; k<2; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( LvlSetCircle_2D );
        }

        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_2d.e" );

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        moris::ge::Circle tCircle( 0.4501, 0.0, 0.0 );

        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tCircle,tPhaseTable, 2);

        xtk::Model tXTKModel(2, tInterpolationMesh.get(), tGeometryEngine);

        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment( EntityRank::NODE, 0 );

        // get meshes
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
        tPropEMod->set_parameters( { {{ 1000000.0 }} } );
        tPropEMod->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
        tPropNu->set_parameters( { {{ 0.3 }} } );
        tPropNu->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 0.0 }, { 0.0 }} } );
        tPropDirichlet->set_val_function( tConstValFunction );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( {{{ 1000.0 } , { 100.0 }}} );
        tPropNeumann->set_val_function( tConstValFunction );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tCMStrucLinIso->set_properties( { tPropEMod, tPropNu } );
        tCMStrucLinIso->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
        tSPDirichletNitsche->set_properties( { tPropEMod }, mtk::Master_Slave::MASTER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGBulk->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGBulk->set_constitutive_models( { tCMStrucLinIso }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
        tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGDirichlet->set_stabilization_parameters( { tSPDirichletNitsche } );
        tIWGDirichlet->set_constitutive_models( { tCMStrucLinIso }, mtk::Master_Slave::MASTER );
        tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

        // create a list of active block-sets
        std::string tInterfaceSideSetName = tEnrIntegMesh.get_interface_side_set_name( 0, 0, 1 );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p1") );
        tSetBulk1.set_set_type( fem::Element_Type::BULK );
        tSetBulk1.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1") );
        tSetBulk2.set_set_type( fem::Element_Type::BULK );
        tSetBulk2.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p1") );
        tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p1") );
        tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetDirichlet;
        tSetInfo( 3 ) = tSetNeumann;

        // create model
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
        //            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 1000;
        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
        tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 5;

        //    tLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Epetra );
        //    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 2;
        //        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
        //        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        //        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        tNonlinearSolver.solve( tNonlinearProblem );

        //        Matrix<DDRMat> tFullSol;
        tNonlinearProblem->get_full_vector()->print();

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = false;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldNameUX = "UX";
        std::string tIntegSolFieldNameUY = "UY";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldNameUX, tIntegSolFieldNameUY};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
        Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );

        //    print(tIntegSolUX,"tIntegSolUX");
        //    print(tIntegSolUY,"tIntegSolUY");

        Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tSTKIntegSolUX(i) = tIntegSolUX(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
            tSTKIntegSolUY(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
        }

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);

        //    Matrix<DDRMat> tFullSol;
        //    tNonlinearSolver.get_full_solution(tFullSol);
        //
        //    print(tFullSol,"tFullSol");

        std::string tMeshOutputFile = "./mdl_exo/stk_xtk_inv_ilu_quad_bspline.e";

        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tIntegMesh1;

        delete tModel;
    }
}

#include <armadillo>


TEST_CASE("Arma Test","[Arma_Test]")
{
    uint tIter = 100;
    arma::Mat<double> tMatrix;
    tMatrix.ones(1,64);

    arma::Mat<double> tMatrixS;
    arma::Mat<double> tMatrixSS;

    arma::Mat<double> tVector;
    tVector.ones(320);

    arma::SpMat<double> tSpMat11;
    arma::SpMat<double> tSpMat22;

//    tMatrixS( arma::span( 0, 0 ), arma::span( 0, 63 ) ) = tMatrix;
//    tMatrixS( arma::span( 1, 1 ), arma::span( 64, 127 ) ) = tMatrix;
//    tMatrixS( arma::span( 2, 2 ), arma::span( 128, 191 ) ) = tMatrix;
//    tMatrixS( arma::span( 3, 3 ), arma::span( 192, 255 ) ) = tMatrix;
//    tMatrixS( arma::span( 4, 4 ), arma::span( 256, 319 ) ) = tMatrix;
//
//    tSpMat11( arma::span( 0, 0 ), arma::span( 0, 63 ) ) = tMatrix;
//    tSpMat11( arma::span( 1, 1 ), arma::span( 64, 127 ) ) = tMatrix;
//    tSpMat11( arma::span( 2, 2 ), arma::span( 128, 191 ) ) = tMatrix;
//    tSpMat11( arma::span( 3, 3 ), arma::span( 192, 255 ) ) = tMatrix;
//    tSpMat11( arma::span( 4, 4 ), arma::span( 256, 319 ) ) = tMatrix;

    tSpMat11.sprandu(1000, 5, 0.2);
    tSpMat22.sprandu(5, 1000, 0.2);

    tMatrixS = arma::mat(tSpMat11);
    tMatrixSS = arma::mat(tSpMat22);

//    sp_mat C = 2*B;
//    sp_mat D = A*C;
//
//    tSpMat11.print();


//    tMatrixSS.eye();

    tic tTimer_Dense;
//    arma::Mat<double> tMatrixR1(320,320);

    for(uint Ik = 0; Ik<tIter; Ik++)
    {
        arma::Mat<double> tMatrixR1;

        tMatrixR1 =  tMatrixS   * tMatrixSS ;
    }

    real tElapsedTime_Dense = tTimer_Dense.toc<moris::chronos::milliseconds>().wall;

    MORIS_LOG( " dense Multiplication on processor %u took %5.3f seconds.\n", ( uint ) par_rank(), ( double ) tElapsedTime_Dense / 1000);

    arma::SpMat<double> tSpMat = arma::sp_mat(tMatrixS);

    arma::SpMat<double> tSpMatS = arma::sp_mat(tMatrixSS);

    tic tTimer_Sparse;
//    arma::Mat<double> tMatrixR2(320,320);
    for(uint Ik = 0; Ik<tIter; Ik++)
    {

        arma::Mat<double> tMatrixR2;
        tMatrixR2 = tSpMat11 * tSpMat22;
    }

    real tElapsedTime_Sparse = tTimer_Sparse.toc<moris::chronos::milliseconds>().wall;

    MORIS_LOG( " Sparse Multiplication on processor %u took %5.3f seconds.\n", ( uint ) par_rank(), ( double ) tElapsedTime_Sparse / 1000);

    tic tTimer_Sparse_2;

//    arma::SpMat<double> tMatrixR3(320,320);

    for(uint Ik = 0; Ik<tIter; Ik++)
    {
        arma::SpMat<double> tMatrixR3;

        tMatrixR3 = tSpMat11 * tSpMat22;
    }

    real tElapsedTime_Sparse_2 = tTimer_Sparse_2.toc<moris::chronos::milliseconds>().wall;

    MORIS_LOG( " Sparse 2 Multiplication on processor %u took %5.3f seconds.\n", ( uint ) par_rank(), ( double ) tElapsedTime_Sparse_2 / 1000);

//    tMatrixS.print();

}

TEST_CASE("2D XTK WITH HMR Struc Interface 3D","[XTK_HMR_Struc_Interface_3D]")
{
    if(par_size()<=1)
    {
        uint tLagrangeMeshIndex = 0;
        std::string tFieldName = "Cylinder";

        uint tSpatialDimension = 3;

         hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

         tParameters.set( "number_of_elements_per_dimension", "22, 8, 2");
         tParameters.set( "domain_dimensions", "6, 2, 1" );
         tParameters.set( "domain_offset", "-3.0, -1.0, -0.5" );
         tParameters.set( "domain_sidesets", "1,2,3,4,5,6" );
         tParameters.set( "lagrange_output_meshes", "0" );

         tParameters.set( "lagrange_orders", "1" );
         tParameters.set( "lagrange_pattern", "0" );
         tParameters.set( "bspline_orders", "1" );
         tParameters.set( "bspline_pattern", "0" );

         tParameters.set( "lagrange_to_bspline", "0" );

         tParameters.set( "truncate_bsplines", 1 );
         tParameters.set( "refinement_buffer", 3 );
         tParameters.set( "staircase_buffer", 3 );
         tParameters.set( "initial_refinement", 0 );

         tParameters.set( "use_multigrid", 0 );
         tParameters.set( "severity_level", 2 );

         hmr::HMR tHMR( tParameters );

        //initial refinement
         tHMR.perform_initial_refinement( 0 );

         std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

       //  create field
       Cell<std::shared_ptr< moris::hmr::Field > > tHMRFields;
       tHMRFields.resize(2);

       // create field
       tHMRFields(0) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
       tHMRFields(1) = tMesh->create_field( "Geom", tLagrangeMeshIndex );

       tHMRFields(0)->evaluate_scalar_function( LevelSetFunction_star1 );
       tHMRFields(1)->evaluate_scalar_function( Plane4MatMDL1 );

       for( uint k=0; k<2; ++k )
       {
           tHMR.flag_surface_elements_on_working_pattern( tHMRFields(0) );
           tHMR.flag_surface_elements_on_working_pattern( tHMRFields(1) );

           tHMR.perform_refinement_based_on_working_pattern( 0 );

           tHMRFields(0)->evaluate_scalar_function( LevelSetFunction_star1 );
           tHMRFields(1)->evaluate_scalar_function( Plane4MatMDL1 );
       }

      tHMR.finalize();

//      tHMR.save_to_exodus( 0, "./xtk_exo/mdl_xtk_hmr_3d.e" );

      std::shared_ptr< hmr::Interpolation_Mesh_HMR > tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

      moris::ge::GEN_Geom_Field tCircleFieldAsGeom(tHMRFields(0));
      moris::ge::GEN_Geom_Field tPlaneFieldAsGeom2(tHMRFields(1));
      moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = {&tCircleFieldAsGeom,&tPlaneFieldAsGeom2};

      moris::ge::GEN_Phase_Table     tPhaseTable (tGeometryVector.size(),  Phase_Table_Structure::EXP_BASE_2);
      moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable,tSpatialDimension);
      xtk::Model           tXTKModel(tSpatialDimension,tInterpMesh.get(),tGeometryEngine);
      tXTKModel.mVerbose = false;

      Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
      tXTKModel.decompose(tDecompositionMethods);
      tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

      // get meshes
      xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
      xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

      // place the pair in mesh manager
      mtk::Mesh_Manager tMeshManager;
      tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

      //------------------------------------------------------------------------------
      // create the properties
      std::shared_ptr< fem::Property > tPropEMod1 = std::make_shared< fem::Property >();
      tPropEMod1->set_parameters( { {{ 1.0 }} } );
      tPropEMod1->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropEMod1bis = std::make_shared< fem::Property >();
      tPropEMod1bis->set_parameters( { {{ 1.0 }} } );
      tPropEMod1bis->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropEMod2 = std::make_shared< fem::Property >();
      tPropEMod2->set_parameters( { {{ 1.0 }} } );
      tPropEMod2->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropEMod2bis = std::make_shared< fem::Property >();
      tPropEMod2bis->set_parameters( { {{ 1.0 }} } );
      tPropEMod2bis->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
      tPropNu->set_parameters( { {{ 0.0 }} } );
      tPropNu->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropNubis = std::make_shared< fem::Property >();
      tPropNubis->set_parameters( { {{ 0.0 }} } );
      tPropNubis->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
      tPropDirichlet->set_parameters( { {{0.0}, {0.0}, {0.0}} } );
      tPropDirichlet->set_val_function( tConstValFunction );

      std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
      tPropNeumann->set_parameters( { {{1.0}, {0.0}, {0.0}} } );
      tPropNeumann->set_val_function( tConstValFunction );

      // define constitutive models
      fem::CM_Factory tCMFactory;

      std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
      tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
      tCMStrucLinIso1->set_properties( { tPropEMod1, tPropNu } );
      tCMStrucLinIso1->set_space_dim( 3 );

      std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1bis = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
      tCMStrucLinIso1bis->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
      tCMStrucLinIso1bis->set_properties( { tPropEMod1bis, tPropNubis } );
      tCMStrucLinIso1bis->set_space_dim( 3 );

      std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
      tCMStrucLinIso2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
      tCMStrucLinIso2->set_properties( { tPropEMod2, tPropNu } );
      tCMStrucLinIso2->set_space_dim( 3 );

      std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2bis = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
      tCMStrucLinIso2bis->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
      tCMStrucLinIso2bis->set_properties( { tPropEMod2bis, tPropNubis } );
      tCMStrucLinIso2bis->set_space_dim( 3 );

      // define stabilization parameters
      fem::SP_Factory tSPFactory;
      std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
      tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
      tSPDirichletNitsche->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );

      std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface1 = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
      tSPNitscheInterface1->set_parameters( { {{ 1.0 }} } );
      tSPNitscheInterface1->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
      tSPNitscheInterface1->set_properties( { tPropEMod2bis }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface1 = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
      tSPMasterWeightInterface1->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
      tSPMasterWeightInterface1->set_properties( { tPropEMod2bis }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface1 = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
      tSPSlaveWeightInterface1->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
      tSPSlaveWeightInterface1->set_properties( { tPropEMod2bis }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface2 = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
      tSPNitscheInterface2->set_parameters( { {{ 1.0 }} } );
      tSPNitscheInterface2->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
      tSPNitscheInterface2->set_properties( { tPropEMod1 }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface2 = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
      tSPMasterWeightInterface2->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
      tSPMasterWeightInterface2->set_properties( { tPropEMod1 }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface2 = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
      tSPSlaveWeightInterface2->set_properties( { tPropEMod2 }, mtk::Master_Slave::MASTER );
      tSPSlaveWeightInterface2->set_properties( { tPropEMod1 }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface3 = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
      tSPNitscheInterface3->set_parameters( { {{ 1.0 }} } );
      tSPNitscheInterface3->set_properties( { tPropEMod1 }, mtk::Master_Slave::MASTER );
      tSPNitscheInterface3->set_properties( { tPropEMod1bis }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface3 = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
      tSPMasterWeightInterface3->set_properties( { tPropEMod1 }, mtk::Master_Slave::MASTER );
      tSPMasterWeightInterface3->set_properties( { tPropEMod1bis }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface3 = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
      tSPSlaveWeightInterface3->set_properties( { tPropEMod1 }, mtk::Master_Slave::MASTER );
      tSPSlaveWeightInterface3->set_properties( { tPropEMod1bis }, mtk::Master_Slave::SLAVE );

      // define the IWGs
      fem::IWG_Factory tIWGFactory;

      std::shared_ptr< fem::IWG > tIWGBulk1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
      tIWGBulk1->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGBulk1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGBulk1->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::MASTER );

      std::shared_ptr< fem::IWG > tIWGBulk2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
      tIWGBulk2->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGBulk2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGBulk2->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );

      std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
      tIWGDirichlet->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGDirichlet->set_stabilization_parameters( { tSPDirichletNitsche } );
      tIWGDirichlet->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );
      tIWGDirichlet->set_properties( { tPropDirichlet }, mtk::Master_Slave::MASTER );

      std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
      tIWGNeumann->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGNeumann->set_properties( { tPropNeumann }, mtk::Master_Slave::MASTER );

      std::shared_ptr< fem::IWG > tIWGInterface1 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
      tIWGInterface1->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGInterface1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGInterface1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},mtk::Master_Slave::SLAVE );
      tIWGInterface1->set_stabilization_parameters( { tSPNitscheInterface1, tSPMasterWeightInterface1, tSPSlaveWeightInterface1 } );
      tIWGInterface1->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );
      tIWGInterface1->set_constitutive_models( { tCMStrucLinIso2bis }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::IWG > tIWGInterface2 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
      tIWGInterface2->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGInterface2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGInterface2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},mtk::Master_Slave::SLAVE );
      tIWGInterface2->set_stabilization_parameters( { tSPNitscheInterface2, tSPMasterWeightInterface2, tSPSlaveWeightInterface2 } );
      tIWGInterface2->set_constitutive_models( { tCMStrucLinIso2 }, mtk::Master_Slave::MASTER );
      tIWGInterface2->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::SLAVE );

      std::shared_ptr< fem::IWG > tIWGInterface3 = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
      tIWGInterface3->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );
      tIWGInterface3->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
      tIWGInterface3->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }},mtk::Master_Slave::SLAVE );
      tIWGInterface3->set_stabilization_parameters( { tSPNitscheInterface3, tSPMasterWeightInterface3, tSPSlaveWeightInterface3 } );
      tIWGInterface3->set_constitutive_models( { tCMStrucLinIso1 }, mtk::Master_Slave::MASTER );
      tIWGInterface3->set_constitutive_models( { tCMStrucLinIso1bis }, mtk::Master_Slave::SLAVE );

      // create a list of active block-sets
      std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,1);
      std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name(0,2);
      std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name(1,3);
      std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name(2,3);

      std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01)<<std::endl;
      std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02)<<std::endl;

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
      tSetBulk3.set_IWGs( { tIWGBulk2 } );

      fem::Set_User_Info tSetBulk4;
      tSetBulk4.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p1") );
      tSetBulk4.set_set_type( fem::Element_Type::BULK );
      tSetBulk4.set_IWGs( { tIWGBulk2 } );

      fem::Set_User_Info tSetBulk5;
      tSetBulk5.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p2") );
      tSetBulk5.set_set_type( fem::Element_Type::BULK );
      tSetBulk5.set_IWGs( { tIWGBulk1 } );

      fem::Set_User_Info tSetBulk6;
      tSetBulk6.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p2") );
      tSetBulk6.set_set_type( fem::Element_Type::BULK );
      tSetBulk6.set_IWGs( { tIWGBulk1 } );

      fem::Set_User_Info tSetBulk7;
      tSetBulk7.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_c_p3") );
      tSetBulk7.set_set_type( fem::Element_Type::BULK );
      tSetBulk7.set_IWGs( { tIWGBulk1 } );

      fem::Set_User_Info tSetBulk8;
      tSetBulk8.set_mesh_index( tEnrIntegMesh.get_block_set_index("HMR_dummy_n_p3") );
      tSetBulk8.set_set_type( fem::Element_Type::BULK );
      tSetBulk8.set_IWGs( { tIWGBulk1 } );

      fem::Set_User_Info tSetDirichlet;
      tSetDirichlet.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_4_n_p2") );
      tSetDirichlet.set_set_type( fem::Element_Type::SIDESET );
      tSetDirichlet.set_IWGs( { tIWGDirichlet } );

      fem::Set_User_Info tSetNeumann;
      tSetNeumann.set_mesh_index( tEnrIntegMesh.get_side_set_index("SideSet_2_n_p3") );
      tSetNeumann.set_set_type( fem::Element_Type::SIDESET );
      tSetNeumann.set_IWGs( { tIWGNeumann } );

      fem::Set_User_Info tSetInterface1;
      tSetInterface1.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName01) );
      tSetInterface1.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
      tSetInterface1.set_IWGs( { tIWGInterface1 } );

      fem::Set_User_Info tSetInterface2;
      tSetInterface2.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName02) );
      tSetInterface2.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
      tSetInterface2.set_IWGs( { tIWGInterface2 } );

      fem::Set_User_Info tSetInterface3;
      tSetInterface3.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName13) );
      tSetInterface3.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
      tSetInterface3.set_IWGs( { tIWGInterface2 } );

      fem::Set_User_Info tSetInterface4;
      tSetInterface4.set_mesh_index( tEnrIntegMesh.get_double_sided_set_index(tDblInterfaceSideSetName23) );
      tSetInterface4.set_set_type( fem::Element_Type::DOUBLE_SIDESET );
      tSetInterface4.set_IWGs( { tIWGInterface3 } );

      // create a cell of set info
      moris::Cell< fem::Set_User_Info > tSetInfo( 14 );
      tSetInfo( 0 ) = tSetBulk1;
      tSetInfo( 1 ) = tSetBulk2;
      tSetInfo( 2 ) = tSetBulk3;
      tSetInfo( 3 ) = tSetBulk4;
      tSetInfo( 4 ) = tSetBulk5;
      tSetInfo( 5 ) = tSetBulk6;
      tSetInfo( 6 ) = tSetBulk7;
      tSetInfo( 7 ) = tSetBulk8;
      tSetInfo( 8 ) = tSetDirichlet;
      tSetInfo( 9 ) = tSetNeumann;
      tSetInfo( 10 ) = tSetInterface1;
      tSetInfo( 11 ) = tSetInterface2;
      tSetInfo( 12 ) = tSetInterface3;
      tSetInfo( 13 ) = tSetInterface4;

      // create model
      mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                             0,
                                             tSetInfo );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::AZTEC_IMPL );
        //            std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( SolverType::PETSC );

        tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
        tLinearSolverAlgorithm->set_param("AZ_max_iter") = 5000;
        tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;
        tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
        tLinearSolverAlgorithm->set_param("AZ_ilut_fill") = 10.0;
        tLinearSolverAlgorithm->set_param("rel_residual") = 1e-8;

//            tLinearSolverAlgorithm->set_param("Use_ML_Prec") = true;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Epetra );
        //    NLA::Nonlinear_Problem * tNonlinearProblem =  new NLA::Nonlinear_Problem( tModel->get_solver_interface(), 0, true, MapType::Petsc );

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 2;
        tNonlinearSolverAlgorithm->set_param("NLA_rel_residual") = 1e-4;
        //        tNonlinearSolverAlgorithm->set_param("NLA_hard_break") = false;
        //        tNonlinearSolverAlgorithm->set_param("NLA_max_lin_solver_restarts") = 2;
        //        tNonlinearSolverAlgorithm->set_param("NLA_rebuild_jacobian") = true;

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;

        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        tNonlinearSolver.solve( tNonlinearProblem );

        std::cout<<" Solution Vector "<<std::endl;
        tNonlinearProblem->get_full_vector()->print();

        // output solution and meshes
        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = false;
        tOutputOptions.mAddClusters = false;

        // add solution field to integration mesh
        std::string tIntegSolFieldNameUX = "UX";
        std::string tIntegSolFieldNameUY = "UY";
        std::string tIntegSolFieldNameUZ = "UZ";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldNameUX, tIntegSolFieldNameUY, tIntegSolFieldNameUZ};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // Write to Integration mesh for visualization
        Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
        Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );
        Matrix<DDRMat> tIntegSolUZ = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UZ );

        Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        Matrix<DDRMat> tSTKIntegSolUZ(tIntegMesh1->get_num_entities(EntityRank::NODE),1);

        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        {
            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
            tSTKIntegSolUX(i) = tIntegSolUX(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
            tSTKIntegSolUY(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
            tSTKIntegSolUZ(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
        }

        // add solution field to integration mesh
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);
        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUZ,EntityRank::NODE,tSTKIntegSolUZ);

        std::string tMeshOutputFile = "./mdl_exo/hmr_xtk_linear_elastic_3D.e";

        tIntegMesh1->create_output_mesh(tMeshOutputFile);

        delete tIntegMesh1;

        delete tModel;
    }
}
}



