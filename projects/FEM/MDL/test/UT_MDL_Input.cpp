/*
 * UT_MDL_Input.cpp
 *
 *  Created on: Oct 23, 2019
 *      Author: noel
 */

#include "catch.hpp"

#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "typedefs.hpp"
#include "paths.hpp"

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

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "cl_FEM_IWG_Factory.hpp"             //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"           //FEM/INT/src
#include "cl_FEM_Model.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

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

#include "cl_GEN_Plane.hpp"

#include "fn_norm.hpp"


namespace moris
{

// define free function for properties
void tPropValConstFunc
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}
void tPropValFunc
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    moris::fem::Field_Interpolator * tFI
    = aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP );
    aPropMatrix = aParameters( 0 ) + aParameters( 1 ) * tFI->val();
}
void tPropDerFunc
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    moris::fem::Field_Interpolator * tFI
    = aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP );
    aPropMatrix = aParameters( 1 ) * tFI->N();
}

// define free function for geometry
moris::real LevelSetPlaneFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    moris::real mXn = 0;
    moris::real mYn = 0;
    moris::real mZn = 1.0;
    moris::real mXc = 1.0;
    moris::real mYc = 1.0;
    moris::real mZc = 1.4;
    return mXn*(aPoint(0)-mXc) + mYn*(aPoint(1)-mYc) + mZn*(aPoint(2)-mZc);
}

/* This UT tests:
* -
* -
*/
TEST_CASE("MDL Input","[MDL_Input]")
{
    if(par_size() == 1)
    {
        std::string tFieldName = "Circle";


        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {1}, {1}, {4} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {2} });
        tParameters.set_domain_offset({ {0.0}, {0.0}, {0.0} });
        tParameters.set_bspline_truncation( true );
        tParameters.set_side_sets({ {5}, {6} });

        tParameters.set_output_meshes( {{ {0} }} );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 2);

        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );

        tField->evaluate_scalar_function( LevelSetPlaneFunction );

        for( uint k=0; k<0; ++k )
        {
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );

            tField->evaluate_scalar_function( LevelSetPlaneFunction );
        }

        tHMR.finalize();

        tHMR.save_to_exodus( 0, "./mdl_exo/mdl_input.e" );

        hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );

        moris::Cell< std::shared_ptr<moris::ge::Geometry> > tGeometryVector(1);
        tGeometryVector(0) = std::make_shared<moris::ge::Plane>(1.0, 1.0, 1.4, 0.0, 0.0, 1.0);

        // Tell the geometry engine about the discrete field mesh and how to interpret phases
        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometryVector;
        ge::Geometry_Engine tGeometryEngine(tInterpMesh, tGeometryEngineParameters);

        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        xtk::Model tXTKModel(tModelDimension,tInterpMesh, &tGeometryEngine);
        tXTKModel.mSameMesh = true;
        tXTKModel.mVerbose = false;

        // Do the cutting
        tXTKModel.decompose(tDecompositionMethods);

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = true;

        // add solution field to integration mesh
        std::string tIntegSolFieldName = "solution";
        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(tInterpMesh, tIntegMesh1);

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
        tPropConductivity->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );
        tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
        tPropConductivity->set_val_function( tPropValFunc );
        tPropConductivity->set_dof_derivative_functions( { tPropDerFunc } );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ 5.0 }} } );
        tPropDirichlet->set_val_function( tPropValConstFunc );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { {{ 20.0 }} } );
        tPropNeumann->set_val_function( tPropValConstFunc );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} ); // FIXME through the factory?
        tCMDiffLinIso->set_property( tPropConductivity, "Conductivity" );
        tCMDiffLinIso->set_local_properties();

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulk = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
        tIWGBulk->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );                          // FIXME through the factory?
        tIWGBulk->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER ); // FIXME through the factory?
        tIWGBulk->set_constitutive_model( tCMDiffLinIso, "Diffusion", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGDirichlet = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );                          // FIXME through the factory?
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER ); // FIXME through the factory?
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );                          // FIXME through the factory?
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER ); // FIXME through the factory?
        tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Master_Slave::MASTER );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "child_0" );
        tSetBulk1.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "parent_0" );
        tSetBulk2.set_IWGs( { tIWGBulk } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_set_name( "iside_g_0_p0_0_p1_1" );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_1" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk1;
        tSetInfo( 1 ) = tSetBulk2;
        tSetInfo( 2 ) = tSetNeumann;
        tSetInfo( 3 ) = tSetDirichlet;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                               tBSplineMeshIndex,
                                               tSetInfo );
        //------------------------------------------------------------------------------

        // clean up
        delete tModel;
        delete tInterpMesh;
    }

}/* END_TEST_CASE */

/* This UT tests:
* - the unpacking of the FEM inputs
* - the building of the FEM Model from these inputs
*/
//TEST_CASE("FEM_MDL_Input","[FEM_MDL_Input]")
//{
//    if(par_size() == 1)
//    {
//        std::string tFieldName = "Circle";
//
//        moris::uint tLagrangeMeshIndex = 0;
//        moris::uint tBSplineMeshIndex = 0;
//
//        moris::hmr::Parameters tParameters;
//
//        tParameters.set_number_of_elements_per_dimension( { {1}, {1}, {4} } );
//        tParameters.set_domain_dimensions({ {1}, {1}, {2} });
//        tParameters.set_domain_offset({ {0.0}, {0.0}, {0.0} });
//        tParameters.set_bspline_truncation( true );
//        tParameters.set_side_sets({ {5}, {6} });
//
//        tParameters.set_output_meshes( { {0} } );
//
//        tParameters.set_lagrange_orders  ( { {1} });
//        tParameters.set_lagrange_patterns({ {0} });
//
//        tParameters.set_bspline_orders   ( { {1} } );
//        tParameters.set_bspline_patterns ( { {0} } );
//
//        tParameters.set_union_pattern( 2 );
//        tParameters.set_working_pattern( 3 );
//
//        tParameters.set_refinement_buffer( 2 );
//        tParameters.set_staircase_buffer( 2);
//
//        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
//        tLagrangeToBSplineMesh( 0 ) = { {0} };
//
//        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );
//
//        hmr::HMR tHMR( tParameters );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//        // create field
//        std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( tFieldName, tLagrangeMeshIndex );
//
//        tField->evaluate_scalar_function( LevelSetPlaneFunction );
//
//        for( uint k=0; k<0; ++k )
//        {
//            tHMR.flag_surface_elements_on_working_pattern( tField );
//            tHMR.perform_refinement_based_on_working_pattern( 0 );
//
//            tField->evaluate_scalar_function( LevelSetPlaneFunction );
//        }
//
//        tHMR.finalize();
//
//        tHMR.save_to_exodus( 0, "./mdl_exo/mdl_input.e" );
//
//        std::shared_ptr< hmr::Interpolation_Mesh_HMR > tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex  );
//
//        moris::ge::Geometry_Field_HMR tFieldAsGeom(tField);
//
//        moris::Cell<ge::GEN_Geometry*> tGeometryVector = {&tFieldAsGeom};
//
//        // Tell the geometry engine about the discrete field mesh and how to interpret phases
//        ge::GEN_Phase_Table tPhaseTable (1);
//        ge::Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);
//
//        // Tell the XTK model that it should decompose with a C_HIERARCHY_TET4, on the same mesh that the level set field is defined on.
//        size_t tModelDimension = 3;
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//
//        xtk::Model tXTKModel(tModelDimension,tIPMesh.get(), &tGeometryEngine);
//        tXTKModel.mSameMesh = true;
//        tXTKModel.mVerbose = false;
//
//        // do the cutting
//        tXTKModel.decompose(tDecompositionMethods);
//
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = true;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldName = "solution";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
//
//        moris::mtk::Integration_Mesh* tIGMesh = tXTKModel.get_output_mesh(tOutputOptions);
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair( tIPMesh.get(), tIGMesh );
//
//        //------------------------------------------------------------------------------
//        // create a cell of cell of parameter list for fem
//        moris::Cell< moris::Cell< ParameterList > > tParameterList( 5 );
//
//        //------------------------------------------------------------------------------
//        // fill the property part of the parameter list
//        uint tNumProperties = 3;
//        tParameterList( 0 ).resize( tNumProperties );
//
//        // create parameter list for property 1
//        tParameterList( 0 )( 0 ) = prm::create_property_parameter_list();
//        tParameterList( 0 )( 0 ).set( "property_name",       "PropertyConductivity" );
//        tParameterList( 0 )( 0 ).set( "dof_dependencies",    "TEMP" );
//        tParameterList( 0 )( 0 ).set( "function_parameters", "1.0" );
//        tParameterList( 0 )( 0 ).set( "value_function",      "Func1" );
//
//        // create parameter list for property 2
//        tParameterList( 0 )( 1 ) = prm::create_property_parameter_list();
//        tParameterList( 0 )( 1 ).set( "property_name",       "PropertyDirichlet" );
//        tParameterList( 0 )( 1 ).set( "dof_dependencies",    "TEMP" );
//        tParameterList( 0 )( 1 ).set( "function_parameters", "5.0" );
//        tParameterList( 0 )( 1 ).set( "value_function",      "Func1" );
//
//        // create parameter list for property 2
//        tParameterList( 0 )( 2 ) = prm::create_property_parameter_list();
//        tParameterList( 0 )( 2 ).set( "property_name",       "PropertyNeumann" );
//        tParameterList( 0 )( 2 ).set( "dof_dependencies",    "TEMP" );
//        tParameterList( 0 )( 2 ).set( "function_parameters", "20.0" );
//        tParameterList( 0 )( 2 ).set( "value_function",      "Func1" );
//
//        //------------------------------------------------------------------------------
//        // fill the constitutive model part of the parameter list
//        uint tNumCMs = 1;
//        tParameterList( 1 ).resize( tNumCMs );
//
//        // create parameter list for constitutive model 1
//        tParameterList( 1 )( 0 ) = prm::create_constitutive_model_parameter_list();
//        tParameterList( 1 )( 0 ).set( "constitutive_name", "CMDiffusion" );
//        tParameterList( 1 )( 0 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
//        tParameterList( 1 )( 0 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP"), std::string("Temperature" ) );
//        tParameterList( 1 )( 0 ).set( "properties",        "PropertyConductivity,Conductivity" );
//
//        //------------------------------------------------------------------------------
//        // fill the stabilization parameter part of the parameter list
//        uint tNumSPs = 1;
//        tParameterList( 2 ).resize( tNumSPs );
//
//        // create parameter list for stabilization parameter 1
//        tParameterList( 2 )( 0 ) = prm::create_stabilization_parameter_parameter_list();
//        tParameterList( 2 )( 0 ).set( "stabilization_name",  "SPDirichlet" );
//        tParameterList( 2 )( 0 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
//        tParameterList( 2 )( 0 ).set( "function_parameters", "1.0" );
//        tParameterList( 2 )( 0 ).set( "master_properties",   "PropertyConductivity,Material" );
//
//        //------------------------------------------------------------------------------
//        // fill the IWG part of the parameter list
//        uint tNumIWGs = 3;
//        tParameterList( 3 ).resize( tNumIWGs );
//
//        // create parameter list for IWG 1
//        tParameterList( 3 )( 0 ) = prm::create_IWG_parameter_list();
//        tParameterList( 3 )( 0 ).set( "IWG_name",                   "SPATIALDIFF_BULK" );
//        tParameterList( 3 )( 0 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
//        tParameterList( 3 )( 0 ).set( "dof_residual",               "TEMP" );
//        tParameterList( 3 )( 0 ).set( "master_dof_dependencies",    "TEMP" );
//        tParameterList( 3 )( 0 ).set( "master_constitutive_models", "CMDiffusion,DiffLinIso" );
//        tParameterList( 3 )( 0 ).set( "mesh_set_names",             "child_0,parent_0" );
//
//        // create parameter list for IWG 2
//        tParameterList( 3 )( 1 ) = prm::create_IWG_parameter_list();
//        tParameterList( 3 )( 1 ).set( "IWG_name",                   "SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE" );
//        tParameterList( 3 )( 1 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) );
//        tParameterList( 3 )( 1 ).set( "dof_residual",               "TEMP" );
//        tParameterList( 3 )( 1 ).set( "master_dof_dependencies",    "TEMP" );
//        tParameterList( 3 )( 1 ).set( "master_properties",          "PropertyDirichlet,Dirichlet" );
//        tParameterList( 3 )( 1 ).set( "master_constitutive_models", "CMDiffusion,DiffLinIso" );
//        tParameterList( 3 )( 1 ).set( "stabilization_parameters",   "SPDirichlet,DirichletNitsche" );
//        tParameterList( 3 )( 1 ).set( "mesh_set_names",             "SideSet_1" );
//
//        // create parameter list for IWG 3
//        tParameterList( 3 )( 2 ) = prm::create_IWG_parameter_list();
//        tParameterList( 3 )( 2 ).set( "IWG_name",                "SPATIALDIFF_NEUMANN" );
//        tParameterList( 3 )( 2 ).set( "IWG_type",                static_cast< uint >( fem::IWG_Type::SPATIALDIFF_NEUMANN ) );
//        tParameterList( 3 )( 2 ).set( "dof_residual",            "TEMP" );
//        tParameterList( 3 )( 2 ).set( "master_dof_dependencies", "TEMP" );
//        tParameterList( 3 )( 2 ).set( "master_properties",       "PropertyNeumann,Neumann" );
//        tParameterList( 3 )( 2 ).set( "mesh_set_names",          "iside_g_0_p0_0_p1_1" );
//
//        //------------------------------------------------------------------------------
//        // path for property function reading
//        std::string tMeshFilePath = moris::get_base_moris_dir();
//        tMeshFilePath = tMeshFilePath + "projects/FEM/INT/test/data/FEM_input_test.so";
//        std::shared_ptr< Library_IO > tLibrary = std::make_shared< Library_IO >( tMeshFilePath );
//
//        // create model
//        fem::FEM_Model * tFEMModel = new fem::FEM_Model( &tMeshManager,
//                                                          tBSplineMeshIndex,
//                                                          tParameterList,
//                                                          tLibrary );
//
//        //------------------------------------------------------------------------------
//        // clean up
//        delete tFEMModel;
//    }
//
//}/* END_TEST_CASE */

}/* END_MORIS_NAMESPACE */




