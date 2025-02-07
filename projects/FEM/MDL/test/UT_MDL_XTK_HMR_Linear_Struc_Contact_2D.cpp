/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_Linear_Struc_Contact_2D.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "moris_typedefs.hpp"

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
#include "cl_FEM_CM_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_SP_Factory.hpp"              //FEM/INT/src
#include "cl_FEM_Set_User_Info.hpp"              //FEM/INT/src

#include "cl_MDL_Model.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
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

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Line.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

#include <functional>

namespace moris
{

//-------------------------------------------------------------------------------------
// Functions for Parameters in FEM
void ConstFunctionVal
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tMValFunctionContact
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                   0.0 },
                   { 0.0,                   aParameters( 0 )( 1 ) }};
}

//-------------------------------------------------------------------------------------

TEST_CASE("2D Linear Stuct Contract","[XTK_HMR_LS_Contact_2D]")
{
    if(par_size()<=1)
    {
        // Geometry Parameters
        moris::real tBlockL   =  1.0; /* Length of the block  (m) */
        moris::real tBlockH   =  1.0; /* Length of the block  (m) */
        moris::real tDomainLX = 2.0; /* Length of full domain in x  (m) */
        moris::real tDomainLY = 2.0; /* Length of full domain in y  (m) */
        moris::real tAngle =   0;
        Matrix<DDRMat> tCenterPoint = {{0.01111,0.0111}}; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/

        //Material Parameters
        moris::real tEa  = 1e6; // Pa
        moris::real tNua = 0.3;
        moris::real tEb  = 1e12; // Pa
        moris::real tNub = 0.3;

        // Boundary Conditions
        moris::real tDirchletX = 0.0; // m
        moris::real tDBCGamma  = 1000.0;

        moris::real tForce =  -10000.0; // N/m (normal to top surface)

        // Mesh Setup
        moris::uint tNumX   = 20; /* Number of elements in x*/
        moris::uint tNumY   = 20; /* Number of elements in y*/
        moris::uint tNumRef = 1;  /* Number of HMR refinements */

        // Files
        std::string tHMRIPMeshFileName = "./mdl_exo/mdl_xtk_hmr_2d.e";
        std::string tEnrIgMeshFileName = "./mdl_exo/contact_enr_ig.e";

        // flags
        bool tVizGhost = false;
        bool tVerboseGeometry = false;
        bool tVizIGMeshBeforeFEM = false;

        // Construct Left Plane

        Matrix<moris::DDRMat> tLeftCenters = {{ (tCenterPoint(0) + 0.01)  - ( tBlockL / 2.0 ) * std::cos(tAngle) ,
                                                (tCenterPoint(1) + 0.01) - tBlockH / 2.0 * std::sin(tAngle) }};
        Matrix<moris::DDRMat> tLeftNormal  = {{-std::cos(tAngle),-std::sin(tAngle)}};
        std::shared_ptr<moris::gen::Line> tLeftPlane = std::make_shared<moris::gen::Line>(tLeftCenters(0), tLeftCenters(1), tLeftNormal(0), tLeftNormal(1));
        auto tLeftPlaneFP = [&tLeftPlane] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tLeftPlane->get_field_value(aCoordinates); }; /*Lambda pointer for class */

        // Construct Right Plane
        Matrix<moris::DDRMat> tRightCenters = {{ (tCenterPoint(0) + 0.01)  + ( tBlockL / 2.0 ) * std::cos(tAngle) ,
                                                 (tCenterPoint(1) + 0.01) + tBlockH / 2.0 * std::sin(tAngle) }};
        Matrix<moris::DDRMat> tRightNormal  = {{std::cos(tAngle),std::sin(tAngle)}};
        std::shared_ptr<moris::gen::Line> tRightPlane = std::make_shared<moris::gen::Line>(tRightCenters(0), tRightCenters(1), tRightNormal(0), tRightNormal(1));
        auto tRightPlaneFP = [&tRightPlane] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tRightPlane->get_field_value(aCoordinates); }; /*Lambda pointer for class */

        // Construct Mid Plane
//        Matrix<moris::DDRMat> tMidCenters = {{ tCenterPoint(0) , tCenterPoint(1)  }};
//        Matrix<moris::DDRMat> tMidNormal  = {{std::cos(tAngle) - 1,1-std::sin(tAngle)}};
//        moris::gen::Plane tMidPlane(tMidCenters(0), tMidCenters(1), tMidNormal(0), tMidNormal(1));
//        auto tMidPlaneFP = [&tMidPlane] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tMidPlane.get_field_value_with_single_coordinate(aCoordinates); }; /*Lambda for class */

        // Construct Top Plane
        Matrix<moris::DDRMat> tTopCenters = {{(tCenterPoint(0) + 0.01) - tBlockH / 2.0 * std::sin(tAngle) ,
                                              (tCenterPoint(1) + 0.01) + tBlockH / 2.0 * std::cos(tAngle)}};
        Matrix<moris::DDRMat> tTopNormal  = {{std::cos(tAngle) - 1,1-std::sin(tAngle)}};
        std::shared_ptr<moris::gen::Line> tTopPlane = std::make_shared<moris::gen::Line>(tTopCenters(0), tTopCenters(1), tTopNormal(0), tTopNormal(1));
        auto tTopPlaneFP = [&tTopPlane] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tTopPlane->get_field_value(aCoordinates); }; /*Lambda for class */

        // Construct Bottom Plane
        Matrix<moris::DDRMat> tBottomCenters = {{(tCenterPoint(0) + 0.01) + tBlockH / 2.0 * std::sin(tAngle) ,
                                                 (tCenterPoint(1) + 0.01) - tBlockH / 2.0 * std::cos(tAngle)}};

        Matrix<moris::DDRMat> tBottomNormal  = {{1 - std::cos(tAngle),std::sin(tAngle) - 1}};
        std::shared_ptr<moris::gen::Line> tBottomPlane = std::make_shared<moris::gen::Line>(tBottomCenters(0), tBottomCenters(1), tBottomNormal(0), tBottomNormal(1));
        auto tBottomPlaneFP = [&tBottomPlane] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tBottomPlane->get_field_value(aCoordinates); }; /*Lambda for class */

        if(tVerboseGeometry)
        {
            std::cout << "\nLeft Plane :" << '\n';
            std::cout << "    xc = " << std::setw( 16 ) << tLeftCenters( 0 ) << std::setw( 16 ) << tLeftCenters( 1 ) << '\n';
            std::cout << "    n  = " << std::setw( 16 ) << tLeftNormal( 0 ) << std::setw( 16 ) << tLeftNormal( 1 ) << '\n';

            std::cout << "\nRight Plane :" << '\n';
            std::cout << "    xc = " << std::setw( 16 ) << tRightCenters( 0 ) << std::setw( 16 ) << tRightCenters( 1 ) << '\n';
            std::cout << "    n  = " << std::setw( 16 ) << tRightNormal( 0 ) << std::setw( 16 ) << tRightNormal( 1 ) << '\n';

            std::cout << "\nBottom Plane :" << '\n';
            std::cout << "    xc = " << std::setw( 16 ) << tBottomCenters( 0 ) << std::setw( 16 ) << tBottomCenters( 1 ) << '\n';
            std::cout << "    n  = " << std::setw( 16 ) << tBottomNormal( 0 ) << std::setw( 16 ) << tBottomNormal( 1 ) << '\n';

            std::cout << "\nTop Plane :" << '\n';
            std::cout << "    xc = " << std::setw( 16 ) << tTopCenters( 0 ) << std::setw( 16 ) << tTopCenters( 1 ) << '\n';
            std::cout << "    n  = " << std::setw( 16 ) << tTopNormal( 0 ) << std::setw( 16 ) << tTopNormal( 1 ) << '\n';
        }

        uint tLagrangeMeshIndex = 0;
        std::string tLeftFieldName   = "LeftPlane";
        std::string tRightFieldName  = "RightPlane";
        std::string tMidFieldName    = "MidPlane";
        std::string tTopFieldName    = "TopPlane";
        std::string tBottomFieldName = "BottomPlane";

        Parameter_List tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", tNumX, tNumY );
        tParameters.set( "domain_dimensions", tDomainLX, tDomainLY );
        tParameters.set( "domain_offset", -tDomainLX / 2, -tDomainLY / 2 );
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
        tParameters.set("use_number_aura",1);

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement();

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //  create field
        std::shared_ptr< moris::hmr::Field > tLeftField   = tMesh->create_field( tLeftFieldName, tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tRightField  = tMesh->create_field( tRightFieldName, tLagrangeMeshIndex );
//        std::shared_ptr< moris::hmr::Field > tMidField    = tMesh->create_field( tMidFieldName, tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tTopField    = tMesh->create_field( tTopFieldName, tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tBottomField = tMesh->create_field( tBottomFieldName, tLagrangeMeshIndex );

        for( uint k=0; k<tNumRef; ++k )
        {
            tLeftField->evaluate_scalar_function( tLeftPlaneFP );
            tRightField->evaluate_scalar_function( tRightPlaneFP );
//            tMidField->evaluate_scalar_function( tMidPlaneFP );
            tTopField->evaluate_scalar_function( tTopPlaneFP );
            tBottomField->evaluate_scalar_function( tBottomPlaneFP );

            tHMR.flag_surface_elements_on_working_pattern( tLeftField );
            tHMR.flag_surface_elements_on_working_pattern( tRightField );
//            tHMR.flag_surface_elements_on_working_pattern( tMidField );
            tHMR.flag_surface_elements_on_working_pattern( tTopField );
            tHMR.flag_surface_elements_on_working_pattern( tBottomField );

            tHMR.perform_refinement_based_on_working_pattern( 0 );

        }

        tLeftField->evaluate_scalar_function( tLeftPlaneFP );
        tRightField->evaluate_scalar_function( tRightPlaneFP );
        tTopField->evaluate_scalar_function( tTopPlaneFP );
        tBottomField->evaluate_scalar_function( tBottomPlaneFP );

        tHMR.finalize();

        tHMR.save_to_exodus( 0, tHMRIPMeshFileName );

        moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        //-----------------------------------------------------------------------------------------------
        Vector< std::shared_ptr< moris::gen::Geometry > > tGeometryVector( 4 );
        tGeometryVector( 0 ) = std::make_shared< gen::Level_Set_Geometry >( tLeftPlane );
        tGeometryVector( 1 ) = std::make_shared< gen::Level_Set_Geometry >( tRightPlane );
        tGeometryVector( 2 ) = std::make_shared< gen::Level_Set_Geometry >( tTopPlane );
        tGeometryVector( 3 ) = std::make_shared< gen::Level_Set_Geometry >( tBottomPlane );

        size_t tModelDimension = 2;
        moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometryVector;
        moris::gen::Geometry_Engine tGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);
        xtk::Model tXTKModel(tModelDimension,tInterpolationMesh,&tGeometryEngine);
        tXTKModel.mVerbose = false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
//        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        // get meshes for FEM
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        if(tVizGhost)
        {
            // access ghost
            xtk::Ghost_Stabilization & tGhostStab = tXTKModel.get_ghost_stabilization();

            for(moris_index iG = 0; iG < 16; iG++)
            {
                tGhostStab.visualize_ghost_on_mesh(iG);
            }
        }

        if(tVizIGMeshBeforeFEM)
        {
            tEnrIntegMesh.deactivate_empty_sets();
            // Write mesh
            moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
            writer.write_mesh("", tEnrIgMeshFileName, "", "temp.exo");

            // Write the fields
            writer.set_time(0.0);
            writer.close_file();

        }

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropEModA = std::make_shared< fem::Property >();
        tPropEModA->set_parameters( { {{ tEa }} } );
        tPropEModA->set_val_function( ConstFunctionVal );

        std::shared_ptr< fem::Property > tPropEModB = std::make_shared< fem::Property >();
        tPropEModB->set_parameters( { {{ tEb }} } );
        tPropEModB->set_val_function( ConstFunctionVal );

        std::shared_ptr< fem::Property > tPropNua = std::make_shared< fem::Property >();
        tPropNua->set_parameters( { {{ tNua }} } );
        tPropNua->set_val_function( ConstFunctionVal );

        std::shared_ptr< fem::Property > tPropNub = std::make_shared< fem::Property >();
        tPropNub->set_parameters( { {{ tNub }} } );
        tPropNub->set_val_function( ConstFunctionVal );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { {{ tDirchletX }, { tDirchletX }} } );
        tPropDirichlet->set_val_function( ConstFunctionVal );

        std::shared_ptr< fem::Property > tPropDirichlet2 = std::make_shared< fem::Property >();
        tPropDirichlet2->set_parameters( { {{ 1.0, 1.0 }} } );
        tPropDirichlet2->set_val_function( tMValFunctionContact );

        std::shared_ptr< fem::Property > tPropTraction = std::make_shared< fem::Property >();
        tPropTraction->set_parameters( {{{ tForce*std::sin(tAngle) } , { tForce*std::cos(tAngle) }}} );
        tPropTraction->set_val_function( ConstFunctionVal );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 =
                tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tCMStrucLinIso1->set_property( tPropEModA, "YoungsModulus" );
        tCMStrucLinIso1->set_property( tPropNua, "PoissonRatio" );
        tCMStrucLinIso1->set_model_type(fem::Model_Type::PLANE_STRESS);
        tCMStrucLinIso1->set_space_dim( 2 );
        tCMStrucLinIso1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso2 =
                tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
        tCMStrucLinIso2->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tCMStrucLinIso2->set_property( tPropEModB, "YoungsModulus" );
        tCMStrucLinIso2->set_property( tPropNub, "PoissonRatio" );
        tCMStrucLinIso2->set_model_type(fem::Model_Type::PLANE_STRESS);
        tCMStrucLinIso2->set_space_dim( 2 );
        tCMStrucLinIso2->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { {{ tDBCGamma }} } );
        tSPDirichletNitsche->set_property( tPropEModB, "Material", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
        tSPNitscheInterface->set_property( tPropEModB, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface->set_property( tPropEModA, "Material", mtk::Leader_Follower::FOLLOWER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGBulkA =
                tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGBulkA->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGBulkA->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGBulkA->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGBulkB =
                tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
        tIWGBulkB->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGBulkB->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGBulkB->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGDirichlet =
                tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWGDirichlet->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGDirichlet->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
        tIWGDirichlet->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
        tIWGDirichlet->set_property( tPropDirichlet2, "Select", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGNeumann =
                tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGNeumann->set_property( tPropTraction, "Traction", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::IWG > tIWGInterface =
                tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
        tIWGInterface->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }},mtk::Leader_Follower::FOLLOWER );
        tIWGInterface->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
        tIWGInterface->set_constitutive_model( tCMStrucLinIso2, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIWGInterface->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );

        // define set info
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p0" );
        tSetBulk1.set_IWGs( { tIWGBulkA } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p0" );
        tSetBulk2.set_IWGs( { tIWGBulkA } );

        fem::Set_User_Info tSetBulk3;
        tSetBulk3.set_mesh_set_name( "HMR_dummy_c_p1" );
        tSetBulk3.set_IWGs( { tIWGBulkB } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk4.set_IWGs( { tIWGBulkB } );

        fem::Set_User_Info tSetBulk5;
        tSetBulk5.set_mesh_set_name( "HMR_dummy_c_p5" );
        tSetBulk5.set_IWGs( { tIWGBulkB } );

        fem::Set_User_Info tSetBulk6;
        tSetBulk6.set_mesh_set_name( "HMR_dummy_n_p5" );
        tSetBulk6.set_IWGs( { tIWGBulkB } );

        fem::Set_User_Info tSetBulk7;
        tSetBulk7.set_mesh_set_name( "HMR_dummy_c_p9" );
        tSetBulk7.set_IWGs( { tIWGBulkB } );

        fem::Set_User_Info tSetBulk8;
        tSetBulk8.set_mesh_set_name( "HMR_dummy_n_p9" );
        tSetBulk8.set_IWGs( { tIWGBulkB } );

        fem::Set_User_Info tSetDirichlet1;
        tSetDirichlet1.set_mesh_set_name( "SideSet_1_n_p9" );
        tSetDirichlet1.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetDirichlet2;
        tSetDirichlet2.set_mesh_set_name( "SideSet_1_c_p9" );
        tSetDirichlet2.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetDirichlet3;
        tSetDirichlet3.set_mesh_set_name( "SideSet_1_c_p1" );
        tSetDirichlet3.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetDirichlet4;
        tSetDirichlet4.set_mesh_set_name( "SideSet_1_n_p1" );
        tSetDirichlet4.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetDirichlet5;
        tSetDirichlet5.set_mesh_set_name( "SideSet_1_c_p5" );
        tSetDirichlet5.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetDirichlet6;
        tSetDirichlet6.set_mesh_set_name( "SideSet_1_n_p5" );
        tSetDirichlet6.set_IWGs( { tIWGDirichlet } );

        std::cout << "tEnrIntegMesh.get_interface_side_set_name(3,9,8) = " << tEnrIntegMesh.get_interface_side_set_name( 3, 9, 8 ) << " Index = " << tEnrIntegMesh.get_set_index_by_name( tEnrIntegMesh.get_interface_side_set_name( 3, 9, 8 ) ) << '\n';
        std::cout << "tEnrIntegMesh.get_interface_side_set_name(3,5,4) = " << tEnrIntegMesh.get_interface_side_set_name( 3, 5, 4 ) << tEnrIntegMesh.get_set_index_by_name( tEnrIntegMesh.get_interface_side_set_name( 3, 5, 4 ) ) << '\n';

        fem::Set_User_Info tSetNeumann1;
        tSetNeumann1.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name(2,0,2) );
        tSetNeumann1.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetNeumann2;
        tSetNeumann2.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name(3,5,4) );
        tSetNeumann2.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetNeumann3;
        tSetNeumann3.set_mesh_set_name( tEnrIntegMesh.get_interface_side_set_name(3,9,8) );
        tSetNeumann3.set_IWGs( { tIWGNeumann } );

        /* This is the contact interface*/
        fem::Set_User_Info tSetInterface1;
        tSetInterface1.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name(0,1) );
        tSetInterface1.set_IWGs( { tIWGInterface } );

        fem::Set_User_Info tSetInterface2;
        tSetInterface2.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name(1,9) );
        tSetInterface2.set_IWGs( { tIWGInterface } );

        fem::Set_User_Info tSetInterface3;
        tSetInterface3.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name(1,5) );
        tSetInterface3.set_IWGs( { tIWGInterface } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 20 );
        tSetInfo( 0 )  = tSetBulk1;
        tSetInfo( 1 )  = tSetBulk2;
        tSetInfo( 2 )  = tSetBulk3;
        tSetInfo( 3 )  = tSetBulk4;
        tSetInfo( 4 )  = tSetBulk5;
        tSetInfo( 5 )  = tSetBulk6;
        tSetInfo( 6 )  = tSetBulk7;
        tSetInfo( 7 )  = tSetBulk8;
        tSetInfo( 8 )  = tSetDirichlet1;
        tSetInfo( 9 )  = tSetDirichlet2;
        tSetInfo( 10 )  = tSetDirichlet3;
        tSetInfo( 11 )  = tSetDirichlet4;
        tSetInfo( 12 )  = tSetDirichlet5;
        tSetInfo( 13 )  = tSetDirichlet6;
        tSetInfo( 14 ) = tSetNeumann1;
        tSetInfo( 15 ) = tSetNeumann2;
        tSetInfo( 16 ) = tSetNeumann3;
        tSetInfo( 17 ) = tSetInterface1;
        tSetInfo( 18 ) = tSetInterface2;
        tSetInfo( 19 ) = tSetInterface3;

        // create model
        mdl::Model * tModel = new mdl::Model( tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Vector< enum MSI::Dof_Type > tDofTypesU( 2 );
        tDofTypesU( 0 ) = MSI::Dof_Type::UX;
        tDofTypesU( 1 ) = MSI::Dof_Type::UY;

        dla::Solver_Factory  tSolFactory;
        Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        tLinearSolverParameterList.set( "AZ_max_iter", 10000 );
        tLinearSolverParameterList.set( "AZ_solver", AZ_gmres );
        tLinearSolverParameterList.set( "AZ_subdomain_solve", AZ_ilu );
        tLinearSolverParameterList.set( "AZ_graph_fill", 10 );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;
        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 3 );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
        //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolverMain;
        tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

        // Create solver database
        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

        tNonlinearSolverMain       .set_solver_warehouse( &tSolverWarehouse );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 3: create time Solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolverMain );

        tsa::Time_Solver tTimeSolver;
        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tTimeSolver.set_dof_type_list( tDofTypesU );

        //------------------------------------------------------------------------------
        tTimeSolver.solve();

        // output solution and meshes

//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = false;
//
//        // add solution field to integration mesh
//        std::string tIntegSolFieldNameUX = "UX";
//        std::string tIntegSolFieldNameUY = "UY";
//        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldNameUX, tIntegSolFieldNameUY};
//
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//
//        // Write to Integration mesh for visualization
//        Matrix<DDRMat> tIntegSolUX = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UX );
//        Matrix<DDRMat> tIntegSolUY = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::UY );
//
//        //    print(tIntegSolUX,"tIntegSolUX");
//        //    print(tIntegSolUY,"tIntegSolUY");
//
//        Matrix<DDRMat> tSTKIntegSolUX(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//        Matrix<DDRMat> tSTKIntegSolUY(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
//
//        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
//        {
//            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//            tSTKIntegSolUX(i) = tIntegSolUX(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//            tSTKIntegSolUY(i) = tIntegSolUY(tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE));
//        }
//
//        // add solution field to integration mesh
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUX,EntityRank::NODE,tSTKIntegSolUX);
//        tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds(tIntegSolFieldNameUY,EntityRank::NODE,tSTKIntegSolUY);
//
//        //    Matrix<DDRMat> tFullSol;
//        //    tNonlinearSolver.get_full_solution(tFullSol);
//        //
//        //    print(tFullSol,"tFullSol");
//

        delete tModel;
        delete tInterpolationMesh;
    }
}
}

