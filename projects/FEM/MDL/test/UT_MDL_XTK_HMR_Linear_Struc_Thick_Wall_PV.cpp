/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_Linear_Struc_Thick_Wall_PV.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
#include "cl_Geom_Field.hpp"
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
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "cl_MDL_Model.hpp"

#include "cl_VIS_Output_Manager.hpp"
#include "cl_VIS_Output_Enums.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Warehouse.hpp"

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

#include "HDF5_Tools.hpp"

#include "fn_norm.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Geom_Field_HMR.hpp"
#include "cl_GEN_Geometry.hpp"

#include "fn_PRM_HMR_Parameters.hpp"
#include "cl_Plane.hpp"

#include <functional>

namespace moris
{

//-------------------------------------------------------------------------------------
// Functions for Parameters in FEM
void ConstFunctionVal
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tMValFunctionContact2( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                            moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                            moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = {{ aParameters( 0 )( 0 ),                   0.0 },
                   { 0.0,                   aParameters( 0 )( 1 ) }};
}

/*!
 * Analytic function for UR
 *https://engineering.purdue.edu/~ce597m/Handouts/Theory%20of%20elasticity%20by%20Timoshenko%20and%20Goodier.pdf
 *  tNu, tE  tA, tB, tp
 */
moris::Matrix< moris::DDRMat >
AnalyticSigr( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
              moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
              moris::fem::Field_Interpolator_Manager         * aFIManager)
{
    // Points coordinate
    Matrix<DDRMat> tCoords = aFIManager->get_IG_geometry_interpolator()->valx();

    // compute R
    moris::real tR = std::pow( tCoords(0) * tCoords(0) + tCoords(1) * tCoords(1),0.5);

    moris::real tNu = aParameters(0)(0);
    moris::real tE  = aParameters(0)(1);
    moris::real tA = aParameters(0)(2);
    moris::real tB = aParameters(0)(3);
    moris::real tPi = aParameters(0)(4);
    moris::real tPo = 0;
    moris::real tC  = tB / tA;

    moris::real tASq = std::pow(tA,2);
    moris::real tBSq = std::pow(tB,2);
    moris::real tRSq = std::pow(tR,2);

    //return {{(tASq * tBSq * (tPo - tPi)) /((tBSq - tASq) * tRSq) + (tPi * tASq - tPo *tBSq)/ (tBSq - tASq) }} ;

return {{(tASq * tPi) / (tBSq - tASq) * ( 1- tBSq/tRSq)}};
}

moris::Matrix< moris::DDRMat >
AnalyticSigTh( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
               moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
               moris::fem::Field_Interpolator_Manager         * aFIManager)
{
    // Points coordinate
    Matrix<DDRMat> tCoords = aFIManager->get_IG_geometry_interpolator()->valx();

    // compute R
    moris::real tR = std::pow( tCoords(0) * tCoords(0) + tCoords(1) * tCoords(1),0.5);

    moris::real tNu = aParameters(0)(0);
    moris::real tE  = aParameters(0)(1);
    moris::real tA = aParameters(0)(2);
    moris::real tB = aParameters(0)(3);
    moris::real tPi = aParameters(0)(4);
    moris::real tPo = 0;
    moris::real tC  = tB / tA;

    moris::real tASq = std::pow(tA,2);
    moris::real tBSq = std::pow(tB,2);
    moris::real tRSq = std::pow(tR,2);

    //return {{- (tASq * tBSq *(tPo - tPi)) / ((tBSq - tASq) * tRSq) + (tPi*tASq - tPo * tBSq)/(tBSq - tASq) }};

     return {{(tASq * tPi) / (tBSq - tASq) * ( 1 + tBSq/tRSq)}};
}

void AnalyticUr( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
                 moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                 moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // Points coordinate
    Matrix<DDRMat> tCoords = aFIManager->get_IG_geometry_interpolator()->valx();

    // compute R
    moris::real tR = std::pow( tCoords(0) * tCoords(0) + tCoords(1) * tCoords(1),0.5);

    // compute R
    moris::real tNu = aParameters(0)(0);
    moris::real tE  = aParameters(0)(1);
    moris::real tA = aParameters(0)(2);
    moris::real tB = aParameters(0)(3);
    moris::real tPi = aParameters(0)(4);
    moris::real tPo = 0;
    moris::real tC  = tB / tA;

    moris::real tASq = std::pow(tA,2);
    moris::real tBSq = std::pow(tB,2);
    moris::real tRSq = std::pow(tR,2);

    moris::real tSigR = (tASq * tPi) / (tBSq - tASq) * ( 1 - tBSq/tRSq);
    moris::real tSigTh = (tASq * tPi) / (tBSq - tASq) * ( 1 + tBSq/tRSq);

    moris::real tUr = tR/tE * (tSigTh - tNu*tSigR);

    aPropMatrix = { { tUr } };
}

bool fOutputVizMesh( moris::tsa::Time_Solver * )
{
    return true;
}

//-------------------------------------------------------------------------------------

TEST_CASE("2D Linear Stuct Thick Walled Pressure Vessel","[XTK_HMR_LS_PV]")
{
    if(par_size()<=1)
    {
//        moris::uint tNumRefLevels  = 8;
//        moris::uint tNumOrders  = 2;
        moris::uint tAdaptiveOrUniform = 1; // 0 adapative 1 uniform

        // flags
        bool tVizGhost           = false;
        bool tGhostInModel       = true;
        bool tVerboseGeometry    = false;
        bool tVizIGMeshBeforeFEM = true;
//
//        for(moris::uint iOrder = 1; iOrder<tNumOrders+1; iOrder++)
//        {
//            for( moris::uint  iRef = 0; iRef<tNumRefLevels; iRef++)
//            {
                // Geometry Parameters
                moris::real tA = 8.251;  // Inner Radius
                moris::real tB = 40.01;  // Outer Radius
                Matrix<DDRMat> tCenterPoint = {{0.0,0.0}}; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
                moris::real tDomainLX = 50.0; /* Length of full domain in x  (m) */
                moris::real tDomainLY = 50.0; /* Length of full domain in y  (m) */

                //Material Parameters
                moris::real tE  = 1000; // Pa
                moris::real tNu = 0.3;

                // Boundary Condition Parameters
                moris::real tDirchletX = 0.0; // m
                moris::real tDirchletY = 0.0; // m
                moris::real tDBCGamma  = 10000.0;

                moris::real tP = -1; // N/m (pressure on interior of cylinder), negative for into body pressure

                // Mesh Setup
                moris::uint tNumX   = 4; /* Number of elements in x*/
                moris::uint tNumY   = 4; /* Number of elements in y*/
                moris::uint tNumRef = 0;  /* Number of HMR refinements */

                moris::uint iRef = 5;

                if(tAdaptiveOrUniform == 0)
                {
                    tNumRef = iRef;
                }
                if(tAdaptiveOrUniform == 1)
                {
                    tNumX = tNumX + iRef*tNumX;
                    tNumY = tNumY + iRef*tNumY;
                }

                moris::uint tOrder  = 1;  /* Lagrange Order and Bspline Order (forced to be same for this example) */

                // Files
                std::string tHMRIPMeshFileName = "./mdl_exo/thick_wall_pv_hmr_ip.e";
                std::string tEnrIgMeshFileName = "./mdl_exo/thick_wall_pv_ig.e";
//                std::string tPath = "/home/doble/Documents/Comps/exofiles_gf_" + std::to_string(tGhostInModel) + "_ar_"+ std::to_string(tAdaptiveOrUniform) + "/";
                std::string tPath = "";

                std::string tFileBase = "thick_wall_pv_ig_o_" + std::to_string(tOrder) + "_nr_" + std::to_string(iRef) ;
                std::string tVizMeshFileName =tFileBase + ".e";

                // field names
                std::string tInnerCircleFieldName = "InnerCircle";
                std::string tOuterCircleFieldName = "OuterCircle";

                // Dof Types
                Cell< MSI::Dof_Type > tDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY };

                // Construct inner circle Plane
                ge::Circle tInnerCircle(tA,tCenterPoint(0),tCenterPoint(1));
                auto tInnerCircleFP = [&tInnerCircle] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tInnerCircle.evaluate_field_value_with_single_coordinate(aCoordinates); }; /*Lambda pointer for class */

                ge::Circle tOuterCircle(tB,tCenterPoint(0),tCenterPoint(1));
                auto tOuterCircleFP = [&tOuterCircle] (moris::Matrix< moris::DDRMat > const & aCoordinates) { return tOuterCircle.evaluate_field_value_with_single_coordinate(aCoordinates); }; /*Lambda pointer for class */

                if(tVerboseGeometry)
                {
                    std::cout<<"\nInner Circle:"<<std::endl;
                    std::cout<<"    xc = "<<std::setw(16)<<tCenterPoint(0)<<std::setw(16)<<tCenterPoint(1)<<std::endl;
                    std::cout<<"    r  = "<<std::setw(16)<<tA<<std::endl;

                    std::cout<<"\nOuter Circle :"<<std::endl;
                    std::cout<<"    xc = "<<std::setw(16)<<tCenterPoint(0)<<std::setw(16)<<tCenterPoint(1)<<std::endl;
                    std::cout<<"    r  = "<<std::setw(16)<<tB<<std::endl;
                }

                uint tLagrangeMeshIndex = 0;

                moris::ParameterList tParameters = prm::create_hmr_parameter_list();

                tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
                tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
                tParameters.set( "domain_offset", "0,0" );
                tParameters.set( "domain_sidesets", "1,2,3,4" );
                tParameters.set( "lagrange_output_meshes", "0" );

                tParameters.set( "lagrange_orders", std::to_string(tOrder) );
                tParameters.set( "lagrange_pattern",std::string( "0" ));
                tParameters.set( "bspline_orders", std::to_string(tOrder));
                tParameters.set( "bspline_pattern",std::string( "0" ));

                tParameters.set( "lagrange_to_bspline",std::string( "0" ));

                tParameters.set( "truncate_bsplines", 1 );
                tParameters.set( "refinement_buffer", 3 );
                tParameters.set( "staircase_buffer", 3 );
                tParameters.set( "initial_refinement", "0" );
                tParameters.set( "initial_refinement_pattern", "0" );

                tParameters.set( "use_multigrid", 0 );
                tParameters.set( "severity_level", 2 );
                tParameters.set("use_number_aura", 1);

                hmr::HMR tHMR( tParameters );

                //initial refinement
                tHMR.perform_initial_refinement();

                std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                //  create field
                std::shared_ptr< moris::hmr::Field > tInnerCircleField = tMesh->create_field( tInnerCircleFieldName, tLagrangeMeshIndex );
                std::shared_ptr< moris::hmr::Field > tOuterCircleField = tMesh->create_field( tOuterCircleFieldName, tLagrangeMeshIndex );

                for( uint k=0; k<tNumRef; ++k )
                {
                    tInnerCircleField->evaluate_scalar_function( tInnerCircleFP );
                    tOuterCircleField->evaluate_scalar_function( tOuterCircleFP );

                    tHMR.flag_surface_elements_on_working_pattern( tInnerCircleField );
                    tHMR.flag_surface_elements_on_working_pattern( tOuterCircleField );

                    tHMR.perform_refinement_based_on_working_pattern( 0 );
                }

                tInnerCircleField->evaluate_scalar_function( tInnerCircleFP );
                tOuterCircleField->evaluate_scalar_function( tOuterCircleFP );

                tHMR.finalize();

//                std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh
//                = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh
                = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

                //-----------------------------------------------------------------------------------------------
                moris::ge::GEN_Geom_Field_HMR tInnerCircleForGE(tInnerCircleField);
                moris::ge::GEN_Geom_Field_HMR tOuterCircleForGE(tOuterCircleField);

                // NOTE the order of this geometry vector is important. If it changes the resulting bulk phase of the output mesh change.
                moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { & tInnerCircleForGE , & tOuterCircleForGE};

                size_t tModelDimension = 2;
                moris::ge::GEN_Geometry_Engine tGeometryEngine(tGeometryVector, tModelDimension);
                xtk::Model tXTKModel(tModelDimension,tInterpolationMesh,&tGeometryEngine);
                tXTKModel.mVerbose = false;

                //Specify decomposition Method and Cut Mesh ---------------------------------------
                Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
                tXTKModel.decompose(tDecompositionMethods);

                tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
                if(tGhostInModel)
                {
                    tXTKModel.construct_face_oriented_ghost_penalization_cells();
                }

                // get meshes for FEM
                xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
                xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

                if(tVizGhost)
                {
                    // access ghost
                    xtk::Ghost_Stabilization & tGhostStab = tXTKModel.get_ghost_stabilization();

                    for(moris_index iG = 0; iG < (moris_index)tPhaseTable.get_num_phases(); iG++)
                    {
                        tGhostStab.visualize_ghost_on_mesh(iG);
                    }
                }

//                if(tVizIGMeshBeforeFEM)
//                {
//                    // Write mesh
//                    Writer_Exodus writer(&tEnrIntegMesh);
//                    writer.write_mesh("", tEnrIgMeshFileName, "", "temp.exo");
//
//                    // Write the fields
//                    writer.set_time(0.0);
//                    writer.close_file();
//                }

                // place the pair in mesh manager
                std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
                tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

                //------------------------------------------------------------------------------
                // create the properties
                moris::Cell< MSI::Dof_Type > tResDofTypes = { MSI::Dof_Type::UX, MSI::Dof_Type::UY };

                std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
                tPropEMod->set_parameters( { {{ tE }} } );
                tPropEMod->set_val_function( ConstFunctionVal );

                std::shared_ptr< fem::Property > tPropNua = std::make_shared< fem::Property >();
                tPropNua->set_parameters( { {{ tNu }} } );
                tPropNua->set_val_function( ConstFunctionVal );

                std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
                tPropDirichlet->set_parameters( { {{ tDirchletX }, { tDirchletX }} } );
                tPropDirichlet->set_val_function( ConstFunctionVal );

                std::shared_ptr< fem::Property > tPropDirichletSelectX = std::make_shared< fem::Property >();
                tPropDirichletSelectX->set_parameters( { {{ 1.0, 0.0 }} } );
                tPropDirichletSelectX->set_val_function( tMValFunctionContact2 );

                std::shared_ptr< fem::Property > tPropDirichletSelectY = std::make_shared< fem::Property >();
                tPropDirichletSelectY->set_parameters( { {{ 0.0, 1.0 }} } );
                tPropDirichletSelectY->set_val_function( tMValFunctionContact2 );

                std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
                tPropNeumann->set_parameters( {{{ tP } , { tP }}} );
                tPropNeumann->set_val_function( ConstFunctionVal );

                std::shared_ptr< fem::Property > tPropUrAnalytic = std::make_shared< fem::Property >();
                tPropUrAnalytic->set_parameters( {{{ tNu } , { tE }, {tA}, {tB},{-tP}}} );
                tPropUrAnalytic->set_val_function( AnalyticUr );

                // define constitutive models
                fem::CM_Factory tCMFactory;

                std::shared_ptr< fem::Constitutive_Model > tCMStrucLinIso1 = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                tCMStrucLinIso1->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
                tCMStrucLinIso1->set_property( tPropEMod, "YoungsModulus" );
                tCMStrucLinIso1->set_property( tPropNua, "PoissonRatio" );
                tCMStrucLinIso1->set_model_type(fem::Model_Type::PLANE_STRESS);
                tCMStrucLinIso1->set_space_dim( 2 );
                tCMStrucLinIso1->set_local_properties();

                // define stabilization parameters
                fem::SP_Factory tSPFactory;
                std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
                tSPDirichletNitsche->set_parameters( { {{ tDBCGamma }} } );
                tSPDirichletNitsche->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::Stabilization_Parameter > tSPGhost = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
                tSPGhost->set_parameters( { {{ 0.01 }} });
                tSPGhost->set_property( tPropEMod, "Material", mtk::Leader_Follower::LEADER );

                // define the IQIs
                fem::IQI_Factory tIQIFactory;

                std::shared_ptr< fem::IQI > tIQIUrAnalytic = tIQIFactory.create_IQI( fem::IQI_Type::PROPERTY );
                tIQIUrAnalytic->set_property( tPropUrAnalytic, "Property", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IQI > tIQIUX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
                tIQIUX->set_quantity_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
                tIQIUX->set_output_type( vis::Output_Type::UX );
                tIQIUX->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
                tIQIUX->set_output_type_index( 0 );

                std::shared_ptr< fem::IQI > tIQIUY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
                tIQIUY->set_quantity_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
                tIQIUY->set_output_type( vis::Output_Type::UY );
                tIQIUY->set_dof_type_list( { tResDofTypes }, mtk::Leader_Follower::LEADER );
                tIQIUY->set_output_type_index( 1 );

                // define the IWGs
                fem::IWG_Factory tIWGFactory;

                std::shared_ptr< fem::IWG > tIWGBulkA = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_BULK );
                tIWGBulkA->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGBulkA->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
                tIWGBulkA->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGDirichletX = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
                tIWGDirichletX->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGDirichletX->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
                tIWGDirichletX->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
                tIWGDirichletX->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
                tIWGDirichletX->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
                tIWGDirichletX->set_property( tPropDirichletSelectX, "Select", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGDirichletY = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
                tIWGDirichletY->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGDirichletY->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
                tIWGDirichletY->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
                tIWGDirichletY->set_constitutive_model( tCMStrucLinIso1, "ElastLinIso", mtk::Leader_Follower::LEADER );
                tIWGDirichletY->set_property( tPropDirichlet, "Dirichlet", mtk::Leader_Follower::LEADER );
                tIWGDirichletY->set_property( tPropDirichletSelectY, "Select", mtk::Leader_Follower::LEADER );

                std::shared_ptr< fem::IWG > tIWGNeumann = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_NEUMANN );
                tIWGNeumann->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
                tIWGNeumann->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
                tIWGNeumann->set_property( tPropNeumann, "Neumann", mtk::Leader_Follower::LEADER );

                // GHOST STABILIZATION SETUP
                std::shared_ptr< fem::IWG > tIWGGhost = tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
                tIWGGhost->set_residual_dof_type( tDofTypes );
                tIWGGhost->set_dof_type_list( { tDofTypes } );
                tIWGGhost->set_dof_type_list( { tDofTypes }, mtk::Leader_Follower::FOLLOWER );
                tIWGGhost->set_stabilization_parameter( tSPGhost, "GhostSP" );

                // define set info
                fem::Set_User_Info tSetBulk1;
                tSetBulk1.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_c_p2") );
                tSetBulk1.set_IWGs( { tIWGBulkA } );
                tSetBulk1.set_IQIs( { tIQIUX, tIQIUY,tIQIUrAnalytic } );

                fem::Set_User_Info tSetBulk2;
                tSetBulk2.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("HMR_dummy_n_p2") );
                tSetBulk2.set_IWGs( { tIWGBulkA } );
                tSetBulk2.set_IQIs( { tIQIUX, tIQIUY,tIQIUrAnalytic } );

                fem::Set_User_Info tSetDirichletXc;
                tSetDirichletXc.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_4_c_p2") );
                tSetDirichletXc.set_IWGs( { tIWGDirichletX } );

                fem::Set_User_Info tSetDirichletXn;
                tSetDirichletXn.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_4_n_p2") );
                tSetDirichletXn.set_IWGs( { tIWGDirichletX } );

                fem::Set_User_Info tSetDirichletYc;
                tSetDirichletYc.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_c_p2") );
                tSetDirichletYc.set_IWGs( { tIWGDirichletY } );

                fem::Set_User_Info tSetDirichletYn;
                tSetDirichletYn.set_mesh_index( tEnrIntegMesh.get_set_index_by_name("SideSet_1_n_p2") );
                tSetDirichletYn.set_IWGs( { tIWGDirichletY } );

                fem::Set_User_Info tSetNeumann1;
                tSetNeumann1.set_mesh_index( tEnrIntegMesh.get_set_index_by_name(tEnrIntegMesh.get_interface_side_set_name(0,2,0)) );
                tSetNeumann1.set_IWGs( { tIWGNeumann } );

                // create a cell of set info
                moris::Cell< fem::Set_User_Info > tSetInfo( 7 );
                tSetInfo( 0 )  = tSetBulk1;
                tSetInfo( 1 )  = tSetBulk2;
                tSetInfo( 2 )  = tSetDirichletXc;
                tSetInfo( 3 )  = tSetDirichletYc;
                tSetInfo( 4 )  = tSetDirichletXn;
                tSetInfo( 5 )  = tSetDirichletYn;
                tSetInfo( 6 )  = tSetNeumann1;

                if(tGhostInModel)
                {
                    xtk::Ghost_Stabilization & tGhostStab = tXTKModel.get_ghost_stabilization();

                    fem::Set_User_Info tSetGhost;
                    tSetGhost.set_mesh_index(  tEnrIntegMesh.get_set_index_by_name(tGhostStab.get_ghost_dbl_side_set_name(2)));
                    tSetGhost.set_IWGs({tIWGGhost});
                    tSetInfo.push_back(tSetGhost);
                }

                // create model
                mdl::Model * tModel = new mdl::Model( tMeshManager,
                                                      0,
                                                      tSetInfo,
                                                      0, false );

                // --------------------------------------------------------------------------------------
                // Define outputs

                vis::Output_Manager tOutputData;

                tOutputData.set_outputs( 0,
                                         vis::VIS_Mesh_Type::STANDARD,
                                         tPath,
                                         tVizMeshFileName,
                                         "./",
                                          "temp.exo",
                                         { "HMR_dummy_c_p2", "HMR_dummy_n_p2"},
                                         { "UX", "UY","UR_Analytic" },
                                         { vis::Field_Type::NODAL,vis::Field_Type::NODAL, vis::Field_Type::NODAL   },
                                         { vis::Output_Type::UX , vis::Output_Type::UY, vis::Output_Type::PROPERTY } );

                tModel->set_output_manager( &tOutputData );

                // --------------------------------------------------------------------------------------

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 1: create linear solver and algorithm
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                moris::Cell< enum MSI::Dof_Type > tDofTypesU( 2 );
                tDofTypesU( 0 ) = MSI::Dof_Type::UX;
                tDofTypesU( 1 ) = MSI::Dof_Type::UY;

                dla::Solver_Factory  tSolFactory;
                std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

//                tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_none;
//                tLinearSolverAlgorithm->set_param("AZ_output") = AZ_none;
//                tLinearSolverAlgorithm->set_param("AZ_max_iter") = 10000;
//                tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres;
//                tLinearSolverAlgorithm->set_param("AZ_subdomain_solve") = AZ_ilu;
//                tLinearSolverAlgorithm->set_param("AZ_graph_fill") = 10;
                //        tLinearSolverAlgorithm->set_param("ml_prec_type") = "SA";
                dla::Linear_Solver tLinSolver;
                tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                // STEP 2: create nonlinear solver and algorithm
                // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                NLA::Nonlinear_Solver_Factory tNonlinFactory;
                std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
                //        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithmMonolythicU = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

                tNonlinearSolverAlgorithm->set_param("NLA_max_iter")   = 5;
                //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_hard_break") = false;
                //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_max_lin_solver_restarts") = 2;
                //        tNonlinearSolverAlgorithmMonolythic->set_param("NLA_rebuild_jacobian") = true;

                tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
                //        tNonlinearSolverAlgorithmMonolythicU->set_linear_solver( &tLinSolver );

                NLA::Nonlinear_Solver tNonlinearSolverMain;
                tNonlinearSolverMain.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

                tNonlinearSolverMain.set_dof_type_list( tDofTypesU );

                // Create solver database
                sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );

                tNonlinearSolverMain.set_solver_warehouse( &tSolverWarehouse );

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
                tTimeSolver.set_output( 0, fOutputVizMesh );

                //------------------------------------------------------------------------------
                tTimeSolver.solve();

                moris::uint tNumDofs = tTimeSolver.get_solver_interface()->get_max_num_global_dofs();
                std::string tHDF5File  = tPath + tFileBase + ".hdf";
                hid_t tHDFId = create_hdf5_file(tHDF5File.c_str());
                herr_t tErr;
                save_scalar_to_hdf5_file(tHDFId,"dofs",tNumDofs,tErr);

                delete tModel;
//            }
//        }
    }

}
}

