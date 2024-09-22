/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_XTK_HMR_4_Material_Bar_Plane_Hole.cpp
 *
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
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
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

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
#include "cl_FEM_SP_Factory.hpp"         //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"         //FEM/INT/src
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
#include "cl_SOL_Warehouse.hpp"

#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Line.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Field.hpp"
#include "fn_norm.hpp"

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
Plane4MatMDL2( const moris::Matrix< moris::DDRMat >& aPoint )
{
    moris::real mXC = 1.4;
    moris::real mYC = 0.0;
    moris::real mNx = 1.0;
    moris::real mNy = 0.0;
    return ( mNx * ( aPoint( 0 ) - mXC ) + mNy * ( aPoint( 1 ) - mYC ) );
}

inline moris::real
Plane4MatMDL3( const moris::Matrix< moris::DDRMat >& aPoint )
{
    moris::real mXC = -1.4;
    moris::real mYC = 0.0;
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
tConstValFunction2MatMDL( moris::Matrix< moris::DDRMat >& aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >&         aParameters,
        moris::fem::Field_Interpolator_Manager*           aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

inline bool
tSolverOutputCriteria_4MatMDL( moris::tsa::Time_Solver* )
{
    return true;
}

inline void
run_hmr_for_multi_mat_model_2d(
        hmr::HMR&                                       aHMR,
        Vector< std::shared_ptr< moris::hmr::Field > >& aFields )
{
    moris_index tLagrangeMeshIndex = 0;

    std::shared_ptr< moris::hmr::Mesh > tMesh = aHMR.create_mesh( tLagrangeMeshIndex );

    aFields.resize( 2 );

    // create field
    aFields( 0 ) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    aFields( 1 ) = tMesh->create_field( "Geom", tLagrangeMeshIndex );

    aFields( 0 )->evaluate_scalar_function( Circle4MatMDL );
    aFields( 1 )->evaluate_scalar_function( Plane4MatMDL1 );

    // FIXME what is the following test about
    if ( false )
    {
        for ( uint k = 0; k < 1; ++k )
        {
            aHMR.flag_surface_elements_on_working_pattern( aFields( 0 ) );
            aHMR.flag_surface_elements_on_working_pattern( aFields( 1 ) );
            aHMR.perform_refinement_based_on_working_pattern( 0 );
            aFields( 0 )->evaluate_scalar_function( Circle4MatMDL );
            aFields( 1 )->evaluate_scalar_function( Plane4MatMDL1 );
        }
    }

    aHMR.finalize();
}

inline moris::real
MultiMat3dPlane( const moris::Matrix< moris::DDRMat >& aPoint )
{

    real mXn = 1.0;
    real mYn = 0.0;
    real mZn = 0.0;
    real mXc = 0.1;
    real mYc = 0.1;
    real mZc = 0.1;

    return mXn * ( aPoint( 0 ) - mXc ) + mYn * ( aPoint( 1 ) - mYc ) + mZn * ( aPoint( 2 ) - mZc );
}

inline moris::real
MultiMat3dCyl( const moris::Matrix< moris::DDRMat >& aPoint )
{
    moris::Matrix< moris::DDRMat > aCenter = { { 0.01 }, { 0.01 }, { 0.0 } };
    moris::Matrix< moris::DDRMat > aAxis   = { { 0.0 }, { 1.0 }, { 0.0 } };
    moris::real                    aRad    = 0.47334;
    moris::real                    aLength = 10;

    MORIS_ASSERT( aCenter.numel() == 3, "Centers need to have length 3" );
    MORIS_ASSERT( aAxis.numel() == 3, "axis need to have length 3" );

    Vector< moris::real > relativePosition = { ( aPoint( 0 ) - aCenter( 0 ) ), ( aPoint( 1 ) - aCenter( 1 ) ), ( aPoint( 2 ) - aCenter( 2 ) ) };
    moris::real           lsFromLeft       = ( relativePosition( 0 ) * ( -aAxis( 0 ) ) + relativePosition( 1 ) * ( -aAxis( 1 ) ) + relativePosition( 2 ) * ( -aAxis( 2 ) ) ) - aLength / 2.0;
    moris::real           lsFromRight      = ( relativePosition( 0 ) * ( aAxis( 0 ) ) + relativePosition( 1 ) * ( aAxis( 1 ) ) + relativePosition( 2 ) * ( aAxis( 2 ) ) ) - aLength / 2.0;

    moris::real           axialCrd  = ( relativePosition( 0 ) * ( aAxis( 0 ) ) + relativePosition( 1 ) * ( aAxis( 1 ) ) + relativePosition( 2 ) * ( aAxis( 2 ) ) );
    Vector< moris::real > radDir    = { ( relativePosition( 0 ) - aAxis( 0 ) * axialCrd ), ( relativePosition( 1 ) - aAxis( 1 ) * axialCrd ), ( relativePosition( 2 ) - aAxis( 2 ) * axialCrd ) };
    moris::real           radDist   = std::pow( radDir( 0 ) * radDir( 0 ) + radDir( 1 ) * radDir( 1 ) + radDir( 2 ) * radDir( 2 ), 0.5 );
    moris::real           lsFromRad = radDist - aRad;

    return std::max( std::max( lsFromLeft, lsFromRight ), lsFromRad );
}

inline real
MultiMat3dCylGeometry(
        const Matrix< DDRMat >& aCoordinates,
        const Vector< real >&   aParameters )
{
    return MultiMat3dCyl( aCoordinates );
}

inline void
run_hmr_for_multi_mat_model_3d(
        hmr::HMR&                                       aHMR,
        Vector< std::shared_ptr< moris::hmr::Field > >& aFields )
{
    moris_index tLagrangeMeshIndex = 0;

    std::shared_ptr< moris::hmr::Mesh > tMesh = aHMR.create_mesh( tLagrangeMeshIndex );

    aFields.resize( 2 );

    // create field
    aFields( 0 ) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    aFields( 0 )->evaluate_scalar_function( MultiMat3dCyl );
    aFields( 1 ) = tMesh->create_field( "Geom", tLagrangeMeshIndex );
    aFields( 1 )->evaluate_scalar_function( MultiMat3dPlane );
    for ( uint k = 0; k < 1; ++k )
    {
        aHMR.flag_surface_elements_on_working_pattern( aFields( 0 ) );
        aHMR.flag_surface_elements_on_working_pattern( aFields( 1 ) );
        aHMR.perform_refinement_based_on_working_pattern( 0 );
        aFields( 0 )->evaluate_scalar_function( MultiMat3dCyl );
        aFields( 1 )->evaluate_scalar_function( MultiMat3dPlane );
    }

    aFields( 0 )->evaluate_scalar_function( MultiMat3dCyl );
    aFields( 1 )->evaluate_scalar_function( MultiMat3dPlane );
    aHMR.finalize();
}

TEST_CASE( "XTK HMR 4 Material Bar Intersected By Plane and Hole", "[XTK_HMR_PLANE_BAR_HOLE_2D]" )
{

    if ( par_size() == 1 )
    {
        std::string tFieldName = "Geometry";

        moris::uint tLagrangeOrder     = 1;
        moris::uint tBsplineOrder      = 1;
        moris::uint tLagrangeMeshIndex = 0;

        moris::hmr::Parameters tParameters;
        tParameters.set_number_of_elements_per_dimension( { { 22 }, { 8 } } );
        tParameters.set_domain_dimensions( { { 6 }, { 2 } } );
        tParameters.set_domain_offset( { { -3.0 }, { -1.0 } } );
        tParameters.set_bspline_truncation( true );
        tParameters.set_output_meshes( { { 0 } } );
        tParameters.set_lagrange_orders( { { tLagrangeOrder } } );
        tParameters.set_lagrange_patterns( { { 0 } } );
        tParameters.set_bspline_orders( { { tBsplineOrder } } );
        tParameters.set_bspline_patterns( { { 0 } } );
        tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );
        tParameters.set_max_refinement_level( 2 );
        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );
        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 2 );
        tParameters.set_lagrange_to_bspline_mesh( { { { 0 } } } );

        hmr::HMR                                       tHMR( tParameters );
        Vector< std::shared_ptr< moris::hmr::Field > > tHMRFields;
        run_hmr_for_multi_mat_model_2d( tHMR, tHMRFields );

        hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        Vector< std::shared_ptr< gen::Geometry > > tGeometryVector( 2 );
        auto                                       tCircle = std::make_shared< gen::Circle >( 0.01, 0.01, 0.47334 );
        auto                                       tPlane  = std::make_shared< gen::Line >( 0.1, 0.1, 1.0, 0.0 );
        tGeometryVector( 0 )                               = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };
        tGeometryVector( 1 )                               = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

        size_t                          tModelDimension = 2;
        gen::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometryVector;
        gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );
        moris::xtk::Model    tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
        tXTKModel.mVerbose = false;

        Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
        tXTKModel.decompose( tDecompositionMethods );

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );

        moris::xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        moris::xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
        tPropConductivity1->set_parameters( { { { 1.0 } } } );
        tPropConductivity1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity1bis = std::make_shared< fem::Property >();
        tPropConductivity1bis->set_parameters( { { { 1.0 } } } );
        tPropConductivity1bis->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
        tPropConductivity2->set_parameters( { { { 0.1 } } } );
        tPropConductivity2->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity2bis = std::make_shared< fem::Property >();
        tPropConductivity2bis->set_parameters( { { { 0.1 } } } );
        tPropConductivity2bis->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { { { 5.0 } } } );
        tPropDirichlet->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { { { 20.0 } } } );
        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { { { 0.0 } } } );
        tPropTempLoad1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
        tPropTempLoad2->set_parameters( { { { 50.0 } } } );
        tPropTempLoad2->set_val_function( tConstValFunction2MatMDL );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
        tCMDiffLinIso1->set_space_dim( 2 );
        tCMDiffLinIso1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1bis =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1bis->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso1bis->set_property( tPropConductivity1bis, "Conductivity" );
        tCMDiffLinIso1bis->set_space_dim( 2 );
        tCMDiffLinIso1bis->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
        tCMDiffLinIso2->set_space_dim( 2 );
        tCMDiffLinIso2->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2bis =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2bis->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso2bis->set_property( tPropConductivity2bis, "Conductivity" );
        tCMDiffLinIso2bis->set_space_dim( 2 );
        tCMDiffLinIso2bis->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory                                 tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
        tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface1 =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface1->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface1->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface1->set_property( tPropConductivity2bis, "Material", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface2 =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface2->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface2->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface2->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface3 =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface3->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface3->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface3->set_property( tPropConductivity1bis, "Material", mtk::Leader_Follower::FOLLOWER );

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

        std::shared_ptr< fem::IWG > tIWGInterface1 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface1->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface1->set_stabilization_parameter( tSPNitscheInterface1, "NitscheInterface" );
        tIWGInterface1->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface1->set_constitutive_model( tCMDiffLinIso2bis, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::IWG > tIWGInterface2 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface2->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface2->set_stabilization_parameter( tSPNitscheInterface2, "NitscheInterface" );
        tIWGInterface2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface2->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::IWG > tIWGInterface3 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface3->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface3->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface3->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface3->set_stabilization_parameter( tSPNitscheInterface3, "NitscheInterface" );
        tIWGInterface3->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface3->set_constitutive_model( tCMDiffLinIso1bis, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        // create the IQIs
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQITEMP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIQITEMP->set_name( "IQI_TEMP" );

        // create a list of active block-sets
        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );
        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 2 );
        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name( 1, 3 );
        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name( 2, 3 );

        // std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName01)<<std::endl;
        // std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName02)<<std::endl;

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
        tSetBulk3.set_IWGs( { tIWGBulk2 } );
        tSetBulk3.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk4.set_IWGs( { tIWGBulk2 } );
        tSetBulk4.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk5;
        tSetBulk5.set_mesh_set_name( "HMR_dummy_c_p2" );
        tSetBulk5.set_IWGs( { tIWGBulk1 } );
        tSetBulk5.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk6;
        tSetBulk6.set_mesh_set_name( "HMR_dummy_n_p2" );
        tSetBulk6.set_IWGs( { tIWGBulk1 } );
        tSetBulk6.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk7;
        tSetBulk7.set_mesh_set_name( "HMR_dummy_c_p3" );
        tSetBulk7.set_IWGs( { tIWGBulk1 } );
        tSetBulk7.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk8;
        tSetBulk8.set_mesh_set_name( "HMR_dummy_n_p3" );
        tSetBulk8.set_IWGs( { tIWGBulk1 } );
        tSetBulk8.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p2" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_set_name( "SideSet_2_n_p3" );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetInterface1;
        tSetInterface1.set_mesh_set_name( tDblInterfaceSideSetName01 );
        tSetInterface1.set_IWGs( { tIWGInterface1 } );

        fem::Set_User_Info tSetInterface2;
        tSetInterface2.set_mesh_set_name( tDblInterfaceSideSetName02 );
        tSetInterface2.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface3;
        tSetInterface3.set_mesh_set_name( tDblInterfaceSideSetName13 );
        tSetInterface3.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface4;
        tSetInterface4.set_mesh_set_name( tDblInterfaceSideSetName23 );
        tSetInterface4.set_IWGs( { tIWGInterface3 } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 14 );
        tSetInfo( 0 )  = tSetBulk1;
        tSetInfo( 1 )  = tSetBulk2;
        tSetInfo( 2 )  = tSetBulk3;
        tSetInfo( 3 )  = tSetBulk4;
        tSetInfo( 4 )  = tSetBulk5;
        tSetInfo( 5 )  = tSetBulk6;
        tSetInfo( 6 )  = tSetBulk7;
        tSetInfo( 7 )  = tSetBulk8;
        tSetInfo( 8 )  = tSetDirichlet;
        tSetInfo( 9 )  = tSetNeumann;
        tSetInfo( 10 ) = tSetInterface1;
        tSetInfo( 11 ) = tSetInterface2;
        tSetInfo( 12 ) = tSetInterface3;
        tSetInfo( 13 ) = tSetInterface4;

        // create model
        mdl::Model* tModel = new mdl::Model( tMeshManager,
                0,
                tSetInfo,
                0,
                false );

        // --------------------------------------------------------------------------------------
        // Define outputs
        vis::Output_Manager tOutputData;
        std::string         tMeshOutputFile = "./mdl_exo/xtk_hmr_bar_plane_hole_bl_2_mat_l" + std::to_string( tLagrangeOrder ) + "_b" + std::to_string( tBsplineOrder ) + ".e";
        tOutputData.set_outputs( 0,
                vis::VIS_Mesh_Type::STANDARD,
                "./",
                tMeshOutputFile,
                "./",
                "temp.exo",
                { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_c_p2", "HMR_dummy_c_p3", "HMR_dummy_n_p0", "HMR_dummy_n_p1", "HMR_dummy_n_p2", "HMR_dummy_n_p3" },
                { "Temperature" },
                { vis::Field_Type::NODAL },
                { "IQI_TEMP" },
                { vis::Analysis_Type::FORWARD } );

        tModel->set_output_manager( &tOutputData );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory tSolFactory;
        Parameter_List      tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );
        tLinearSolverParameterList.set( "AZ_orthog", 1 );
        tLinearSolverParameterList.set( "AZ_solver", AZ_gmres_condnum );
        tLinearSolverParameterList.set( "AZ_precond", AZ_dom_decomp );
        tLinearSolverParameterList.set( "AZ_ilut_fill", 10.0 );
        tLinearSolverParameterList.set( "AZ_max_iter", 100 );
        tLinearSolverParameterList.set( "rel_residual", 1e-6 );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        Parameter_List                tNonlinearSolverParameterList = prm::create_nonlinear_algorithm_parameter_list();
        tNonlinearSolverParameterList.set( "NLA_max_iter", 2 );
        tNonlinearSolverParameterList.set( "NLA_hard_break", false );
        tNonlinearSolverParameterList.set( "NLA_max_lin_solver_restarts", 2 );
        tNonlinearSolverParameterList.set( "NLA_rebuild_jacobian", true );
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver( tNonlinearSolverParameterList );

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

        tTimeSolver.set_output( 0, tSolverOutputCriteria_4MatMDL );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        //
        //        // verify solution
        //        moris::Matrix<DDRMat> tGoldSolution =
        //        {{ 5.000000e+00},
        //         { 2.500000e+01},
        //         { 4.500000e+01},
        //         { 6.500000e+01},
        //         { 5.000000e+00},
        //         { 2.500000e+01},
        //         { 4.500000e+01},
        //         { 6.500000e+01}};
        //
        Matrix< DDRMat > tFullSol;
        tTimeSolver.get_full_solution( tFullSol );

        delete tModel;
        delete tInterpMesh;
    }
}

TEST_CASE( "XTK HMR 4 Material Bar Intersected By Plane and Hole 3D", "[XTK_HMR_PLANE_BAR_HOLE_3D]" )
{

    if ( par_size() == 1 )
    {
        std::string tFieldName = "Geometry";

        moris::uint tLagrangeOrder     = 1;
        moris::uint tBsplineOrder      = 1;
        moris::uint tLagrangeMeshIndex = 0;

        moris::hmr::Parameters tParameters;
        tParameters.set_number_of_elements_per_dimension( { { 22 }, { 8 }, { 8 } } );
        tParameters.set_domain_dimensions( { { 6 }, { 2 }, { 2 } } );
        tParameters.set_domain_offset( { { -3.0 }, { -1.0 }, { -1 } } );
        tParameters.set_bspline_truncation( true );
        tParameters.set_output_meshes( { { 0 } } );
        tParameters.set_lagrange_orders( { { tLagrangeOrder } } );
        tParameters.set_lagrange_patterns( { { 0 } } );
        tParameters.set_bspline_orders( { { tBsplineOrder } } );
        tParameters.set_bspline_patterns( { { 0 } } );
        tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );
        tParameters.set_max_refinement_level( 2 );
        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );
        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 2 );
        tParameters.set_lagrange_to_bspline_mesh( { { { 0 } } } );

        hmr::HMR                                       tHMR( tParameters );
        Vector< std::shared_ptr< moris::hmr::Field > > tHMRFields;
        run_hmr_for_multi_mat_model_3d( tHMR, tHMRFields );

        hmr::Interpolation_Mesh_HMR* tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        Vector< std::shared_ptr< gen::Geometry > > tGeometryVector( 2 );
        auto                                       tUserDefinedField = std::make_shared< gen::User_Defined_Field >( &( MultiMat3dCylGeometry ) );
        auto                                       tPlane            = std::make_shared< gen::Plane >( 0.1, 0.1, 0.1, 1.0, 0.0, 0.0 );
        tGeometryVector( 0 )                                         = { std::make_shared< gen::Level_Set_Geometry >( tUserDefinedField ) };
        tGeometryVector( 1 )                                         = { std::make_shared< gen::Level_Set_Geometry >( tPlane ) };

        size_t                          tModelDimension = 3;
        gen::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometryVector;
        gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );
        moris::xtk::Model    tXTKModel( tModelDimension, tInterpMesh, &tGeometryEngine );
        tXTKModel.mVerbose = false;

        Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
        tXTKModel.decompose( tDecompositionMethods );

        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0 );
        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        // Write mesh
        moris::xtk::Enriched_Integration_Mesh& tEnrIgMesh = tXTKModel.get_enriched_integ_mesh( 0 );

        moris_index tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set( 6, "ghost_ss_0" );
        tEnrIgMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_0", mtk::CellTopology::HEX8 );

        tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set( 7, "ghost_ss_1" );
        tEnrIgMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_1", mtk::CellTopology::HEX8 );

        tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set( 8, "ghost_ss_2" );
        tEnrIgMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_2", mtk::CellTopology::HEX8 );

        tSSIndex = tEnrIgMesh.create_side_set_from_dbl_side_set( 9, "ghost_ss_3" );
        tEnrIgMesh.create_block_set_from_cells_of_side_set( tSSIndex, "ghost_bs_3", mtk::CellTopology::HEX8 );

        moris::mtk::Writer_Exodus writer( &tEnrIgMesh );
        writer.write_mesh( "", "./mdl_exo/xtk_hmr_bar_hole_integ_ghost.exo", "", "temp.exo" );

        // Write the fields
        writer.set_time( 0.0 );
        writer.close_file();

        moris::xtk::Enriched_Interpolation_Mesh& tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        moris::xtk::Enriched_Integration_Mesh&   tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
        tMeshManager->register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

        //------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropConductivity1 = std::make_shared< fem::Property >();
        tPropConductivity1->set_parameters( { { { 1.0 } } } );
        tPropConductivity1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity1bis = std::make_shared< fem::Property >();
        tPropConductivity1bis->set_parameters( { { { 1.0 } } } );
        tPropConductivity1bis->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity2 = std::make_shared< fem::Property >();
        tPropConductivity2->set_parameters( { { { 0.1 } } } );
        tPropConductivity2->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropConductivity2bis = std::make_shared< fem::Property >();
        tPropConductivity2bis->set_parameters( { { { 0.1 } } } );
        tPropConductivity2bis->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropDirichlet = std::make_shared< fem::Property >();
        tPropDirichlet->set_parameters( { { { 5.0 } } } );
        tPropDirichlet->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropNeumann = std::make_shared< fem::Property >();
        tPropNeumann->set_parameters( { { { 20.0 } } } );
        tPropNeumann->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad1 = std::make_shared< fem::Property >();
        tPropTempLoad1->set_parameters( { { { 0.0 } } } );
        tPropTempLoad1->set_val_function( tConstValFunction2MatMDL );

        std::shared_ptr< fem::Property > tPropTempLoad2 = std::make_shared< fem::Property >();
        tPropTempLoad2->set_parameters( { { { 50.0 } } } );
        tPropTempLoad2->set_val_function( tConstValFunction2MatMDL );

        // define constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso1->set_property( tPropConductivity1, "Conductivity" );
        tCMDiffLinIso1->set_space_dim( 3 );
        tCMDiffLinIso1->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso1bis =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso1bis->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso1bis->set_property( tPropConductivity1bis, "Conductivity" );
        tCMDiffLinIso1bis->set_space_dim( 3 );
        tCMDiffLinIso1bis->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2 =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso2->set_property( tPropConductivity2, "Conductivity" );
        tCMDiffLinIso2->set_space_dim( 3 );
        tCMDiffLinIso2->set_local_properties();

        std::shared_ptr< fem::Constitutive_Model > tCMDiffLinIso2bis =
                tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
        tCMDiffLinIso2bis->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tCMDiffLinIso2bis->set_property( tPropConductivity2bis, "Conductivity" );
        tCMDiffLinIso2bis->set_space_dim( 3 );
        tCMDiffLinIso2bis->set_local_properties();

        // define stabilization parameters
        fem::SP_Factory                                 tSPFactory;
        std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
                tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tSPDirichletNitsche->set_parameters( { { { 1.0 } } } );
        tSPDirichletNitsche->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface1 =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface1->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface1->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface1->set_property( tPropConductivity2bis, "Material", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface2 =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface2->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface2->set_property( tPropConductivity2, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface2->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface3 =
                tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
        tSPNitscheInterface3->set_parameters( { { { 1.0 } } } );
        tSPNitscheInterface3->set_property( tPropConductivity1, "Material", mtk::Leader_Follower::LEADER );
        tSPNitscheInterface3->set_property( tPropConductivity1bis, "Material", mtk::Leader_Follower::FOLLOWER );

        // create the IQIs
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQITEMP =
                tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQITEMP->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQITEMP->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
        tIQITEMP->set_name( "IQI_TEMP" );

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

        std::shared_ptr< fem::IWG > tIWGInterface1 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface1->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface1->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface1->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface1->set_stabilization_parameter( tSPNitscheInterface1, "NitscheInterface" );
        tIWGInterface1->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface1->set_constitutive_model( tCMDiffLinIso2bis, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::IWG > tIWGInterface2 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface2->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface2->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface2->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface2->set_stabilization_parameter( tSPNitscheInterface2, "NitscheInterface" );
        tIWGInterface2->set_constitutive_model( tCMDiffLinIso2, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface2->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        std::shared_ptr< fem::IWG > tIWGInterface3 =
                tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE );
        tIWGInterface3->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface3->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
        tIWGInterface3->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::FOLLOWER );
        tIWGInterface3->set_stabilization_parameter( tSPNitscheInterface3, "NitscheInterface" );
        tIWGInterface3->set_constitutive_model( tCMDiffLinIso1, "Diffusion", mtk::Leader_Follower::LEADER );
        tIWGInterface3->set_constitutive_model( tCMDiffLinIso1bis, "Diffusion", mtk::Leader_Follower::FOLLOWER );

        // create a list of active block-sets
        std::string tDblInterfaceSideSetName01 = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 );
        std::string tDblInterfaceSideSetName02 = tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 2 );
        std::string tDblInterfaceSideSetName13 = tEnrIntegMesh.get_dbl_interface_side_set_name( 1, 3 );
        std::string tDblInterfaceSideSetName23 = tEnrIntegMesh.get_dbl_interface_side_set_name( 2, 3 );

        // std::cout<<"tDblInterfaceSideSetName01 = "<<tDblInterfaceSideSetName01<<" | Index = "<<tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName01)<<std::endl;
        // std::cout<<"tDblInterfaceSideSetName02 = "<<tDblInterfaceSideSetName02<<" | Index = "<<tEnrIntegMesh.get_set_index_by_name(tDblInterfaceSideSetName02)<<std::endl;

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
        tSetBulk3.set_IWGs( { tIWGBulk2 } );
        tSetBulk3.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk4;
        tSetBulk4.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk4.set_IWGs( { tIWGBulk2 } );
        tSetBulk4.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk5;
        tSetBulk5.set_mesh_set_name( "HMR_dummy_c_p2" );
        tSetBulk5.set_IWGs( { tIWGBulk1 } );
        tSetBulk5.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk6;
        tSetBulk6.set_mesh_set_name( "HMR_dummy_n_p2" );
        tSetBulk6.set_IWGs( { tIWGBulk1 } );
        tSetBulk6.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk7;
        tSetBulk7.set_mesh_set_name( "HMR_dummy_c_p3" );
        tSetBulk7.set_IWGs( { tIWGBulk1 } );
        tSetBulk7.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetBulk8;
        tSetBulk8.set_mesh_set_name( "HMR_dummy_n_p3" );
        tSetBulk8.set_IWGs( { tIWGBulk1 } );
        tSetBulk8.set_IQIs( { tIQITEMP } );

        fem::Set_User_Info tSetDirichlet;
        tSetDirichlet.set_mesh_set_name( "SideSet_4_n_p2" );
        tSetDirichlet.set_IWGs( { tIWGDirichlet } );

        fem::Set_User_Info tSetNeumann;
        tSetNeumann.set_mesh_set_name( "SideSet_2_n_p3" );
        tSetNeumann.set_IWGs( { tIWGNeumann } );

        fem::Set_User_Info tSetInterface1;
        tSetInterface1.set_mesh_set_name( tDblInterfaceSideSetName01 );
        tSetInterface1.set_IWGs( { tIWGInterface1 } );

        fem::Set_User_Info tSetInterface2;
        tSetInterface2.set_mesh_set_name( tDblInterfaceSideSetName02 );
        tSetInterface2.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface3;
        tSetInterface3.set_mesh_set_name( tDblInterfaceSideSetName13 );
        tSetInterface3.set_IWGs( { tIWGInterface2 } );

        fem::Set_User_Info tSetInterface4;
        tSetInterface4.set_mesh_set_name( tDblInterfaceSideSetName23 );
        tSetInterface4.set_IWGs( { tIWGInterface3 } );

        // create a cell of set info
        Vector< fem::Set_User_Info > tSetInfo( 14 );
        tSetInfo( 0 )  = tSetBulk1;
        tSetInfo( 1 )  = tSetBulk2;
        tSetInfo( 2 )  = tSetBulk3;
        tSetInfo( 3 )  = tSetBulk4;
        tSetInfo( 4 )  = tSetBulk5;
        tSetInfo( 5 )  = tSetBulk6;
        tSetInfo( 6 )  = tSetBulk7;
        tSetInfo( 7 )  = tSetBulk8;
        tSetInfo( 8 )  = tSetDirichlet;
        tSetInfo( 9 )  = tSetNeumann;
        tSetInfo( 10 ) = tSetInterface1;
        tSetInfo( 11 ) = tSetInterface2;
        tSetInfo( 12 ) = tSetInterface3;
        tSetInfo( 13 ) = tSetInterface4;

        // create model
        mdl::Model* tModel = new mdl::Model( tMeshManager,
                0,
                tSetInfo,
                0,
                false );

        Vector< enum MSI::Dof_Type > tDofTypes1( 1, MSI::Dof_Type::TEMP );

        // --------------------------------------------------------------------------------------
        // Define outputs
        vis::Output_Manager tOutputData;
        std::string         tMeshOutputFile = "xtk_hmr_bar_plane_hole_3d_l" + std::to_string( tLagrangeOrder ) + "_b" + std::to_string( tBsplineOrder ) + ".e";
        tOutputData.set_outputs( 0,
                vis::VIS_Mesh_Type::STANDARD,
                "./",
                tMeshOutputFile,
                "./",
                "temp.exo",
                { "HMR_dummy_c_p0", "HMR_dummy_c_p1", "HMR_dummy_c_p2", "HMR_dummy_c_p3", "HMR_dummy_n_p0", "HMR_dummy_n_p1", "HMR_dummy_n_p2", "HMR_dummy_n_p3" },
                { "Temperature" },
                { vis::Field_Type::NODAL },
                { "IQI_TEMP" },
                { vis::Analysis_Type::FORWARD } );

        tModel->set_output_manager( &tOutputData );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 1: create linear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        dla::Solver_Factory tSolFactory;
        Parameter_List      tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_aztec();
        tLinearSolverParameterList.set( "AZ_diagnostics", AZ_none );
        tLinearSolverParameterList.set( "AZ_output", AZ_none );
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm = tSolFactory.create_solver( tLinearSolverParameterList );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 2: create nonlinear solver and algorithm
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        NLA::Nonlinear_Solver_Factory               tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm = tNonlinFactory.create_nonlinear_solver();

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

        tTimeSolver.set_output( 0, tSolverOutputCriteria_4MatMDL );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // STEP 4: Solve and check
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        tTimeSolver.solve();

        ////
        ////        // verify solution
        ////        moris::Matrix<DDRMat> tGoldSolution =
        ////        {{ 5.000000e+00},
        ////         { 2.500000e+01},
        ////         { 4.500000e+01},
        ////         { 6.500000e+01},
        ////         { 5.000000e+00},
        ////         { 2.500000e+01},
        ////         { 4.500000e+01},
        ////         { 6.500000e+01}};
        ////
        Matrix< DDRMat > tFullSol;
        tTimeSolver.get_full_solution( tFullSol );

        // print(tFullSol,"tFullSol");

        //        // Declare the fields related to enrichment strategy in output options
        //        // output solution and meshes
        //        moris::xtk::Output_Options tOutputOptions;
        //        tOutputOptions.mAddNodeSets = false;
        //        tOutputOptions.mAddSideSets = true;
        //        tOutputOptions.mAddClusters = false;
        //
        //        // add solution field to integration mesh
        //        std::string tIntegSolFieldName = "solution";
        //        tOutputOptions.mRealNodeExternalFieldNames = {tIntegSolFieldName};
        //
        //        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
        //
        //        // Write to Integration mesh for visualization
        //        Matrix<DDRMat> tIntegSol = tModel->get_solution_for_integration_mesh_output( MSI::Dof_Type::TEMP );
        //
        //
        //        Matrix<DDRMat> tSTKIntegSol(tIntegMesh1->get_num_entities(EntityRank::NODE),1);
        //
        //        for(moris::uint i = 0; i < tIntegMesh1->get_num_entities(EntityRank::NODE); i++)
        //        {
        //            moris::moris_id tID = tIntegMesh1->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
        //            moris::moris_index tMyIndex = tEnrIntegMesh.get_loc_entity_ind_from_entity_glb_id(tID,EntityRank::NODE);
        //
        //            tSTKIntegSol(i) = tIntegSol(tMyIndex);
        //        }
        //

        delete tModel;
        delete tInterpMesh;
    }
}
