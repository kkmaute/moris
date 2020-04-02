/*
 * UT_MDL_Fluid_Benchmark.cpp
 *
 *  Created on: Mar 28, 2020
 *      Author: noel
 */

#include "catch.hpp"

#include "cl_Geom_Field.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
#include "cl_Plane.hpp"
//PRM
#include "cl_PRM_HMR_Parameters.hpp"
#include "cl_PRM_SOL_Parameters.hpp"
//MTK/src
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
//XTK/src
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Ghost_Stabilization.hpp"
//HMR/src
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Element.hpp"
#include "cl_HMR_Factory.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Parameters.hpp"
//FEM/INT/src
#include "cl_FEM_Element_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
//FEM/MDL/src
#include "cl_MDL_Model.hpp"
//FEM/VIS/src
#include "cl_VIS_Factory.hpp"
#include "cl_VIS_Visualization_Mesh.hpp"
#include "cl_VIS_Output_Manager.hpp"
//FEM/MSI/src
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_MSI_Model_Solver_Interface.hpp"
//SOL/src
#include "cl_SOL_Warehouse.hpp"
//SOL/DLA/src
#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver.hpp"
//SOL/NLA/src
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_NLA_Nonlinear_Solver.hpp"
#include "cl_NLA_Nonlinear_Problem.hpp"
//SOL/TSA/src
#include "cl_TSA_Time_Solver_Factory.hpp"
#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Time_Solver.hpp"
//GEN
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_Geometry.hpp"
#include <functional>
#include "../../../GEN/GEN_MAIN/src/geometry/cl_GEN_Geom_Field_HMR.hpp"

namespace moris
{

//-------------------------------------------------------------------------------------
void ConstFuncVal_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void InletVelocityFunc_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // unpack parameters
    real tRadiusChannel = aParameters( 0 )( 0 );
    real tYChannel      = aParameters( 1 )( 0 );

    // get position in space
    real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

//    std::cout<<"x"<<aFIManager->get_IP_geometry_interpolator()->valx()( 0 )<<std::endl;
//    std::cout<<"y"<<aFIManager->get_IP_geometry_interpolator()->valx()( 1 )<<std::endl;

    // set size for aPropMatrix
    aPropMatrix.set_size( 2, 1, 0.0 );

    // velocity along x direction
    aPropMatrix( 0 ) = - ( tY - ( tYChannel + tRadiusChannel ) )
                       * ( tY - ( tYChannel - tRadiusChannel ) ) / ( 2.0 * std::pow( tRadiusChannel, 2.0 ) );
}

void FSVelocityFunc_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // get space dim
    uint tSpaceDim = aFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

    // set size for aPropMatrix
    aPropMatrix.set_size( tSpaceDim, 1, 0.0 );
}

void InletPressureFunc_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // set size for aPropMatrix
    aPropMatrix.set_size( 1, 1, aParameters( 0 )( 0 ) );
}


bool tSolverOutputCriteria_MDLFluidBench( moris::tsa::Time_Solver * )
{
    return true;
}

//-------------------------------------------------------------------------------------
TEST_CASE("MDL_Fluid_Benchmark_Immersed_Inlet_Velocity","[MDL_Fluid_Benchmark_Immersed_Inlet_Velocity]")
{
    if( par_size() <= 1 )
    {
        // Geometry Parameters
        moris::real tDomainLX = 10.0;                   /* Length of full domain in x (m) */
        moris::real tDomainLY = 10.0;                   /* Length of full domain in y (m) */
        Matrix<DDRMat> tCenterPoint = { { 0.37, 0.26 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
        moris::real tPlaneBottom = -0.5;                  /* y bottom plane (m) */
        moris::real tPlaneTop    =  0.5;                  /* y top plane    (m) */
        moris::real tPlaneLeft   = -2.5;                /* x left plane   (m) */
        moris::real tPlaneRight  =  2.5;                /* x right plane  (m) */
        moris::real tChannelRadius = ( tPlaneTop - tPlaneBottom ) / 2.0; /* channel radius  (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0; /* Inlet pressure  () */
        moris::real tOutletPressure =  0.0; /* Outlet pressure () */
        moris::real tGammaNitsche   = 1000.0;  /* Penalty for Dirichlet BC */

        // Mesh Setup
        moris::uint tNumX   = 200; /* Number of elements in x*/
        moris::uint tNumY   = 200; /* Number of elements in y*/
        moris::uint tNumRef = 1;    /* Number of HMR refinements */
        moris::uint tOrder  = 1;    /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint tLagrangeMeshIndex = 0;
        std::string tPlaneBottomFieldName = "PlaneBottom";
        std::string tPlaneTopFieldName    = "PlaneTop";
        std::string tPlaneLeftFieldName   = "PlaneLeft";
        std::string tPlaneRightFieldName  = "PlaneRight";

        ParameterList tParameters = prm::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2+tCenterPoint(0)) + "," + std::to_string(-tDomainLY/2+tCenterPoint(1)) );
        tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
        tParameters.set( "lagrange_output_meshes", std::string("0") );

        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", 0 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //  create field
        std::shared_ptr< moris::hmr::Field > tPlaneBottomField = tMesh->create_field( tPlaneBottomFieldName, tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tPlaneTopField    = tMesh->create_field( tPlaneTopFieldName,    tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tPlaneLeftField   = tMesh->create_field( tPlaneLeftFieldName,   tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tPlaneRightField  = tMesh->create_field( tPlaneRightFieldName,  tLagrangeMeshIndex );

        for( uint k = 0; k < tNumRef; k++ )
        {
            moris::ge::Plane< 2 > tPlane00( {{ 0.0, tPlaneBottom }}, {{ 0.0, 1.0 }} );
            moris::ge::Plane< 2 > tPlane01( {{ 0.0, tPlaneTop }},    {{ 0.0, 1.0 }} );
            moris::ge::Plane< 2 > tPlane02( {{ tPlaneLeft, 0.0 }},   {{ 1.0, 0.0 }} );
            moris::ge::Plane< 2 > tPlane03( {{ tPlaneRight, 0.0 }},  {{ 1.0, 0.0 }} );
            moris::Cell< moris::ge::GEN_Geometry* > tGeomVec = { &tPlane00, &tPlane01, &tPlane02, &tPlane03 };

            moris::ge::GEN_Phase_Table     tPhaseTable( tGeomVec.size(),  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine tGENGeometryEngine( tGeomVec, tPhaseTable, 2 );

            moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );

            uint tNumIPNodes = tMesh->get_num_nodes();
            Matrix< DDRMat > tFieldData0( tNumIPNodes,1 );
            Matrix< DDRMat > tFieldData1( tNumIPNodes,1 );
            Matrix< DDRMat > tFieldData2( tNumIPNodes,1 );
            Matrix< DDRMat > tFieldData3( tNumIPNodes,1 );

            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumIPNodes );
            Matrix< DDRMat > tCoords( tNumIPNodes, 2 );
            for( uint i = 0; i < tNumIPNodes; i++ )
            {
                tCoords.set_row( i, tMesh->get_mtk_vertex(i).get_coords() );
            }

            tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );

            for( uint i = 0; i < tNumIPNodes; i++ )
            {
                tFieldData0( i ) = tGENGeometryEngine.get_entity_phase_val( i, 0 );
                tFieldData1( i ) = tGENGeometryEngine.get_entity_phase_val( i, 1 );
                tFieldData2( i ) = tGENGeometryEngine.get_entity_phase_val( i, 2 );
                tFieldData3( i ) = tGENGeometryEngine.get_entity_phase_val( i, 3 );
            }

            tHMR.based_on_field_put_elements_on_queue( tFieldData0, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData1, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData2, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData3, tLagrangeMeshIndex );

            tHMR.perform_refinement_based_on_working_pattern( 0, false );
        }
        tHMR.finalize();
//        tHMR.save_to_exodus( 0, tHMRIPMeshFileName );

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        //-----------------------------------------------------------------------------------------------

        moris::ge::Plane< 2 > tPlane0( {{ 0.0, tPlaneBottom }}, {{ 0.0, 1.0 }} );
        moris::ge::Plane< 2 > tPlane1( {{ 0.0, tPlaneTop }},    {{ 0.0, 1.0 }} );
        moris::ge::Plane< 2 > tPlane2( {{ tPlaneLeft, 0.0 }},   {{ 1.0, 0.0 }} );
        moris::ge::Plane< 2 > tPlane3( {{ tPlaneRight, 0.0 }},  {{ 1.0, 0.0 }} );

        // NOTE the order of this geometry vector is important.
        // If it changes the resulting bulk phase of the output mesh change.
        moris::Cell<moris::ge::GEN_Geometry*> tGeomVec0 = { &tPlane0, &tPlane1, &tPlane2, &tPlane3 };

        size_t tModelDimension = 2;
        moris::ge::GEN_Phase_Table     tPhaseTable0( tGeomVec0.size(), Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::GEN_Geometry_Engine tGENGeometryEngine0( tGeomVec0, tPhaseTable0, tModelDimension );

        // --------------------------------------------------------------------------------------
        xtk::Model tXTKModel(tModelDimension,tInterpolationMesh.get(),&tGENGeometryEngine0);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
//        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // output integration mesh
        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
        std::string tOutputFile = "./mdl_exo/MDL_Fluid_Benchmark_Immersed_Inlet_Velocity.exo";
        tIntegMesh1->create_output_mesh( tOutputFile );

        // get meshes for FEM
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create for fem
        // --------------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropInletVelocity = std::make_shared< fem::Property >();
        tPropInletVelocity->set_parameters( { {{ tChannelRadius }}, {{ 0.0 }} } );
        tPropInletVelocity->set_val_function( InletVelocityFunc_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFSVelocity = std::make_shared< fem::Property >();
        tPropFSVelocity->set_val_function( FSVelocityFunc_MDLFluidBench );

        // create constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMFluid
        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
        tCMFluid->set_property( tPropFluidDensity, "Density" );
        tCMFluid->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_parameters( { {{ 36.0 }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGVelocityBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGPressureBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGInletVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGInletVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGInletVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletVelocity->set_property( tPropInletVelocity, "Dirichlet" );
        tIWGInletVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGInletVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGInletPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGInletPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGInletPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletPressure->set_property( tPropInletVelocity, "Dirichlet" );
        tIWGInletPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        std::shared_ptr< fem::IWG > tIWGFSVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGFSVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGFSVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSVelocity->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGFSVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGFSPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGFSPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGFSPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSPressure->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVX->set_output_type( vis::Output_Type::VX );
        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVX->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVY->set_output_type( vis::Output_Type::VY );
        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVY->set_output_type_index( 1 );

        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIP->set_output_type( vis::Output_Type::P );
        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIQIP->set_output_type_index( 0 );

        // create set info
        // --------------------------------------------------------------------------------------
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p10" );
        tSetBulk1.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk1.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p10" );
        tSetBulk2.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk2.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );

        // Fluid/solid bottom
        fem::Set_User_Info tSetFSBottom;
        tSetFSBottom.set_mesh_set_name( "iside_g_0_b0_10_b1_2" );
        tSetFSBottom.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid top
        fem::Set_User_Info tSetFSTop;
        tSetFSTop.set_mesh_set_name( "iside_g_1_b0_10_b1_14" );
        tSetFSTop.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Inlet
        fem::Set_User_Info tSetInlet;
        tSetInlet.set_mesh_set_name( "iside_g_2_b0_10_b1_8" );
        tSetInlet.set_IWGs( { tIWGInletVelocity, tIWGInletPressure } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 5 );
        tSetInfo( 0 )  = tSetBulk1;
        tSetInfo( 1 )  = tSetBulk2;
        tSetInfo( 2 )  = tSetInlet;
        tSetInfo( 3 )  = tSetFSBottom;
        tSetInfo( 4 )  = tSetFSTop;

        // create model
        // --------------------------------------------------------------------------------------
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // define outputs
        // --------------------------------------------------------------------------------------
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                                 vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                                 "MDL_Fluid_Benchmark_Immersed_Inlet_Velocity_Output.exo",
                                 { "HMR_dummy_c_p10", "HMR_dummy_n_p10" },
                                 { "VX", "VY", "P" },
                                 { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                                 { vis::Output_Type::VX,  vis::Output_Type::VY, vis::Output_Type::P } );
        tModel->set_output_manager( &tOutputData );

        // --------------------------------------------------------------------------------------
        // define linear solver and algorithm
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
        = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

 //       tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_all;
 //       tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
 //       tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // --------------------------------------------------------------------------------------
        // define nonlinear solver and algorithm
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // --------------------------------------------------------------------------------------
        // define time solver and algorithm
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
        = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
        tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );

        tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );

        // --------------------------------------------------------------------------------------
        // solve and check
        tTimeSolver.solve();

//        // create linear solver and algorithm
//        // --------------------------------------------------------------------------------------
//        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
//
//        moris::Cell< moris::Cell< moris::ParameterList > > tParameterlist( 7 );
//        for( uint Ik = 0; Ik < 7; Ik ++)
//        {
//            tParameterlist( Ik ).resize(1);
//        }
//
//        tParameterlist( 0 )(0) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS2_IMPL );
////        tParameterlist( 0 )(0).set("AZ_diagnostics"    , AZ_none  );
////        tParameterlist( 0 )(0).set("AZ_output"         , AZ_none  );
////        tParameterlist( 0 )(0).set("AZ_max_iter"       , 10000    );
////        tParameterlist( 0 )(0).set("AZ_solver"         , AZ_gmres );
////        tParameterlist( 0 )(0).set("AZ_subdomain_solve", AZ_ilu   );
////        tParameterlist( 0 )(0).set("AZ_graph_fill"     , 10       );
////        tParameterlist( 0 )(0).set("Use_ML_Prec"       ,  true    );
//
//        tParameterlist( 1 )(0) = moris::prm::create_linear_solver_parameter_list();
//        tParameterlist( 2 )(0) = moris::prm::create_nonlinear_algorithm_parameter_list();
//        tParameterlist( 3 )(0) = moris::prm::create_nonlinear_solver_parameter_list();
//        tParameterlist( 3 )(0).set("NLA_DofTypes"      , std::string("VX,VY;P") );
//
//        tParameterlist( 4 )(0) = moris::prm::create_time_solver_algorithm_parameter_list();
//        tParameterlist( 5 )(0) = moris::prm::create_time_solver_parameter_list();
//        tParameterlist( 5 )(0).set("TSA_DofTypes"      , std::string("VX,VY;P") );
//
//        tParameterlist( 6 )(0) = moris::prm::create_solver_warehouse_parameterlist();
//
//        tSolverWarehouse.set_parameterlist( tParameterlist );
//
//        tSolverWarehouse.initialize();
//
//        tsa::Time_Solver * tTimeSolver = tSolverWarehouse.get_main_time_solver();
//
//        tTimeSolver->set_output( 0, tSolverOutputCriteria_MDLFluidBench );
//
//        tTimeSolver->solve();


//         //print full solution for debug
//         Matrix<DDRMat> tFullSol;
//         tTimeSolver->get_full_solution(tFullSol);
//         print(tFullSol,"tFullSol");

        // clean up
        //------------------------------------------------------------------------------
        delete tIntegMesh1;
        delete tModel;
    }
}

//-------------------------------------------------------------------------------------
TEST_CASE("MDL_Fluid_Benchmark_Immersed_Inlet_Pressure","[MDL_Fluid_Benchmark_Immersed_Inlet_Pressure]")
{
    if( par_size() <= 1 )
    {
        // Geometry Parameters
        moris::real tDomainLX = 10.0;                   /* Length of full domain in x (m) */
        moris::real tDomainLY = 10.0;                   /* Length of full domain in y (m) */
        Matrix<DDRMat> tCenterPoint = { { 0.37, 0.26 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
        moris::real tPlaneBottom = -0.5;                  /* y bottom plane (m) */
        moris::real tPlaneTop    =  0.5;                  /* y top plane    (m) */
        moris::real tPlaneLeft   = -2.5;                /* x left plane   (m) */
        moris::real tPlaneRight  =  2.5;                /* x right plane  (m) */
        moris::real tChannelRadius = ( tPlaneTop - tPlaneBottom ) / 2.0; /* channel radius  (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0; /* Inlet pressure  () */
        moris::real tOutletPressure =  0.0; /* Outlet pressure () */
        moris::real tGammaNitsche   = 1000.0;  /* Penalty for Dirichlet BC */

        // Mesh Setup
        moris::uint tNumX   = 400; /* Number of elements in x*/
        moris::uint tNumY   = 400; /* Number of elements in y*/
        moris::uint tNumRef = 1;    /* Number of HMR refinements */
        moris::uint tOrder  = 1;    /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint tLagrangeMeshIndex = 0;
        std::string tPlaneBottomFieldName = "PlaneBottom";
        std::string tPlaneTopFieldName    = "PlaneTop";
        std::string tPlaneLeftFieldName   = "PlaneLeft";
        std::string tPlaneRightFieldName  = "PlaneRight";

        ParameterList tParameters = prm::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2+tCenterPoint(0)) + "," + std::to_string(-tDomainLY/2+tCenterPoint(1)) );
        tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
        tParameters.set( "lagrange_output_meshes", std::string("0") );

        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", 0 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //  create field
        std::shared_ptr< moris::hmr::Field > tPlaneBottomField = tMesh->create_field( tPlaneBottomFieldName, tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tPlaneTopField    = tMesh->create_field( tPlaneTopFieldName,    tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tPlaneLeftField   = tMesh->create_field( tPlaneLeftFieldName,   tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tPlaneRightField  = tMesh->create_field( tPlaneRightFieldName,  tLagrangeMeshIndex );

        for( uint k = 0; k < tNumRef; k++ )
        {
            moris::ge::Plane< 2 > tPlane00( {{ 0.0, tPlaneBottom }}, {{ 0.0, 1.0 }} );
            moris::ge::Plane< 2 > tPlane01( {{ 0.0, tPlaneTop }},    {{ 0.0, 1.0 }} );
            moris::ge::Plane< 2 > tPlane02( {{ tPlaneLeft, 0.0 }},   {{ 1.0, 0.0 }} );
            moris::ge::Plane< 2 > tPlane03( {{ tPlaneRight, 0.0 }},  {{ 1.0, 0.0 }} );
            moris::Cell< moris::ge::GEN_Geometry* > tGeomVec = { &tPlane00, &tPlane01, &tPlane02, &tPlane03 };

            moris::ge::GEN_Phase_Table     tPhaseTable( tGeomVec.size(),  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine tGENGeometryEngine( tGeomVec, tPhaseTable, 2 );

            moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );

            uint tNumIPNodes = tMesh->get_num_nodes();
            Matrix< DDRMat > tFieldData0( tNumIPNodes,1 );
            Matrix< DDRMat > tFieldData1( tNumIPNodes,1 );
            Matrix< DDRMat > tFieldData2( tNumIPNodes,1 );
            Matrix< DDRMat > tFieldData3( tNumIPNodes,1 );

            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumIPNodes );
            Matrix< DDRMat > tCoords( tNumIPNodes, 2 );
            for( uint i = 0; i < tNumIPNodes; i++ )
            {
                tCoords.set_row( i, tMesh->get_mtk_vertex(i).get_coords() );
            }

            tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );

            for( uint i = 0; i < tNumIPNodes; i++ )
            {
                tFieldData0( i ) = tGENGeometryEngine.get_entity_phase_val( i, 0 );
                tFieldData1( i ) = tGENGeometryEngine.get_entity_phase_val( i, 1 );
                tFieldData2( i ) = tGENGeometryEngine.get_entity_phase_val( i, 2 );
                tFieldData3( i ) = tGENGeometryEngine.get_entity_phase_val( i, 3 );
            }

            tHMR.based_on_field_put_elements_on_queue( tFieldData0, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData1, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData2, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData3, tLagrangeMeshIndex );

            tHMR.perform_refinement_based_on_working_pattern( 0, false );
        }
        tHMR.finalize();
//        tHMR.save_to_exodus( 0, tHMRIPMeshFileName );

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        //-----------------------------------------------------------------------------------------------

        moris::ge::Plane< 2 > tPlane0( {{ 0.0, tPlaneBottom }}, {{ 0.0, 1.0 }} );
        moris::ge::Plane< 2 > tPlane1( {{ 0.0, tPlaneTop }},    {{ 0.0, 1.0 }} );
        moris::ge::Plane< 2 > tPlane2( {{ tPlaneLeft, 0.0 }},   {{ 1.0, 0.0 }} );
        moris::ge::Plane< 2 > tPlane3( {{ tPlaneRight, 0.0 }},  {{ 1.0, 0.0 }} );

        // NOTE the order of this geometry vector is important.
        // If it changes the resulting bulk phase of the output mesh change.
        moris::Cell<moris::ge::GEN_Geometry*> tGeomVec0 = { &tPlane0, &tPlane1, &tPlane2, &tPlane3 };

        size_t tModelDimension = 2;
        moris::ge::GEN_Phase_Table     tPhaseTable0( tGeomVec0.size(), Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::GEN_Geometry_Engine tGENGeometryEngine0( tGeomVec0, tPhaseTable0, tModelDimension );

        // --------------------------------------------------------------------------------------
        xtk::Model tXTKModel(tModelDimension,tInterpolationMesh.get(),&tGENGeometryEngine0);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
//        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // output integration mesh
        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
        std::string tOutputFile = "./mdl_exo/MDL_Fluid_Benchmark_Immersed_Inlet_Pressure.exo";
        tIntegMesh1->create_output_mesh( tOutputFile );

        // get meshes for FEM
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create for fem
        // --------------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropInletPressure = std::make_shared< fem::Property >();
        tPropInletPressure->set_parameters( { {{ tInletPressure }} } );
        tPropInletPressure->set_val_function( InletPressureFunc_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFSVelocity = std::make_shared< fem::Property >();
        tPropFSVelocity->set_val_function( FSVelocityFunc_MDLFluidBench );

        // create constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMFluid
        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
        tCMFluid->set_property( tPropFluidDensity, "Density" );
        tCMFluid->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_parameters( { {{ 36.0 }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGVelocityBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGPressureBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGInletPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        tIWGInletPressure->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGInletPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletPressure->set_property( tPropInletPressure, "Pressure" );

        std::shared_ptr< fem::IWG > tIWGFSVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGFSVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGFSVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSVelocity->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGFSVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGFSPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGFSPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGFSPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSPressure->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVX->set_output_type( vis::Output_Type::VX );
        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVX->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVY->set_output_type( vis::Output_Type::VY );
        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVY->set_output_type_index( 1 );

        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIP->set_output_type( vis::Output_Type::P );
        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIQIP->set_output_type_index( 0 );

        // create set info
        // --------------------------------------------------------------------------------------
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p10" );
        tSetBulk1.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk1.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p10" );
        tSetBulk2.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk2.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );

        // Fluid/solid bottom
        fem::Set_User_Info tSetFSBottom;
        tSetFSBottom.set_mesh_set_name( "iside_g_0_b0_10_b1_2" );
        tSetFSBottom.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid top
        fem::Set_User_Info tSetFSTop;
        tSetFSTop.set_mesh_set_name( "iside_g_1_b0_10_b1_14" );
        tSetFSTop.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Inlet
        fem::Set_User_Info tSetInlet;
        tSetInlet.set_mesh_set_name( "iside_g_2_b0_10_b1_8" );
        tSetInlet.set_IWGs( { tIWGInletPressure } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 5 );
        tSetInfo( 0 )  = tSetBulk1;
        tSetInfo( 1 )  = tSetBulk2;
        tSetInfo( 2 )  = tSetInlet;
        tSetInfo( 3 )  = tSetFSBottom;
        tSetInfo( 4 )  = tSetFSTop;

        // create model
        // --------------------------------------------------------------------------------------
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // define outputs
        // --------------------------------------------------------------------------------------
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                                 vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                                 "MDL_Fluid_Benchmark_Immersed_Inlet_Pressure_Output.exo",
                                 { "HMR_dummy_c_p10", "HMR_dummy_n_p10" },
                                 { "VX", "VY", "P" },
                                 { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                                 { vis::Output_Type::VX,  vis::Output_Type::VY, vis::Output_Type::P } );
        tModel->set_output_manager( &tOutputData );

        // --------------------------------------------------------------------------------------
        // define linear solver and algorithm
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
        = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

 //       tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_all;
 //       tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
 //       tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // --------------------------------------------------------------------------------------
        // define nonlinear solver and algorithm
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // --------------------------------------------------------------------------------------
        // define time solver and algorithm
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
        = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
        tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );

        tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );

        // --------------------------------------------------------------------------------------
        // solve and check
        tTimeSolver.solve();

//        // create linear solver and algorithm
//        // --------------------------------------------------------------------------------------
//        sol::SOL_Warehouse tSolverWarehouse( tModel->get_solver_interface() );
//
//        moris::Cell< moris::Cell< moris::ParameterList > > tParameterlist( 7 );
//        for( uint Ik = 0; Ik < 7; Ik ++)
//        {
//            tParameterlist( Ik ).resize(1);
//        }
//
//        tParameterlist( 0 )(0) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS2_IMPL );
////        tParameterlist( 0 )(0).set("AZ_diagnostics"    , AZ_none  );
////        tParameterlist( 0 )(0).set("AZ_output"         , AZ_none  );
////        tParameterlist( 0 )(0).set("AZ_max_iter"       , 10000    );
////        tParameterlist( 0 )(0).set("AZ_solver"         , AZ_gmres );
////        tParameterlist( 0 )(0).set("AZ_subdomain_solve", AZ_ilu   );
////        tParameterlist( 0 )(0).set("AZ_graph_fill"     , 10       );
////        tParameterlist( 0 )(0).set("Use_ML_Prec"       ,  true    );
//
//        tParameterlist( 1 )(0) = moris::prm::create_linear_solver_parameter_list();
//        tParameterlist( 2 )(0) = moris::prm::create_nonlinear_algorithm_parameter_list();
//        tParameterlist( 3 )(0) = moris::prm::create_nonlinear_solver_parameter_list();
//        tParameterlist( 3 )(0).set("NLA_DofTypes"      , std::string("VX,VY;P") );
//
//        tParameterlist( 4 )(0) = moris::prm::create_time_solver_algorithm_parameter_list();
//        tParameterlist( 5 )(0) = moris::prm::create_time_solver_parameter_list();
//        tParameterlist( 5 )(0).set("TSA_DofTypes"      , std::string("VX,VY;P") );
//
//        tParameterlist( 6 )(0) = moris::prm::create_solver_warehouse_parameterlist();
//
//        tSolverWarehouse.set_parameterlist( tParameterlist );
//
//        tSolverWarehouse.initialize();
//
//        tsa::Time_Solver * tTimeSolver = tSolverWarehouse.get_main_time_solver();
//
//        tTimeSolver->set_output( 0, tSolverOutputCriteria_MDLFluidBench );
//
//        tTimeSolver->solve();


//         //print full solution for debug
//         Matrix<DDRMat> tFullSol;
//         tTimeSolver->get_full_solution(tFullSol);
//         print(tFullSol,"tFullSol");

        // clean up
        //------------------------------------------------------------------------------
        delete tIntegMesh1;
        delete tModel;
    }
}

//-------------------------------------------------------------------------------------
TEST_CASE("MDL_Fluid_Benchmark_Conform_Inlet_Velocity","[MDL_Fluid_Benchmark_Conform_Inlet_Velocity]")
{
    if( par_size() <= 1 )
    {
        // Geometry Parameters
        moris::real tDomainLX      = 5.0; /* Length of full domain in x (m) */
        moris::real tDomainLY      = 1.0; /* Length of full domain in y (m) */
        moris::real tChannelRadius = 0.5; /* channel radius  (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0; /* Inlet pressure  () */
        moris::real tOutletPressure =  0.0; /* Outlet pressure () */
        moris::real tGammaNitsche   = 1000.0; /* Penalty for Dirichlet BC */

        // Mesh Setup
        moris::uint tNumX   = 200; /* Number of elements in x*/
        moris::uint tNumY   = 40; /* Number of elements in y*/
        moris::uint tNumRef = 0;    /* Number of HMR refinements */
        moris::uint tOrder  = 1;    /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint tLagrangeMeshIndex = 0;

        ParameterList tParameters = prm::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) );
        tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
        tParameters.set( "lagrange_output_meshes", std::string("0") );

        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", 0 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
        tHMR.finalize();

        // construct a mesh manager for the fem
        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIGMesh = tHMR.create_integration_mesh(1, 0, *tIPMesh );

       // place the pair in mesh manager
       mtk::Mesh_Manager tMeshManager;
       tMeshManager.register_mesh_pair( tIPMesh.get(), tIGMesh.get() );

        // create for fem
        // --------------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropInletVelocity = std::make_shared< fem::Property >();
        tPropInletVelocity->set_parameters( { {{ tChannelRadius }}, {{ 0.0 }} } );
        tPropInletVelocity->set_val_function( InletVelocityFunc_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFSVelocity = std::make_shared< fem::Property >();
        tPropFSVelocity->set_val_function( FSVelocityFunc_MDLFluidBench );

        // create constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMFluid
        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
        tCMFluid->set_property( tPropFluidDensity, "Density" );
        tCMFluid->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_parameters( { {{ 36.0 }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGVelocityBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGPressureBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGInletVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGInletVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGInletVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletVelocity->set_property( tPropInletVelocity, "Dirichlet" );
        tIWGInletVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGInletVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGFSVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGFSVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGFSVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSVelocity->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGFSVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGInletPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGInletPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGInletPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletPressure->set_property( tPropInletVelocity, "Dirichlet" );
        tIWGInletPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        std::shared_ptr< fem::IWG > tIWGFSPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGFSPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGFSPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSPressure->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVX->set_output_type( vis::Output_Type::VX );
        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVX->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVY->set_output_type( vis::Output_Type::VY );
        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVY->set_output_type_index( 1 );

        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIP->set_output_type( vis::Output_Type::P );
        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIQIP->set_output_type_index( 0 );

        // create set info
        // --------------------------------------------------------------------------------------
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_set_name( "HMR_dummy" );
        tSetBulk.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );

        // Fluid/solid velocity bottom
        fem::Set_User_Info tSetFSBottom;
        tSetFSBottom.set_mesh_index( 3 );
        tSetFSBottom.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid velocity top
        fem::Set_User_Info tSetFSTop;
        tSetFSTop.set_mesh_index( 1 );
        tSetFSTop.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Inlet velocity
        fem::Set_User_Info tSetInlet;
        tSetInlet.set_mesh_index( 4 );
        tSetInlet.set_IWGs( { tIWGInletVelocity, tIWGInletPressure } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetInlet;
        tSetInfo( 2 ) = tSetFSBottom;
        tSetInfo( 3 ) = tSetFSTop;

        // create model
        // --------------------------------------------------------------------------------------
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // define outputs
        // --------------------------------------------------------------------------------------
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                                 vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                                 "MDL_Fluid_Benchmark_Conform_Inlet_Velocity_Output.exo",
                                 { "HMR_dummy" },
                                 { "VX", "VY", "P" },
                                 { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                                 { vis::Output_Type::VX,  vis::Output_Type::VY, vis::Output_Type::P } );
        tModel->set_output_manager( &tOutputData );

        // --------------------------------------------------------------------------------------
        // define linear solver and algorithm
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
        = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // --------------------------------------------------------------------------------------
        // define nonlinear solver and algorithm
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

        // --------------------------------------------------------------------------------------
        // define time solver and algorithm
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
        = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;

        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

        sol::SOL_Warehouse tSolverWarehouse;

        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
        tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );

        tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );

        // --------------------------------------------------------------------------------------
        // solve and check
        tTimeSolver.solve();

//         //print full solution for debug
//         Matrix<DDRMat> tFullSol;
//         tTimeSolver->get_full_solution(tFullSol);
//         print(tFullSol,"tFullSol");

        // clean up
        //------------------------------------------------------------------------------
//        delete tIntegMesh1;
        delete tModel;
    }
}

//-------------------------------------------------------------------------------------
TEST_CASE("MDL_Fluid_Benchmark_Conform_Inlet_Pressure","[MDL_Fluid_Benchmark_Conform_Inlet_Pressure]")
{
    if( par_size() <= 1 )
    {
        // Geometry Parameters
        moris::real tDomainLX      = 5.0; /* Length of full domain in x (m) */
        moris::real tDomainLY      = 1.0; /* Length of full domain in y (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0; /* Inlet pressure  () */
        moris::real tOutletPressure =  0.0; /* Outlet pressure () */
        moris::real tGammaNitsche   = 1000.0; /* Penalty for Dirichlet BC */

        // Mesh Setup
        moris::uint tNumX   = 200; /* Number of elements in x*/
        moris::uint tNumY   = 40; /* Number of elements in y*/
        moris::uint tNumRef = 0;    /* Number of HMR refinements */
        moris::uint tOrder  = 1;    /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint tLagrangeMeshIndex = 0;

        ParameterList tParameters = prm::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) );
        tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
        tParameters.set( "lagrange_output_meshes", std::string("0") );

        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", 0 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
        tHMR.finalize();

        // construct a mesh manager for the fem
        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIGMesh = tHMR.create_integration_mesh(1, 0, *tIPMesh );

       // place the pair in mesh manager
       mtk::Mesh_Manager tMeshManager;
       tMeshManager.register_mesh_pair( tIPMesh.get(), tIGMesh.get() );

        // create for fem
        // --------------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFSVelocity = std::make_shared< fem::Property >();
        tPropFSVelocity->set_val_function( FSVelocityFunc_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropInletPressure = std::make_shared< fem::Property >();
        tPropInletPressure->set_parameters( { {{ tInletPressure }} } );
        tPropInletPressure->set_val_function( InletPressureFunc_MDLFluidBench );

        // create constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMFluid
        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
        tCMFluid->set_property( tPropFluidDensity, "Density" );
        tCMFluid->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_parameters( { {{ 36.0 }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGVelocityBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGPressureBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::P }, { MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGInletPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        tIWGInletPressure->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGInletPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletPressure->set_property( tPropInletPressure, "Pressure" );

        std::shared_ptr< fem::IWG > tIWGFSVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGFSVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGFSVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSVelocity->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGFSVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGFSPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGFSPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGFSPressure->set_dof_type_list( {{ MSI::Dof_Type::P }, { MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tIWGFSPressure->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVX->set_output_type( vis::Output_Type::VX );
        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVX->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVY->set_output_type( vis::Output_Type::VY );
        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVY->set_output_type_index( 1 );

        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIP->set_output_type( vis::Output_Type::P );
        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIQIP->set_output_type_index( 0 );

        // create set info
        // --------------------------------------------------------------------------------------
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_set_name( "HMR_dummy" );
        tSetBulk.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );

        // Fluid/solid velocity bottom
        fem::Set_User_Info tSetFSBottom;
        tSetFSBottom.set_mesh_index( 3 );
        tSetFSBottom.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid velocity top
        fem::Set_User_Info tSetFSTop;
        tSetFSTop.set_mesh_index( 1 );
        tSetFSTop.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Inlet velocity
        fem::Set_User_Info tSetInlet;
        tSetInlet.set_mesh_index( 4 );
        tSetInlet.set_IWGs( { tIWGInletPressure } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetInlet;
        tSetInfo( 2 ) = tSetFSBottom;
        tSetInfo( 3 ) = tSetFSTop;

        // create model
        // --------------------------------------------------------------------------------------
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // define outputs
        // --------------------------------------------------------------------------------------
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                                 vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                                 "MDL_Fluid_Benchmark_Conform_Inlet_Pressure_Output.exo",
                                 { "HMR_dummy" },
                                 { "VX", "VY", "P" },
                                 { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                                 { vis::Output_Type::VX,  vis::Output_Type::VY, vis::Output_Type::P } );
        tModel->set_output_manager( &tOutputData );

        // --------------------------------------------------------------------------------------
        // define linear solver and algorithm
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
        = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // --------------------------------------------------------------------------------------
        // define nonlinear solver and algorithm
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
        tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );

        // --------------------------------------------------------------------------------------
        // define time solver and algorithm
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
        = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;
        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
        tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );

        sol::SOL_Warehouse tSolverWarehouse;
        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        // --------------------------------------------------------------------------------------
        // solve and check
        tTimeSolver.solve();

//         //print full solution for debug
//         Matrix<DDRMat> tFullSol;
//         tTimeSolver->get_full_solution(tFullSol);
//         print(tFullSol,"tFullSol");

        // clean up
        //------------------------------------------------------------------------------
//        delete tIntegMesh1;
        delete tModel;
    }
}

//-------------------------------------------------------------------------------------

TEST_CASE("MDL_Fluid_Benchmark_Conform_Inlet_Pressure_3D","[MDL_Fluid_Benchmark_Conform_Inlet_Pressure_3D]")
{
    if( par_size() <= 1 )
    {
        // Geometry Parameters
        moris::real tDomainLX      = 5.0; /* Length of full domain in x (m) */
        moris::real tDomainLY      = 1.0; /* Length of full domain in y (m) */
        moris::real tDomainLZ      = 1.0; /* Length of full domain in z (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tInletPressure  = 20.0; /* Inlet pressure  () */
        moris::real tOutletPressure =  0.0; /* Outlet pressure () */
        moris::real tGammaNitsche   = 1000.0; /* Penalty for Dirichlet BC */

        // Mesh Setup
        moris::uint tNumX   = 50; /* Number of elements in x*/
        moris::uint tNumY   = 10; /* Number of elements in y*/
        moris::uint tNumZ   = 10; /* Number of elements in z*/
        moris::uint tNumRef = 0;    /* Number of HMR refinements */
        moris::uint tOrder  = 1;    /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint tLagrangeMeshIndex = 0;

        ParameterList tParameters = prm::create_hmr_parameter_list();
        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY) + "," + std::to_string(tNumZ));
        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) + "," + std::to_string(tDomainLZ) );
        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2) + "," + std::to_string(-tDomainLY/2) + "," + std::to_string(-tDomainLZ/2) );
        tParameters.set( "domain_sidesets", std::string("1,2,3,4,5,6") );
        tParameters.set( "lagrange_output_meshes", std::string("0") );

        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", 0 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
        tHMR.finalize();

        // construct a mesh manager for the fem
        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tIPMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIGMesh = tHMR.create_integration_mesh(1, 0, *tIPMesh );

       // place the pair in mesh manager
       mtk::Mesh_Manager tMeshManager;
       tMeshManager.register_mesh_pair( tIPMesh.get(), tIGMesh.get() );

        // create for fem
        // --------------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFSVelocity = std::make_shared< fem::Property >();
        tPropFSVelocity->set_val_function( FSVelocityFunc_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropInletPressure = std::make_shared< fem::Property >();
        tPropInletPressure->set_parameters( { {{ tInletPressure }} } );
        tPropInletPressure->set_val_function( InletPressureFunc_MDLFluidBench );

        // create constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMFluid
        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }} );
        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
        tCMFluid->set_property( tPropFluidDensity, "Density" );
        tCMFluid->set_space_dim( 3 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_parameters( { {{ 36.0 }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }}, mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGVelocityBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGPressureBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::P }, { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }}, mtk::Master_Slave::MASTER );
        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGInletPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
        tIWGInletPressure->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGInletPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGInletPressure->set_property( tPropInletPressure, "Pressure" );

        std::shared_ptr< fem::IWG > tIWGFSVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGFSVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGFSVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGFSVelocity->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGFSVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGFSPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGFSPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGFSPressure->set_dof_type_list( {{ MSI::Dof_Type::P }, { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }}, mtk::Master_Slave::MASTER );
        tIWGFSPressure->set_property( tPropFSVelocity, "Dirichlet" );
        tIWGFSPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVX->set_output_type( vis::Output_Type::VX );
        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } }, mtk::Master_Slave::MASTER );
        tIQIVX->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVY->set_output_type( vis::Output_Type::VY );
        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } }, mtk::Master_Slave::MASTER );
        tIQIVY->set_output_type_index( 1 );

        std::shared_ptr< fem::IQI > tIQIVZ = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVZ->set_output_type( vis::Output_Type::VZ );
        tIQIVZ->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } }, mtk::Master_Slave::MASTER );
        tIQIVZ->set_output_type_index( 2 );

        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIP->set_output_type( vis::Output_Type::P );
        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIQIP->set_output_type_index( 0 );

        // create set info
        // --------------------------------------------------------------------------------------
        fem::Set_User_Info tSetBulk;
        tSetBulk.set_mesh_set_name( "HMR_dummy" );
        tSetBulk.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk.set_IQIs( { tIQIVX, tIQIVY, tIQIVZ, tIQIP } );

        // Fluid/solid velocity bottom
        fem::Set_User_Info tSetFSBottom;
        tSetFSBottom.set_mesh_index( 3 );
        tSetFSBottom.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid velocity top
        fem::Set_User_Info tSetFSTop;
        tSetFSTop.set_mesh_index( 1 );
        tSetFSTop.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid velocity top
        fem::Set_User_Info tSetFSBack;
        tSetFSBack.set_mesh_index( 5 );
        tSetFSBack.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Fluid/solid velocity top
        fem::Set_User_Info tSetFSFront;
        tSetFSFront.set_mesh_index( 6 );
        tSetFSFront.set_IWGs( { tIWGFSVelocity, tIWGFSPressure } );

        // Inlet velocity
        fem::Set_User_Info tSetInlet;
        tSetInlet.set_mesh_index( 4 );
        tSetInlet.set_IWGs( { tIWGInletPressure } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 6 );
        tSetInfo( 0 ) = tSetBulk;
        tSetInfo( 1 ) = tSetInlet;
        tSetInfo( 2 ) = tSetFSBottom;
        tSetInfo( 3 ) = tSetFSTop;
        tSetInfo( 4 ) = tSetFSBack;
        tSetInfo( 5 ) = tSetFSFront;

        // create model
        // --------------------------------------------------------------------------------------
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // define outputs
        // --------------------------------------------------------------------------------------
        vis::Output_Manager tOutputData;
        tOutputData.set_outputs( 0,
                                 vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                                 "MDL_Fluid_Benchmark_Conform_Inlet_Pressure_3D_Output.exo",
                                 { "HMR_dummy" },
                                 { "VX", "VY", "VZ", "P" },
                                 { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                                 { vis::Output_Type::VX, vis::Output_Type::VY, vis::Output_Type::VZ, vis::Output_Type::P } );
        tModel->set_output_manager( &tOutputData );

        // --------------------------------------------------------------------------------------
        // define linear solver and algorithm
        dla::Solver_Factory  tSolFactory;
        std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
        = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

        dla::Linear_Solver tLinSolver;

        tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

        // --------------------------------------------------------------------------------------
        // define nonlinear solver and algorithm
        NLA::Nonlinear_Solver_Factory tNonlinFactory;
        std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
        = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

        tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

        NLA::Nonlinear_Solver tNonlinearSolver;
        tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
        tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ, MSI::Dof_Type::P } );

        // --------------------------------------------------------------------------------------
        // define time solver and algorithm
        tsa::Time_Solver_Factory tTimeSolverFactory;
        std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
        = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
        tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

        tsa::Time_Solver tTimeSolver;
        tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
        tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ, MSI::Dof_Type::P } );
        tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );

        sol::SOL_Warehouse tSolverWarehouse;
        tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
        tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
        tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

        // --------------------------------------------------------------------------------------
        // solve and check
        tTimeSolver.solve();

//         //print full solution for debug
//         Matrix<DDRMat> tFullSol;
//         tTimeSolver->get_full_solution(tFullSol);
//         print(tFullSol,"tFullSol");

        // clean up
        //------------------------------------------------------------------------------
//        delete tIntegMesh1;
        delete tModel;
    }
}

void RVelocityFunc_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // velocity magnitude
    real tVelocity = aParameters( 0 )( 0 );

    //
    Matrix< DDRMat > tNormal;
    aFIManager->get_IG_geometry_interpolator()->get_normal( tNormal );

    Matrix< DDRMat > tTangent( 2, 1, 0.0 );
    tTangent( 0 ) =   tNormal( 1 );
    tTangent( 1 ) = - tNormal( 0 );

    // set size for aPropMatrix
    aPropMatrix = tVelocity * tTangent;
}

//TEST_CASE("MDL_Fluid_Benchmark_Tube","[MDL_Fluid_Benchmark_Tube]")
//{
//    if(par_size()<=1)
//    {
//        // Geometry Parameters
//        moris::real tDomainLX = 10.0;                   /* Length of full domain in x (m) */
//        moris::real tDomainLY = 10.0;                   /* Length of full domain in y (m) */
//        Matrix<DDRMat> tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
//        moris::real tRInner = 0.55;                     /* Inner circle radius (m) */
//        moris::real tROuter = 1.05;                     /* Outer circle radius (m) */
//
//        //Material Parameters
//        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
//        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */
//
//        // Boundary Conditions
//        moris::real tRInVelocity  = 20.0; /* Inlet pressure  () */
//        moris::real tROutVelocity =  0.0; /* Outlet pressure () */
//        moris::real tGammaNitsche   = 1000.0;  /* Penalty for Dirichlet BC */
//
//        // Mesh Setup
//        moris::uint tNumX   = 20; /* Number of elements in x*/
//        moris::uint tNumY   = 20; /* Number of elements in y*/
//        moris::uint tNumRef = 2;  /* Number of HMR refinements */
//        moris::uint tOrder  = 1;  /* Lagrange Order and Bspline Order (forced to be same for this example) */
//
//        uint tLagrangeMeshIndex = 0;
//        std::string tOuterFieldName   = "Outercircle";
//        std::string tInnerFieldName   = "Innercircle";
//        ParameterList tParameters = prm::create_hmr_parameter_list();
//
//        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
//        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
//        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2+tCenterPoint(0)) + "," + std::to_string(-tDomainLY/2+tCenterPoint(1)) );
//        tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
//        tParameters.set( "lagrange_output_meshes", std::string("0") );
//
//        tParameters.set( "lagrange_orders", std::string("1") );
//        tParameters.set( "lagrange_pattern", std::string("0") );
//        tParameters.set( "bspline_orders", std::string("1") );
//        tParameters.set( "bspline_pattern", std::string("0") );
//
//        tParameters.set( "lagrange_to_bspline", std::string("0") );
//
//        tParameters.set( "truncate_bsplines", 1 );
//        tParameters.set( "refinement_buffer", 3 );
//        tParameters.set( "staircase_buffer", 3 );
//        tParameters.set( "initial_refinement", 2 );
//
//        tParameters.set( "use_multigrid", 0 );
//        tParameters.set( "severity_level", 2 );
//        tParameters.set( "use_number_aura", 0 );
//
//        hmr::HMR tHMR( tParameters );
//
//        //initial refinement
//        tHMR.perform_initial_refinement( 0 );
//
//        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//
//        //  create field
//        std::shared_ptr< moris::hmr::Field > tOuterField = tMesh->create_field( tOuterFieldName, tLagrangeMeshIndex );
//        std::shared_ptr< moris::hmr::Field > tInnerField = tMesh->create_field( tInnerFieldName, tLagrangeMeshIndex );
//
//        for( uint k=0; k<tNumRef; ++k )
//        {
//            moris::ge::Circle tCircle0( tROuter, tCenterPoint( 0 ), tCenterPoint( 1 ) );
//            moris::ge::Circle tCircle1( tRInner, tCenterPoint( 0 ), tCenterPoint( 1 ) );
//            moris::Cell<moris::ge::GEN_Geometry*> tGeomVec = { &tCircle0, &tCircle1 };
//
//            moris::ge::GEN_Phase_Table     tPhaseTable( tGeomVec.size(),  Phase_Table_Structure::EXP_BASE_2 );
//            moris::ge::GEN_Geometry_Engine tGENGeometryEngine( tGeomVec, tPhaseTable,2 );
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
//        tHMR.finalize();
////        tHMR.save_to_exodus( 0, tHMRIPMeshFileName );
//
//        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);
//
//        //-----------------------------------------------------------------------------------------------
//
//        moris::ge::Circle tCircle2( tROuter, tCenterPoint( 0 ), tCenterPoint( 1 ) );
//        moris::ge::Circle tCircle3( tRInner, tCenterPoint( 0 ), tCenterPoint( 1 ) );
//
//        // NOTE the order of this geometry vector is important. If it changes the resulting bulk phase of the output mesh change.
//        moris::Cell<moris::ge::GEN_Geometry*> tGeomVec0 = { &tCircle2, &tCircle3 };
//
//          size_t tModelDimension = 2;
//          moris::ge::GEN_Phase_Table         tPhaseTable0( tGeomVec0.size(),  Phase_Table_Structure::EXP_BASE_2 );
//          moris::ge::GEN_Geometry_Engine     tGENGeometryEngine0( tGeomVec0, tPhaseTable0, tModelDimension );
//
//          // --------------------------------------------------------------------------------------
//          xtk::Model tXTKModel(tModelDimension,tInterpolationMesh.get(),&tGENGeometryEngine0);
//          tXTKModel.mVerbose = true;
//
//        //Specify decomposition Method and Cut Mesh ---------------------------------------
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//        tXTKModel.decompose(tDecompositionMethods);
//
//        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
////        tXTKModel.construct_face_oriented_ghost_penalization_cells();
//
//        xtk::Output_Options tOutputOptions;
//        tOutputOptions.mAddNodeSets = false;
//        tOutputOptions.mAddSideSets = true;
//        tOutputOptions.mAddClusters = false;
//
//        // output integration mesh
//        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
//        std::string tOutputFile = "./mdl_exo/MDL_Fluid_Benchmark_Tube.exo";
//        tIntegMesh1->create_output_mesh(tOutputFile);
//
//        // get meshes for FEM
//        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();
//
//        // place the pair in mesh manager
//        mtk::Mesh_Manager tMeshManager;
//        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);
//
//        // create for fem
//        // --------------------------------------------------------------------------------------
//        // create the properties
//        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
//        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
//        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );
//
//        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
//        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
//        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );
//
//        std::shared_ptr< fem::Property > tPropRInVelocity = std::make_shared< fem::Property >();
//        tPropRInVelocity->set_parameters( { {{ tRInVelocity }} } );
//        tPropRInVelocity->set_val_function( RVelocityFunc_MDLFluidBench );
//
//        std::shared_ptr< fem::Property > tPropROutVelocity = std::make_shared< fem::Property >();
//        tPropROutVelocity->set_parameters( { {{ tROutVelocity }} } );
//        tPropROutVelocity->set_val_function( RVelocityFunc_MDLFluidBench );
//
//        // create constitutive models
//        fem::CM_Factory tCMFactory;
//
//        std::shared_ptr< fem::Constitutive_Model > tCMFluid
//        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
//        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
//        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
//        tCMFluid->set_property( tPropFluidDensity, "Density" );
//        tCMFluid->set_space_dim( 2 );
//
//        // define stabilization parameters
//        fem::SP_Factory tSPFactory;
//
//        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
//        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
//        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
//        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
//        tSPIncFlow->set_parameters( { {{ 36.0 }} } );
//
//        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
//        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
//        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
//        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
//        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
//        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );
//
//        // define the IWGs
//        fem::IWG_Factory tIWGFactory;
//
//        std::shared_ptr< fem::IWG > tIWGVelocityBulk
//        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
//        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
//        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
//        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
//        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );
//
//        std::shared_ptr< fem::IWG > tIWGPressureBulk
//        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
//        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
//        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
//        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
//        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );
//
//        std::shared_ptr< fem::IWG > tIWGRInVelocity
//        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
//        tIWGRInVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
//        tIWGRInVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tIWGRInVelocity->set_property( tPropRInVelocity, "Dirichlet" );
//        tIWGRInVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
//        tIWGRInVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );
//
//        std::shared_ptr< fem::IWG > tIWGRInPressure
//        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
//        tIWGRInPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
//        tIWGRInPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tIWGRInPressure->set_property( tPropRInVelocity, "Dirichlet" );
//        tIWGRInPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
//
//        std::shared_ptr< fem::IWG > tIWGROutVelocity
//        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
//        tIWGROutVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
//        tIWGROutVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tIWGROutVelocity->set_property( tPropROutVelocity, "Dirichlet" );
//        tIWGROutVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
//        tIWGROutVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );
//
//        std::shared_ptr< fem::IWG > tIWGROutPressure
//        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
//        tIWGROutPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
//        tIWGROutPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
//        tIWGROutPressure->set_property( tPropROutVelocity, "Dirichlet" );
//        tIWGROutPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
//
//        // create the IQIs
//        // --------------------------------------------------------------------------------------
//        fem::IQI_Factory tIQIFactory;
//
//        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
//        tIQIVX->set_output_type( vis::Output_Type::VX );
//        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
//        tIQIVX->set_output_type_index( 0 );
//
//        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
//        tIQIVY->set_output_type( vis::Output_Type::VY );
//        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
//        tIQIVY->set_output_type_index( 1 );
//
//        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
//        tIQIP->set_output_type( vis::Output_Type::P );
//        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
//        tIQIP->set_output_type_index( 0 );
//
//        // create set info
//        // --------------------------------------------------------------------------------------
//        fem::Set_User_Info tSetBulk1;
//        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p1" );
//        tSetBulk1.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
//        tSetBulk1.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );
//
//        fem::Set_User_Info tSetBulk2;
//        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p1" );
//        tSetBulk2.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
//        tSetBulk2.set_IQIs( { tIQIVX, tIQIVY, tIQIP } );
//
//        // At RIn
//        fem::Set_User_Info tSetRIn;
//        tSetRIn.set_mesh_set_name( "iside_g_1_b0_1_b1_0" );
//        tSetRIn.set_IWGs( { tIWGRInVelocity, tIWGRInPressure } );
//
//        // At ROut
//        fem::Set_User_Info tSetROut;
//        tSetROut.set_mesh_set_name( "iside_g_0_b0_1_b1_3" );
//        tSetROut.set_IWGs( { tIWGROutVelocity, tIWGROutPressure } );
//
//
//        // create a cell of set info
//        moris::Cell< fem::Set_User_Info > tSetInfo( 4 );
//        tSetInfo( 0 )  = tSetBulk1;
//        tSetInfo( 1 )  = tSetBulk2;
//        tSetInfo( 2 )  = tSetRIn;
//        tSetInfo( 3 )  = tSetROut;
//
//        // create model
//        // --------------------------------------------------------------------------------------
//        mdl::Model * tModel = new mdl::Model( &tMeshManager,
//                                              0,
//                                              tSetInfo,
//                                              0, false );
//
//        // define outputs
//       // --------------------------------------------------------------------------------------
//       vis::Output_Manager tOutputData;
//       tOutputData.set_outputs( 0,
//                                vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
//                                "MDL_Fluid_Benchmark_Tube_Output.exo",
//                                { "HMR_dummy_c_p1", "HMR_dummy_n_p1" },
//                                { "VX", "VY", "P" },
//                                { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
//                                { vis::Output_Type::VX,  vis::Output_Type::VY, vis::Output_Type::P } );
//       tModel->set_output_manager( &tOutputData );
//
//       // --------------------------------------------------------------------------------------
//       // define linear solver and algorithm
//       dla::Solver_Factory  tSolFactory;
//       std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
//       = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );
//
////       tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_all;
////       tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
////       tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;
//
//       dla::Linear_Solver tLinSolver;
//
//       tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );
//
//       // --------------------------------------------------------------------------------------
//       // define nonlinear solver and algorithm
//       NLA::Nonlinear_Solver_Factory tNonlinFactory;
//       std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
//       = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );
//
//       tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );
//
//       NLA::Nonlinear_Solver tNonlinearSolver;
//       tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );
//
//       // --------------------------------------------------------------------------------------
//       // define time solver and algorithm
//       tsa::Time_Solver_Factory tTimeSolverFactory;
//       std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
//       = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );
//
//       tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );
//
//       tsa::Time_Solver tTimeSolver;
//
//       tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );
//
//       sol::SOL_Warehouse tSolverWarehouse;
//
//       tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());
//
//       tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
//       tTimeSolver.set_solver_warehouse( &tSolverWarehouse );
//
//       tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
//       tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
//
//       tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );
//
//       // --------------------------------------------------------------------------------------
//       // solve and check
//       tTimeSolver.solve();
//
////         //print full solution for debug
////         Matrix<DDRMat> tFullSol;
////         tTimeSolver.get_full_solution(tFullSol);
////         print(tFullSol,"tFullSol");
//
//        // clean up
//        //------------------------------------------------------------------------------
//        delete tIntegMesh1;
//        delete tModel;
//    }
//}

void ImposedVelocityFunc_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // velocity magnitude
    real tA = aParameters( 0 )( 0 );
    real tB = aParameters( 1 )( 0 );

    // radius
    real tx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
    real ty = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
    real tR = std::pow( std::pow( tx, 2.0 ) + std::pow( ty, 2.0 ), 0.5 );

    // radial velocity
    real tVTheta = tA * tR + tB / tR;

    // set size for aPropMatrix
    aPropMatrix.set_size( 2, 1, 0.0 );
    aPropMatrix( 0 ) = - tVTheta * ty / tR;
    aPropMatrix( 1 ) =   tVTheta * tx / tR;
}

void AnalyticdVelocitydxFunc_MDLFluidBench
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    // velocity magnitude
    real tA = aParameters( 0 )( 0 );
    real tB = aParameters( 1 )( 0 );

    // radius
    real tx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
    real ty = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
    real tR = std::pow( std::pow( tx, 2.0 ) + std::pow( ty, 2.0 ), 0.5 );

    // set size for aPropMatrix
    aPropMatrix.set_size( 2, 2, 0.0 );
    aPropMatrix( 0, 0 ) = 2.0 * tB * tx * ty / std::pow( tR, 4.0 );
    aPropMatrix( 1, 0 ) = - tA - tB / std::pow( tR, 2.0 ) + 2.0 * tB * std::pow( ty, 2.0 )/ std::pow( tR, 4.0 );

    aPropMatrix( 0, 1 ) = tA + tB / std::pow( tR, 2.0 ) - 2.0 * tB * std::pow( tx, 2.0 ) / std::pow( tR, 4.0);
    aPropMatrix( 1, 1 ) = - 2.0 * tB * tx * ty / std::pow( tR, 2.0 );
}


TEST_CASE("MDL_Fluid_Benchmark_Radial_Couette_Flow","[MDL_Fluid_Benchmark_Radial_Couette_Flow]")
{
    if(par_size()<=1)
    {
        // Geometry Parameters
        moris::real tDomainLX = 10.0;                   /* Length of full domain in x (m) */
        moris::real tDomainLY = 10.0;                   /* Length of full domain in y (m) */
        Matrix<DDRMat> tShift = { { 0.01, 0.01 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
        Matrix<DDRMat> tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the block (intentionally off 0.0,0.0 to prevent interface at node)*/
        moris::real tRIn = 1.0;                     /* Inner circle radius (m) */
        moris::real tROut = 2.0;                     /* Outer circle radius (m) */

        //Material Parameters
        moris::real tFluidDensity   = 1.0; /* Fluid density   () */
        moris::real tFluidViscosity = 1.0; /* Fluid viscosity () */

        // Boundary Conditions
        moris::real tRInVelocity  = 20.0; /* Inlet pressure  () */
        moris::real tROutVelocity =  0.0; /* Outlet pressure () */
        moris::real tGammaNitsche   = 1000.0;  /* Penalty for Dirichlet BC */
        moris::real tGammaGPmu   = 0.05;  /* Penalty for ghost viscosity */
        moris::real tGammaGPu    = 0.05;  /* Penalty for ghost velocity */
        moris::real tGammaGPp    = 0.5;  /* Penalty for ghost pressure */

        moris::real tOmegaIn  = -5.0;
        moris::real tOmegaOut =  5.0;
        moris::real tKn       =  0.0;
        moris::real tb = 1.0 / tRIn + 2.0 * tKn / std::pow( tRIn, 2.0 );
        moris::real td = 1.0 / tROut - 2.0 * tKn / std::pow( tROut, 2.0 );
        moris::real tB = ( tOmegaOut * tROut - tOmegaIn * tRIn * tROut / tRIn ) / ( td -tb * tROut / tRIn );
        moris::real tA = ( tOmegaIn * tRIn - tb * tB ) / tRIn;

        // Mesh Setup
        moris::uint tNumX   = 20; /* Number of elements in x*/
        moris::uint tNumY   = 20; /* Number of elements in y*/
        moris::uint tNumRef = 0;  /* Number of HMR refinements */
        moris::uint tOrder  = 1;  /* Lagrange Order and Bspline Order (forced to be same for this example) */

        uint tLagrangeMeshIndex = 0;
        std::string tOuterFieldName   = "Outercircle";
        std::string tInnerFieldName   = "Innercircle";
        ParameterList tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", std::to_string(tNumX) + "," + std::to_string(tNumY));
        tParameters.set( "domain_dimensions", std::to_string(tDomainLX) + "," + std::to_string(tDomainLY) );
        tParameters.set( "domain_offset", std::to_string(-tDomainLX/2+tShift(0)) + "," + std::to_string(-tDomainLY/2+tShift(1)) );
        tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
        tParameters.set( "lagrange_output_meshes", std::string("0") );

        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 3 );
        tParameters.set( "staircase_buffer", 3 );
        tParameters.set( "initial_refinement", 2 );

        tParameters.set( "use_multigrid", 0 );
        tParameters.set( "severity_level", 2 );
        tParameters.set( "use_number_aura", 0 );

        hmr::HMR tHMR( tParameters );

        //initial refinement
        tHMR.perform_initial_refinement( 0 );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        //  create field
        std::shared_ptr< moris::hmr::Field > tOuterField = tMesh->create_field( tOuterFieldName, tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Field > tInnerField = tMesh->create_field( tInnerFieldName, tLagrangeMeshIndex );

        for( uint k=0; k<tNumRef; ++k )
        {
            moris::ge::Circle tCircle0( tROut, tCenterPoint( 0 ), tCenterPoint( 1 ) );
            moris::ge::Circle tCircle1( tRIn, tCenterPoint( 0 ), tCenterPoint( 1 ) );
            moris::Cell<moris::ge::GEN_Geometry*> tGeomVec = { &tCircle0, &tCircle1 };

            moris::ge::GEN_Phase_Table     tPhaseTable( tGeomVec.size(),  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine tGENGeometryEngine( tGeomVec, tPhaseTable,2 );

            moris_index tMeshIndex = tGENGeometryEngine.register_mesh( tMesh );

            uint tNumIPNodes = tMesh->get_num_nodes();
            Matrix<DDRMat> tFieldData( tNumIPNodes,1 );
            Matrix<DDRMat> tFieldData0( tNumIPNodes,1 );

            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumIPNodes );
            Matrix< DDRMat > tCoords( tNumIPNodes, 2 );
            for( uint i = 0; i < tNumIPNodes; i++ )
            {
                tCoords.set_row( i, tMesh->get_mtk_vertex(i).get_coords() );
            }

            tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );

            for(uint i=0; i<tNumIPNodes; i++)
            {
                tFieldData( i )  = tGENGeometryEngine.get_entity_phase_val( i, 0 );
                tFieldData0( i ) = tGENGeometryEngine.get_entity_phase_val( i, 1 );
            }

            tHMR.based_on_field_put_elements_on_queue( tFieldData, tLagrangeMeshIndex );
            tHMR.based_on_field_put_elements_on_queue( tFieldData0, tLagrangeMeshIndex );

            tHMR.perform_refinement_based_on_working_pattern( 0, false );
        }
        tHMR.finalize();
//        tHMR.save_to_exodus( 0, tHMRIPMeshFileName );

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

        //-----------------------------------------------------------------------------------------------

        moris::ge::Circle tCircle2( tROut, tCenterPoint( 0 ), tCenterPoint( 1 ) );
        moris::ge::Circle tCircle3( tRIn, tCenterPoint( 0 ), tCenterPoint( 1 ) );

        // NOTE the order of this geometry vector is important. If it changes the resulting bulk phase of the output mesh change.
        moris::Cell<moris::ge::GEN_Geometry*> tGeomVec0 = { &tCircle2, &tCircle3 };

          size_t tModelDimension = 2;
          moris::ge::GEN_Phase_Table         tPhaseTable0( tGeomVec0.size(),  Phase_Table_Structure::EXP_BASE_2 );
          moris::ge::GEN_Geometry_Engine     tGENGeometryEngine0( tGeomVec0, tPhaseTable0, tModelDimension );

          // --------------------------------------------------------------------------------------
          xtk::Model tXTKModel(tModelDimension,tInterpolationMesh.get(),&tGENGeometryEngine0);
          tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
        tXTKModel.construct_face_oriented_ghost_penalization_cells();

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = true;
        tOutputOptions.mAddClusters = false;

        // output integration mesh
        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
        std::string tOutputFile = "./mdl_exo/MDL_Fluid_Benchmark_Radial_Couette.exo";
        tIntegMesh1->create_output_mesh(tOutputFile);

        // get meshes for FEM
        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

        // create for fem
        // --------------------------------------------------------------------------------------
        // create the properties
        std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
        tPropFluidDensity->set_parameters( { {{ tFluidDensity }} } );
        tPropFluidDensity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
        tPropFluidViscosity->set_parameters( { {{ tFluidViscosity }} } );
        tPropFluidViscosity->set_val_function( ConstFuncVal_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropImposedVelocity = std::make_shared< fem::Property >();
        tPropImposedVelocity->set_parameters( { {{ tA }}, {{ tB }} } );
        tPropImposedVelocity->set_val_function( ImposedVelocityFunc_MDLFluidBench );

        std::shared_ptr< fem::Property > tPropAnalyticdVelocitydx = std::make_shared< fem::Property >();
        tPropAnalyticdVelocitydx->set_parameters( { {{ tA }}, {{ tB }} } );
        tPropAnalyticdVelocitydx->set_val_function( AnalyticdVelocitydxFunc_MDLFluidBench );

        // create constitutive models
        fem::CM_Factory tCMFactory;

        std::shared_ptr< fem::Constitutive_Model > tCMFluid
        = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCMFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
        tCMFluid->set_property( tPropFluidViscosity, "Viscosity" );
        tCMFluid->set_property( tPropFluidDensity, "Density" );
        tCMFluid->set_space_dim( 2 );

        // define stabilization parameters
        fem::SP_Factory tSPFactory;

        std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow
        = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
        tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPIncFlow->set_parameters( { {{ 36.0 }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche
        = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
        tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
        tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPNitsche->set_parameters( { {{ tGammaNitsche }} } );

        std::shared_ptr< fem::Stabilization_Parameter > tSPViscousGhost
        = tSPFactory.create_SP( fem::Stabilization_Type::VISCOUS_GHOST );
        tSPViscousGhost->set_parameters( {{{ tGammaGPmu }} });
        tSPViscousGhost->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPConvectiveGhost
        = tSPFactory.create_SP( fem::Stabilization_Type::CONVECTIVE_GHOST );
        tSPConvectiveGhost->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPConvectiveGhost->set_parameters( {{{ tGammaGPu }} });
        tSPConvectiveGhost->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::Stabilization_Parameter > tSPPressureGhost
        = tSPFactory.create_SP( fem::Stabilization_Type::PRESSURE_GHOST );
        tSPPressureGhost->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
        tSPPressureGhost->set_parameters( { {{ tGammaGPp }}, {{ 1.0 }} });
        tSPPressureGhost->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
        tSPPressureGhost->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        std::shared_ptr< fem::IWG > tIWGVelocityBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
        tIWGVelocityBulk->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGVelocityBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGVelocityBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGVelocityBulk->set_property( tPropFluidDensity, "Density" );
        tIWGVelocityBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGPressureBulk
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
        tIWGPressureBulk->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGPressureBulk->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGPressureBulk->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGPressureBulk->set_property( tPropFluidDensity, "Density" );
        tIWGPressureBulk->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

        std::shared_ptr< fem::IWG > tIWGDirichletVelocity
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE );
        tIWGDirichletVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
        tIWGDirichletVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGDirichletVelocity->set_property( tPropImposedVelocity, "Dirichlet" );
        tIWGDirichletVelocity->set_constitutive_model( tCMFluid, "IncompressibleFluid" );
        tIWGDirichletVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

        std::shared_ptr< fem::IWG > tIWGDirichletPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_NITSCHE );
        tIWGDirichletPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGDirichletPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
        tIWGDirichletPressure->set_property( tPropImposedVelocity, "Dirichlet" );
        tIWGDirichletPressure->set_constitutive_model( tCMFluid, "IncompressibleFluid" );

        std::shared_ptr< fem::IWG > tIWGGPViscous
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VISCOUS_VELOCITY_GHOST );
        tIWGGPViscous->set_residual_dof_type( { MSI::Dof_Type::VX, MSI::Dof_Type::VY } );
        tIWGGPViscous->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIWGGPViscous->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::SLAVE );
        tIWGGPViscous->set_stabilization_parameter( tSPViscousGhost, "ViscousGhost" );

        std::shared_ptr< fem::IWG > tIWGGPConvective
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST );
        tIWGGPConvective->set_residual_dof_type( { MSI::Dof_Type::VX, MSI::Dof_Type::VY } );
        tIWGGPConvective->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIWGGPConvective->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::SLAVE );
        tIWGGPConvective->set_stabilization_parameter( tSPConvectiveGhost, "ConvectiveGhost" );

        std::shared_ptr< fem::IWG > tIWGGPPressure
        = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_GHOST );
        tIWGGPPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
        tIWGGPPressure->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIWGGPPressure->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::SLAVE );
        tIWGGPPressure->set_stabilization_parameter( tSPPressureGhost, "PressureGhost" );

        // create the IQIs
        // --------------------------------------------------------------------------------------
        fem::IQI_Factory tIQIFactory;

        std::shared_ptr< fem::IQI > tIQIVX = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVX->set_output_type( vis::Output_Type::VX );
        tIQIVX->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVX->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIVY = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIVY->set_output_type( vis::Output_Type::VY );
        tIQIVY->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIVY->set_output_type_index( 1 );

        std::shared_ptr< fem::IQI > tIQIP = tIQIFactory.create_IQI( fem::IQI_Type::DOF );
        tIQIP->set_output_type( vis::Output_Type::P );
        tIQIP->set_dof_type_list( { { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
        tIQIP->set_output_type_index( 0 );

        std::shared_ptr< fem::IQI > tIQIL2 = tIQIFactory.create_IQI( fem::IQI_Type::L2_ERROR_ANALYTIC );
        tIQIL2->set_output_type( vis::Output_Type::L2_ERROR_ANALYTIC );
        tIQIL2->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIL2->set_property( tPropImposedVelocity, "L2Check", mtk::Master_Slave::MASTER );

        std::shared_ptr< fem::IQI > tIQIH1 = tIQIFactory.create_IQI( fem::IQI_Type::H1_ERROR_ANALYTIC );
        tIQIH1->set_output_type( vis::Output_Type::H1_ERROR_ANALYTIC );
        tIQIH1->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
        tIQIH1->set_property( tPropAnalyticdVelocitydx, "H1Check", mtk::Master_Slave::MASTER );

        // create set info
        // --------------------------------------------------------------------------------------
        fem::Set_User_Info tSetBulk1;
        tSetBulk1.set_mesh_set_name( "HMR_dummy_c_p1" );
        tSetBulk1.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk1.set_IQIs( { tIQIVX, tIQIVY, tIQIP, tIQIL2, tIQIH1 } );

        fem::Set_User_Info tSetBulk2;
        tSetBulk2.set_mesh_set_name( "HMR_dummy_n_p1" );
        tSetBulk2.set_IWGs( { tIWGVelocityBulk, tIWGPressureBulk } );
        tSetBulk2.set_IQIs( { tIQIVX, tIQIVY, tIQIP, tIQIL2, tIQIH1 } );

        // At RIn
        fem::Set_User_Info tSetRIn;
        tSetRIn.set_mesh_set_name( "iside_g_1_b0_1_b1_0" );
        tSetRIn.set_IWGs( { tIWGDirichletVelocity, tIWGDirichletPressure } );

        // At ROut
        fem::Set_User_Info tSetROut;
        tSetROut.set_mesh_set_name( "iside_g_0_b0_1_b1_3" );
        tSetROut.set_IWGs( { tIWGDirichletVelocity, tIWGDirichletPressure } );

        // At outer interface
        fem::Set_User_Info tSetGPROut;
        tSetGPROut.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 1, 3 ) );
        tSetGPROut.set_IWGs( { tIWGGPViscous, tIWGGPConvective, tIWGGPPressure } );

        // At inner interface
        fem::Set_User_Info tSetGPRIn;
        tSetGPRIn.set_mesh_set_name( tEnrIntegMesh.get_dbl_interface_side_set_name( 0, 1 ) );
        tSetGPRIn.set_IWGs( { tIWGGPViscous, tIWGGPConvective, tIWGGPPressure } );

        // create a cell of set info
        moris::Cell< fem::Set_User_Info > tSetInfo( 6 );
        tSetInfo( 0 )  = tSetBulk1;
        tSetInfo( 1 )  = tSetBulk2;
        tSetInfo( 2 )  = tSetRIn;
        tSetInfo( 3 )  = tSetROut;
        tSetInfo( 4 )  = tSetGPROut;
        tSetInfo( 5 )  = tSetGPRIn;

        // create model
        // --------------------------------------------------------------------------------------
        mdl::Model * tModel = new mdl::Model( &tMeshManager,
                                              0,
                                              tSetInfo,
                                              0, false );

        // define outputs
       // --------------------------------------------------------------------------------------
       vis::Output_Manager tOutputData;
       tOutputData.set_outputs( 0,
                                vis::VIS_Mesh_Type::STANDARD, //OVERLAPPING_INTERFACE
                                "MDL_Fluid_Benchmark_Radial_Couette_Output.exo",
                                { "HMR_dummy_c_p1", "HMR_dummy_n_p1" },
                                { "VX", "VY", "P", "L2", "H1" },
                                { vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL, vis::Field_Type::NODAL },
                                { vis::Output_Type::VX,  vis::Output_Type::VY, vis::Output_Type::P, vis::Output_Type::L2_ERROR_ANALYTIC, vis::Output_Type::H1_ERROR_ANALYTIC } );
       tModel->set_output_manager( &tOutputData );

       // --------------------------------------------------------------------------------------
       // define linear solver and algorithm
       dla::Solver_Factory  tSolFactory;
       std::shared_ptr< dla::Linear_Solver_Algorithm > tLinearSolverAlgorithm
       = tSolFactory.create_solver( sol::SolverType::AMESOS_IMPL );

//       tLinearSolverAlgorithm->set_param("AZ_diagnostics") = AZ_all;
//       tLinearSolverAlgorithm->set_param("AZ_output") = AZ_all;
//       tLinearSolverAlgorithm->set_param("AZ_solver") = AZ_gmres_condnum;

       dla::Linear_Solver tLinSolver;

       tLinSolver.set_linear_algorithm( 0, tLinearSolverAlgorithm );

       // --------------------------------------------------------------------------------------
       // define nonlinear solver and algorithm
       NLA::Nonlinear_Solver_Factory tNonlinFactory;
       std::shared_ptr< NLA::Nonlinear_Algorithm > tNonlinearSolverAlgorithm
       = tNonlinFactory.create_nonlinear_solver( NLA::NonlinearSolverType::NEWTON_SOLVER );

       tNonlinearSolverAlgorithm->set_linear_solver( &tLinSolver );

       NLA::Nonlinear_Solver tNonlinearSolver;
       tNonlinearSolver.set_nonlinear_algorithm( tNonlinearSolverAlgorithm, 0 );

       // --------------------------------------------------------------------------------------
       // define time solver and algorithm
       tsa::Time_Solver_Factory tTimeSolverFactory;
       std::shared_ptr< tsa::Time_Solver_Algorithm > tTimeSolverAlgorithm
       = tTimeSolverFactory.create_time_solver( tsa::TimeSolverType::MONOLITHIC );

       tTimeSolverAlgorithm->set_nonlinear_solver( &tNonlinearSolver );

       tsa::Time_Solver tTimeSolver;

       tTimeSolver.set_time_solver_algorithm( tTimeSolverAlgorithm );

       sol::SOL_Warehouse tSolverWarehouse;

       tSolverWarehouse.set_solver_interface(tModel->get_solver_interface());

       tNonlinearSolver.set_solver_warehouse( &tSolverWarehouse );
       tTimeSolver.set_solver_warehouse( &tSolverWarehouse );

       tNonlinearSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );
       tTimeSolver.set_dof_type_list( { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::P } );

       tTimeSolver.set_output( 0, tSolverOutputCriteria_MDLFluidBench );

       // --------------------------------------------------------------------------------------
       // solve and check
       tTimeSolver.solve();

//         //print full solution for debug
//         Matrix<DDRMat> tFullSol;
//         tTimeSolver.get_full_solution(tFullSol);
//         print(tFullSol,"tFullSol");

        // clean up
        //------------------------------------------------------------------------------
        delete tIntegMesh1;
        delete tModel;
    }
}
} /* end_moris_namespace */
