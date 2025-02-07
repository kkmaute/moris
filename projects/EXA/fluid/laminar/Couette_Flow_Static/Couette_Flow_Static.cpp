/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Couette_Flow_Static.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    std::string sFluid      = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string sFluidGhost = "ghost_p2";
    std::string sIn         = "iside_b0_2_b1_0";
    std::string sOut        = "iside_b0_2_b1_3";

    std::string sInterfaces = sIn + "," + sOut;

    std::string sAllDofTypes      = "VX,VY;P";
    std::string sVelocityDofTypes = "VX,VY";
    std::string sPressureDofTypes = "P";

    // Geometry Parameters
    Matrix< DDRMat > tCenterPoint = { { 0.0, 0.0 } }; /* Center point of the circle */
    moris::real      tRIn         = 1.0;              /* Inner circle radius (m) */
    moris::real      tROut        = 2.0;              /* Outer circle radius (m) */

    // Velocity Parameters
    moris::real tOmegaIn  = -5.0;
    moris::real tOmegaOut = 5.0;
    moris::real tKn       = 0.0;
    moris::real tb        = 1.0 / tRIn + 2.0 * tKn / std::pow( tRIn, 2.0 );
    moris::real td        = 1.0 / tROut - 2.0 * tKn / std::pow( tROut, 2.0 );
    moris::real tB        = ( tOmegaOut * tROut - tOmegaIn * tRIn * tROut / tRIn ) / ( td - tb * tROut / tRIn );
    moris::real tA        = ( tOmegaIn * tRIn - tb * tB ) / tRIn;

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    ImposedVelocityFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // radius
        moris::real tx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        moris::real ty = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
        moris::real tR = std::pow( std::pow( tx, 2.0 ) + std::pow( ty, 2.0 ), 0.5 );

        // radial velocity
        moris::real tVTheta = tA * tR + tB / tR;

        // set size for aPropMatrix
        aPropMatrix.set_size( 2, 1, 0.0 );
        aPropMatrix( 0 ) = -tVTheta * ty / tR;
        aPropMatrix( 1 ) = tVTheta * tx / tR;
    }

    void
    AnalyticdVelocitydxFunc(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // radius
        real tx = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real ty = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );
        real tR = std::pow( std::pow( tx, 2.0 ) + std::pow( ty, 2.0 ), 0.5 );

        real tx2 = std::pow( tx, 2.0 );
        real ty2 = std::pow( ty, 2.0 );
        real tR2 = std::pow( tR, 2.0 );
        real tR4 = std::pow( tR, 4.0 );

        // set size for aPropMatrix
        aPropMatrix.set_size( 2, 2, 0.0 );
        aPropMatrix( 0, 0 ) = 2.0 * tB * tx * ty / tR4;
        aPropMatrix( 1, 0 ) = -tA - tB / tR2 + 2.0 * tB * ty2 / tR4;

        aPropMatrix( 0, 1 ) = tA + tB / tR2 - 2.0 * tB * tx2 / tR4;
        aPropMatrix( 1, 1 ) = -2.0 * tB * tx * ty / tR4;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 22, 22 );
        aParameterLists.set( "domain_dimensions", 5.0, 5.0 );
        aParameterLists.set( "domain_offset", "-2.5+0.25,-2.5+0.25" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 3 );
        aParameterLists.set( "initial_refinement", "0" );
        aParameterLists.set( "initial_refinement_pattern", "0" );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 0 );

        aParameterLists.set( "adaptive_refinement_level", 2 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "radius", 1.0 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        aParameterLists.set( "radius", 2.0 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // fill the property part of the parameter list

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropViscosity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletInU" );
        aParameterLists.set( "value_function", "ImposedVelocityFunc" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletOutU" );
        aParameterLists.set( "value_function", "ImposedVelocityFunc" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropdUdx" );
        aParameterLists.set( "value_function", "AnalyticdVelocitydxFunc" );

        //------------------------------------------------------------------------------
        // fill the constitutive model part of the parameter list

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMFluid" );
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::FLUID_INCOMPRESSIBLE ) ;
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );
        aParameterLists.set( "properties", "PropViscosity,Viscosity;PropDensity,Density" );

        //------------------------------------------------------------------------------
        // fill the stabilization parameter part of the parameter list

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPIncFlow" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::INCOMPRESSIBLE_FLOW ) ;
        aParameterLists.set( "function_parameters", "36.0" );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY;P", "Velocity,Pressure" ) );

        // create parameter list for stabilization parameter 2
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichletNitscheU" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters", "100.0/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );

        // create parameter list for stabilization parameter 3
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPViscosity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::VISCOUS_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity" );

        // create parameter list for stabilization parameter 4
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPVelocity" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::CONVECTIVE_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.05" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );

        // create parameter list for stabilization parameter 5
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGPPressure" );
        aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::PRESSURE_GHOST ) ;
        aParameterLists.set( "function_parameters", "0.005/1.0" );
        aParameterLists.set( "leader_dof_dependencies", std::pair< std::string, std::string >( "VX,VY", "Velocity" ) );
        aParameterLists.set( "leader_properties", "PropViscosity,Viscosity;PropDensity,Density" );

        //------------------------------------------------------------------------------
        // fill the IWG part of the parameter list

        // create parameter list for IWG 1
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGVelocityBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", sFluid );

        // create parameter list for IWG 2
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGPressureBulk" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPIncFlow,IncompressibleFlow" );
        aParameterLists.set( "mesh_set_names", sFluid );

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropDirichletInU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", sIn );

        // create parameter list for IWG 4
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropDirichletInU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", sIn );

        // create parameter list for IWG 5
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGOutVelocity" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropDirichletOutU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "stabilization_parameters", "SPDirichletNitscheU,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", sOut );

        // create parameter list for IWG 6
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGOutPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "leader_properties", "PropDirichletOutU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMFluid,IncompressibleFluid" );
        aParameterLists.set( "mesh_set_names", sOut );

        // create parameter list for IWG 7
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPViscous" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "stabilization_parameters", "SPGPViscosity,GhostSP" );
        aParameterLists.set( "mesh_set_names", sFluidGhost );

        // create parameter list for IWG 8
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPConvective" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "stabilization_parameters", "SPGPVelocity,GhostSP" );
        aParameterLists.set( "mesh_set_names", sFluidGhost );

        // create parameter list for IWG 9
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGPPressure" );
        aParameterLists.set( "IWG_type",  fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
        aParameterLists.set( "dof_residual", "P" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "follower_dof_dependencies", "VX,VY;P" );
        aParameterLists.set( "stabilization_parameters", "SPGPPressure,GhostSP" );
        aParameterLists.set( "mesh_set_names", sFluidGhost );

        //------------------------------------------------------------------------------
        // fill the IQI part of the parameter list

        // create parameter list for IQI 1
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVX" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", sFluid );

        // create parameter list for IQI 2
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVY" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", sFluid );

        // create parameter list for IQI 3
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkP" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity", "P" );
        aParameterLists.set( "leader_dof_dependencies", "P" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", sFluid );

        // create parameter list for IQI 4
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::L2_ERROR_ANALYTIC ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropDirichletInU,L2Check" );
        aParameterLists.set( "mesh_set_names", sFluid );

        // create parameter list for IQI 5
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkH1Error" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::H1_ERROR_ANALYTIC ) ;
        aParameterLists.set( "dof_quantity", "VX,VY" );
        aParameterLists.set( "leader_dof_dependencies", "VX,VY" );
        aParameterLists.set( "leader_properties", "PropdUdx,H1Check" );
        aParameterLists.set( "mesh_set_names", sFluid );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", 1e-06 );
        aParameterLists.set( "NLA_relaxation_parameter", 1.0 );
        aParameterLists.set( "NLA_max_iter", 20 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "VX,VY;P" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "VX,VY;P" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "VX,1E-4;VY,1E-4;P,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Couette_Flow_Static.exo" ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", sFluid + "," + sInterfaces );
        aParameterLists.set( "Field_Names", "VX,VY,P,L2_ERROR_ANALYTIC,H1_ERROR_ANALYTIC" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkVX,IQIBulkVY,IQIBulkP,IQIBulkL2Error,IQIBulkH1Error" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif

